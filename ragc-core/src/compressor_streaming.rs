// Streaming AGC Compressor
// Memory-efficient implementation that flushes groups incrementally

use crate::{
    genome_io::{GenomeIO, parse_sample_from_header},
    kmer::{Kmer, KmerMode},
    lz_diff::LZDiff,
    segment::split_at_splitters_with_size,
    segment_compression::compress_segment,
    splitters::determine_splitters,
};
use anyhow::{Context, Result};
use ragc_common::{
    stream_delta_name, stream_ref_name, Archive, CollectionV3, Contig, AGC_FILE_MAJOR, AGC_FILE_MINOR,
    CONTIG_SEPARATOR,
};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

/// Configuration for the streaming compressor
#[derive(Debug, Clone)]
pub struct StreamingCompressorConfig {
    pub kmer_length: u32,
    pub segment_size: u32,
    pub min_match_len: u32,
    pub verbosity: u32,
    /// Flush groups to disk when they reach this size (0 = only flush at finalize)
    pub group_flush_threshold: usize,
    /// Flush all groups after processing this many contigs (0 = never auto-flush)
    pub periodic_flush_interval: usize,
}

impl Default for StreamingCompressorConfig {
    fn default() -> Self {
        StreamingCompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            verbosity: 1,
            group_flush_threshold: 100, // Flush after 100 segments per group
            periodic_flush_interval: 500, // Flush all groups after every 500 contigs
        }
    }
}

/// A segment with its metadata
#[derive(Debug, Clone)]
struct SegmentInfo {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    data: Contig,
    is_rev_comp: bool,
}

/// Segment group identified by flanking k-mers
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct SegmentGroupKey {
    kmer_front: u64,
    kmer_back: u64,
}

/// Group metadata for tracking what's been written
struct GroupMetadata {
    group_id: u32,
    stream_id: usize,              // Delta stream ID
    ref_stream_id: Option<usize>,  // Reference stream ID (for LZ groups)
    reference: Option<Contig>,     // First segment (for LZ encoding)
    ref_written: bool,             // Whether reference has been written to archive
    segments_written: usize,       // Number of segments written to archive (not including buffered)
    pending_segments: Vec<SegmentInfo>, // Buffered segments waiting for complete pack
    is_flushed: bool,
}

/// Streaming AGC Compressor
pub struct StreamingCompressor {
    config: StreamingCompressorConfig,
    archive: Archive,
    collection: CollectionV3,

    // Segment storage (temporary, flushed periodically)
    segment_groups: HashMap<SegmentGroupKey, Vec<SegmentInfo>>,

    // Group tracking
    group_metadata: HashMap<SegmentGroupKey, GroupMetadata>,
    next_group_id: u32,

    // Statistics
    total_bases_processed: usize,
    total_segments: usize,
    total_groups_flushed: usize,
    contigs_since_flush: usize,
}

impl StreamingCompressor {
    /// Create a new streaming compressor
    pub fn new(archive_path: &str, config: StreamingCompressorConfig) -> Result<Self> {
        let mut archive = Archive::new_writer();
        archive
            .open(archive_path)
            .context("Failed to create archive")?;

        let mut collection = CollectionV3::new();
        collection.set_config(config.segment_size, config.kmer_length, None);
        collection.prepare_for_compression(&mut archive)?;

        Ok(StreamingCompressor {
            config,
            archive,
            collection,
            segment_groups: HashMap::new(),
            group_metadata: HashMap::new(),
            next_group_id: 0,
            total_bases_processed: 0,
            total_segments: 0,
            total_groups_flushed: 0,
            contigs_since_flush: 0,
        })
    }

    /// Add a FASTA file to the archive (streaming mode)
    pub fn add_fasta_file(&mut self, sample_name: &str, fasta_path: &Path) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Processing sample: {sample_name} from {fasta_path:?}");
        }

        let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
            .context("Failed to open FASTA file")?;

        let mut contigs_processed = 0;

        // Read contigs with conversion (ASCII -> numeric)
        while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
            self.add_contig(sample_name, &contig_name, sequence)?;
            contigs_processed += 1;

            // Periodic flush of all groups
            if self.config.periodic_flush_interval > 0 &&
               self.contigs_since_flush >= self.config.periodic_flush_interval {
                self.flush_all_groups()?;
            }

            if self.config.verbosity > 1 && contigs_processed % 100 == 0 {
                println!("  Processed {contigs_processed} contigs, {} groups in memory, {} flushed",
                    self.segment_groups.len(), self.total_groups_flushed);
            }
        }

        if self.config.verbosity > 0 {
            println!("  Processed {contigs_processed} contigs from {sample_name}");
        }

        Ok(())
    }

    /// Detect if a FASTA file contains multiple samples in headers
    /// (format: >sample#haplotype#chromosome)
    pub fn detect_multi_sample_fasta(fasta_path: &Path) -> Result<bool> {
        let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
            .context("Failed to open FASTA for detection")?;

        // Read first few headers to check format
        for _ in 0..5 {
            if let Some((header, _, _, _)) = reader.read_contig_with_sample()? {
                let parts: Vec<&str> = header.split('#').collect();
                if parts.len() >= 3 {
                    return Ok(true); // Multi-sample format detected
                }
            } else {
                break;
            }
        }

        Ok(false) // No multi-sample format detected
    }

    /// Add a multi-sample FASTA file with splitter-based segmentation
    /// This handles files where samples are encoded in headers (>sample#hap#chr format)
    pub fn add_multi_sample_fasta_with_splitters(&mut self, fasta_path: &Path) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("=== Processing multi-sample FASTA (grouping by sample names in headers) ===");
            println!("Input: {:?}", fasta_path);
            println!();
        }

        // Pass 1: Collect reference contigs for splitter finding
        if self.config.verbosity > 0 {
            println!("=== Pass 1: Finding splitter k-mers from reference ===");
        }

        let mut reference_contigs = Vec::new();
        let mut reference_sample = String::new();

        {
            let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                .context("Failed to open multi-sample FASTA")?;

            let mut contig_count = 0;
            while let Some((_full_header, sample_name, _contig_name, sequence)) = reader.read_contig_with_sample()? {
                if !sequence.is_empty() {
                    if reference_sample.is_empty() {
                        reference_sample = sample_name.clone();
                        if self.config.verbosity > 0 {
                            println!("Using {} as reference to find splitters...", reference_sample);
                        }
                    }

                    // Only collect contigs from the first sample
                    if sample_name == reference_sample {
                        reference_contigs.push(sequence);
                        contig_count += 1;
                    }
                }
            }

            if self.config.verbosity > 0 {
                println!("Collected {} reference contigs from {}", contig_count, reference_sample);
            }
        }

        // Find splitters
        let splitters = determine_splitters(
            &reference_contigs,
            self.config.kmer_length as usize,
            self.config.segment_size as usize
        );

        if self.config.verbosity > 0 {
            println!("Found {} actually-used splitter k-mers", splitters.len());
            println!();
            println!("=== Pass 2: Segmenting and compressing all samples ===");
        }

        drop(reference_contigs); // Free memory

        // Pass 2: Re-read file and process all samples
        let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
            .context("Failed to re-open multi-sample FASTA")?;

        let mut samples_seen = std::collections::HashSet::new();
        let mut contigs_processed = 0;

        while let Some((_full_header, sample_name, contig_name, sequence)) = reader.read_contig_with_sample()? {
            if sequence.is_empty() {
                continue;
            }

            // Track unique samples
            if samples_seen.insert(sample_name.clone()) && self.config.verbosity > 0 {
                println!("Processing sample: {}", sample_name);
            }

            // Register in collection
            self.collection.register_sample_contig(&sample_name, &contig_name)?;

            // Split at splitters
            let segments = split_at_splitters_with_size(
                &sequence,
                &splitters,
                self.config.kmer_length as usize,
                self.config.segment_size as usize
            );

            if self.config.verbosity > 1 {
                println!("  Contig {}/{} split into {} segments", sample_name, contig_name, segments.len());
            }

            // Add each segment
            for (seg_idx, segment) in segments.iter().enumerate() {
                self.total_bases_processed += segment.data.len();

                if self.config.verbosity > 2 {
                    println!("    Segment {}: front_kmer={}, back_kmer={}, len={}",
                        seg_idx, segment.front_kmer, segment.back_kmer, segment.data.len());
                }

                self.add_segment_with_kmers(
                    &sample_name,
                    &contig_name,
                    seg_idx,
                    segment.data.clone(),
                    segment.front_kmer,
                    segment.back_kmer,
                    false, // is_rev_comp
                )?;
            }

            contigs_processed += 1;
            self.contigs_since_flush += 1;

            // Periodic flush
            if self.config.periodic_flush_interval > 0 &&
               self.contigs_since_flush >= self.config.periodic_flush_interval {
                self.flush_all_groups()?;
            }

            if self.config.verbosity > 1 && contigs_processed % 100 == 0 {
                println!("  Processed {} contigs, {} groups in memory, {} flushed",
                    contigs_processed, self.segment_groups.len(), self.total_groups_flushed);
            }
        }

        if self.config.verbosity > 0 {
            println!("Processed {} unique samples, {} total contigs", samples_seen.len(), contigs_processed);
        }

        Ok(())
    }

    /// Add multiple FASTA files with splitter-based segmentation (two-pass)
    ///
    /// This is the real AGC algorithm:
    /// Pass 1: Read all genomes → find singleton k-mers (splitters)
    /// Pass 2: Re-read genomes → segment at splitters → compress
    ///
    /// This achieves much better compression than treating each contig as one segment!
    pub fn add_fasta_files_with_splitters(&mut self, fasta_paths: &[(String, &Path)]) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("=== Pass 1: Finding splitter k-mers across all genomes ===");
        }

        // Pass 1: Collect sequences from FIRST genome only to find splitters
        // This is the reference-based approach: use splitters from first genome
        // to segment all genomes at the same positions
        let mut reference_contigs = Vec::new();

        if let Some((ref_sample_name, ref_fasta_path)) = fasta_paths.first() {
            if self.config.verbosity > 0 {
                println!("Using {} as reference to find splitters...", ref_sample_name);
            }

            let mut reader = GenomeIO::<Box<dyn Read>>::open(ref_fasta_path)
                .context(format!("Failed to open reference {}", ref_sample_name))?;

            while let Some((_contig_name, sequence)) = reader.read_contig_converted()? {
                if !sequence.is_empty() {
                    reference_contigs.push(sequence);
                }
            }

            if self.config.verbosity > 0 {
                println!("Collected {} reference contigs, finding singleton k-mers...", reference_contigs.len());
            }
        } else {
            anyhow::bail!("No input files provided");
        }

        // Find splitters (three-pass algorithm):
        // 1. Find singletons (candidates)
        // 2. Scan reference to find which candidates are actually used
        // 3. Return only actually-used splitters
        let splitters = determine_splitters(
            &reference_contigs,
            self.config.kmer_length as usize,
            self.config.segment_size as usize
        );

        if self.config.verbosity > 0 {
            println!("Found {} actually-used splitter k-mers", splitters.len());
            println!();
            println!("=== Pass 2: Segmenting and compressing genomes ===");
        }

        // Free memory from reference_contigs
        drop(reference_contigs);

        // Pass 2: Re-read files and segment using splitters
        for (sample_name, fasta_path) in fasta_paths {
            if self.config.verbosity > 0 {
                println!("Processing sample: {} from {:?}", sample_name, fasta_path);
            }

            let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                .context(format!("Failed to re-open {}", sample_name))?;

            let mut contigs_processed = 0;

            while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
                if sequence.is_empty() {
                    continue;
                }

                // Register in collection
                self.collection.register_sample_contig(sample_name, &contig_name)?;

                // Split at splitter positions with minimum segment size
                let segments = split_at_splitters_with_size(
                    &sequence,
                    &splitters,
                    self.config.kmer_length as usize,
                    self.config.segment_size as usize
                );

                if self.config.verbosity > 1 {
                    println!("  Contig {} split into {} segments", contig_name, segments.len());
                }

                // Add each segment
                for (seg_idx, segment) in segments.iter().enumerate() {
                    self.total_bases_processed += segment.data.len();

                    // CRITICAL: Check for segments < k (will cause C++ AGC errors)
                    if seg_idx > 0 && segment.data.len() < self.config.kmer_length as usize {
                        eprintln!("WARNING: Segment {} of contig {} has size {} < k={} bytes!",
                            seg_idx, contig_name, segment.data.len(), self.config.kmer_length);
                        eprintln!("  Sample: {}, front_kmer={}, back_kmer={}",
                            sample_name, segment.front_kmer, segment.back_kmer);
                        eprintln!("  This will cause 'Corrupted archive!' errors in C++ AGC");
                    }

                    if self.config.verbosity > 2 {
                        println!("    Contig {}, Segment {}: front_kmer={}, back_kmer={}, len={}",
                            contig_name, seg_idx, segment.front_kmer, segment.back_kmer, segment.data.len());
                    }

                    self.add_segment_with_kmers(
                        sample_name,
                        &contig_name,
                        seg_idx,
                        segment.data.clone(),
                        segment.front_kmer,
                        segment.back_kmer,
                        false, // is_rev_comp
                    )?;
                }

                contigs_processed += 1;
                self.contigs_since_flush += 1;

                // Periodic flush
                if self.config.periodic_flush_interval > 0 &&
                   self.contigs_since_flush >= self.config.periodic_flush_interval {
                    self.flush_all_groups()?;
                }

                if self.config.verbosity > 1 && contigs_processed % 100 == 0 {
                    println!("  Processed {} contigs, {} groups in memory, {} flushed",
                        contigs_processed, self.segment_groups.len(), self.total_groups_flushed);
                }
            }

            if self.config.verbosity > 0 {
                println!("  Processed {} contigs from {}", contigs_processed, sample_name);
            }
        }

        Ok(())
    }

    /// Add a contig to the archive
    pub fn add_contig(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        sequence: Contig,
    ) -> Result<()> {
        if sequence.is_empty() {
            return Ok(());
        }

        self.total_bases_processed += sequence.len();
        self.contigs_since_flush += 1;

        // Register in collection
        self.collection
            .register_sample_contig(sample_name, contig_name)?;

        // For now, treat entire contig as a single segment
        self.add_segment(sample_name, contig_name, 0, sequence)?;

        Ok(())
    }

    /// Add a segment to the compressor
    #[allow(clippy::needless_range_loop)]
    fn add_segment(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_part_no: usize,
        segment: Contig,
    ) -> Result<()> {
        // Extract flanking k-mers
        let k = self.config.kmer_length;

        let (kmer_front, kmer_back) = if segment.len() >= (k * 2) as usize {
            // Extract front k-mer
            let mut front = Kmer::new(k, KmerMode::Canonical);
            for i in 0..(k as usize) {
                if segment[i] > 3 {
                    front.reset();
                    break;
                }
                front.insert(segment[i] as u64);
            }
            let front_kmer_val = if front.is_full() { front.data() } else { 0 };

            // Extract back k-mer
            let mut back = Kmer::new(k, KmerMode::Canonical);
            let start = segment.len() - (k as usize);
            for i in 0..(k as usize) {
                if segment[start + i] > 3 {
                    back.reset();
                    break;
                }
                back.insert(segment[start + i] as u64);
            }
            let back_kmer_val = if back.is_full() { back.data() } else { 0 };

            (front_kmer_val, back_kmer_val)
        } else {
            (0u64, 0u64)
        };

        self.add_segment_with_kmers(
            sample_name,
            contig_name,
            seg_part_no,
            segment,
            kmer_front,
            kmer_back,
            false, // is_rev_comp
        )
    }

    /// Add a segment with known k-mers
    #[allow(clippy::too_many_arguments)]
    fn add_segment_with_kmers(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_part_no: usize,
        segment: Contig,
        kmer_front: u64,
        kmer_back: u64,
        is_rev_comp: bool,
    ) -> Result<()> {
        let key = SegmentGroupKey {
            kmer_front,
            kmer_back,
        };

        let seg_info = SegmentInfo {
            sample_name: sample_name.to_string(),
            contig_name: contig_name.to_string(),
            seg_part_no,
            data: segment,
            is_rev_comp,
        };

        // Add to group
        self.segment_groups.entry(key.clone()).or_default().push(seg_info);
        self.total_segments += 1;

        // Check if this group should be flushed
        if self.config.group_flush_threshold > 0 {
            let group_size = self.segment_groups.get(&key).map(|v| v.len()).unwrap_or(0);
            if group_size >= self.config.group_flush_threshold {
                self.flush_group(&key)?;
            }
        }

        Ok(())
    }

    /// Flush a single group to the archive
    fn flush_group(&mut self, key: &SegmentGroupKey) -> Result<()> {
        let segments = match self.segment_groups.remove(key) {
            Some(segs) if !segs.is_empty() => segs,
            _ => return Ok(()), // Nothing to flush
        };

        if self.config.verbosity > 1 {
            println!("  Flushing group with {} segments", segments.len());
        }

        // Get or create group metadata
        let metadata = self.group_metadata.entry(key.clone()).or_insert_with(|| {
            let group_id = self.next_group_id;
            self.next_group_id += 1;

            if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                eprintln!("DEBUG: Creating group {} for key front_kmer={}, back_kmer={}",
                    group_id, key.kmer_front, key.kmer_back);
            }

            let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;

            // Register delta stream (always)
            let stream_name = stream_delta_name(archive_version, group_id);
            let stream_id = self.archive.register_stream(&stream_name);

            // Register ref stream for LZ groups (groups >= 16)
            const NO_RAW_GROUPS: u32 = 16;
            let ref_stream_id = if group_id >= NO_RAW_GROUPS {
                let ref_stream_name = stream_ref_name(archive_version, group_id);
                Some(self.archive.register_stream(&ref_stream_name))
            } else {
                None
            };

            GroupMetadata {
                group_id,
                stream_id,
                ref_stream_id,
                reference: None,
                ref_written: false,
                segments_written: 0,
                pending_segments: Vec::new(),
                is_flushed: false,
            }
        });

        let group_id = metadata.group_id;
        let stream_id = metadata.stream_id;

        // Use LZ encoding for groups >= 16 (matches C++ AGC default)
        // Groups 0-15 are raw-only for better random access
        // Groups 16+ use LZ encoding with reference for better compression
        const NO_RAW_GROUPS: u32 = 16;
        let use_lz_encoding = group_id >= NO_RAW_GROUPS;

        // Remember if we need to skip first segment (before we write ref and set ref_written=true)
        let skip_first_segment = use_lz_encoding && !metadata.ref_written;

        // For LZ groups: write reference segment to separate ref stream
        if skip_first_segment && !segments.is_empty() {
            let ref_segment = &segments[0];

            // Store as reference for LZ encoding
            metadata.reference = Some(ref_segment.data.clone());

            // Write to ref stream
            if let Some(ref_stream_id) = metadata.ref_stream_id {
                // Try compressing the reference
                let compressed = compress_segment(&ref_segment.data)?;

                if compressed.len() + 1 < ref_segment.data.len() {
                    let mut compressed_with_marker = compressed;
                    compressed_with_marker.push(0); // Marker byte: 0 = standard ZSTD
                    self.archive.add_part(ref_stream_id, &compressed_with_marker, ref_segment.data.len() as u64)?;
                } else {
                    // Uncompressed
                    self.archive.add_part(ref_stream_id, &ref_segment.data, 0)?;
                }

                // Register reference in collection (in_group_id = 0)
                self.collection.add_segment_placed(
                    &ref_segment.sample_name,
                    &ref_segment.contig_name,
                    ref_segment.seg_part_no,
                    group_id,
                    0, // Reference is always at position 0
                    ref_segment.is_rev_comp,
                    ref_segment.data.len() as u32,
                )?;

                metadata.ref_written = true;

                if self.config.verbosity > 2 {
                    eprintln!("DEBUG: Group {} wrote reference segment to ref stream, len={}",
                        group_id, ref_segment.data.len());
                }
            }
        }

        const PACK_CARDINALITY: usize = 50;

        // Determine which segments to pack into delta stream
        // IMPORTANT: Only skip first segment if this is the first flush for an LZ group!
        let mut new_delta_segments: Vec<SegmentInfo> = if skip_first_segment {
            // Skip first segment (just written to ref stream above)
            segments.into_iter().skip(1).collect()
        } else {
            // Include all segments for raw-only groups OR if ref already written
            segments
        };

        // Add new segments to pending buffer
        metadata.pending_segments.append(&mut new_delta_segments);

        // Only write complete packs (multiples of PACK_CARDINALITY)
        let num_complete_packs = metadata.pending_segments.len() / PACK_CARDINALITY;
        let segments_to_write = num_complete_packs * PACK_CARDINALITY;

        if segments_to_write == 0 {
            // Not enough for a full pack yet, keep buffering
            if self.config.verbosity > 2 {
                eprintln!("Group {} buffering {} segments (need {} for complete pack)",
                    group_id, metadata.pending_segments.len(), PACK_CARDINALITY);
            }
            return Ok(());
        }

        // Split off segments to write now
        let segments_to_pack: Vec<SegmentInfo> = metadata.pending_segments.drain(..segments_to_write).collect();

        if self.config.verbosity > 2 {
            eprintln!("Group {} writing {} complete packs ({} segments), {} segments remain buffered",
                group_id, num_complete_packs, segments_to_write, metadata.pending_segments.len());
        }

        // Pack and write delta segments in complete packs
        let mut pack_count = 0;
        for pack_segments in segments_to_pack.chunks(PACK_CARDINALITY) {
            let mut packed_data = Vec::new();

            for (idx_in_pack, seg_info) in pack_segments.iter().enumerate() {
                // For LZ groups: in_group_id starts at 1 (reference is 0)
                // For raw groups: in_group_id starts at 0
                let global_in_group_id = if use_lz_encoding {
                    metadata.segments_written + idx_in_pack + 1
                } else {
                    metadata.segments_written + idx_in_pack
                };

                if (group_id < 50 || group_id == 93) && idx_in_pack < 3 {
                    eprintln!("REGISTER_NORMAL: group={}, in_group_id={}, idx_in_pack={}, segments_written={}",
                        group_id, global_in_group_id, idx_in_pack, metadata.segments_written);
                }

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} adding segment: contig={}, seg_part={}, in_group_id={}, data_len={}",
                        group_id, seg_info.contig_name, seg_info.seg_part_no, global_in_group_id, seg_info.data.len());
                }

                let contig_data = if use_lz_encoding {
                    // LZ-encoded segment (all delta segments are encoded)
                    let mut lz_diff = LZDiff::new(self.config.min_match_len);
                    if let Some(ref reference) = metadata.reference {
                        lz_diff.prepare(reference);
                        lz_diff.encode(&seg_info.data)
                    } else {
                        seg_info.data.clone()
                    }
                } else {
                    // Raw segment
                    seg_info.data.clone()
                };

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} contig_data len after encoding: {}", group_id, contig_data.len());
                }

                packed_data.extend_from_slice(&contig_data);
                packed_data.push(CONTIG_SEPARATOR);

                // Register in collection
                if self.config.verbosity > 2 && group_id == 801 {
                    eprintln!("REGISTER: group={}, in_group_id={}, pack_id would be {}, contig={}, seg_part={}",
                        group_id, global_in_group_id, (global_in_group_id - 1) / 50, seg_info.contig_name, seg_info.seg_part_no);
                }
                self.collection.add_segment_placed(
                    &seg_info.sample_name,
                    &seg_info.contig_name,
                    seg_info.seg_part_no,
                    group_id,
                    global_in_group_id as u32,
                    seg_info.is_rev_comp,
                    seg_info.data.len() as u32,
                )?;
            }

            // Try compressing
            let compressed = compress_segment(&packed_data)?;

            // C++ AGC logic: only use compression if it's beneficial
            // If compressed + marker >= uncompressed, write uncompressed instead
            if compressed.len() + 1 < packed_data.len() {
                // Compression is beneficial
                // Marker byte at the END (C++ AGC format)
                let mut compressed_with_marker = compressed;
                compressed_with_marker.push(0); // Marker byte: 0 = standard ZSTD

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} pack: packed_len={}, compressed_len={}, writing compressed",
                        group_id, packed_data.len(), compressed_with_marker.len());
                }

                self.archive.add_part(stream_id, &compressed_with_marker, packed_data.len() as u64)?;
            } else {
                // Compression not beneficial, write uncompressed
                // CRITICAL: Do NOT add marker byte for uncompressed data!
                // C++ AGC expects raw data when uncompressed_size = 0

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} pack: packed_len={}, writing uncompressed (no marker)",
                        group_id, packed_data.len());
                }

                self.archive.add_part(stream_id, &packed_data, 0)?;
            }

            pack_count += 1;
            metadata.segments_written += pack_segments.len();
        }

        if group_id < 50 || group_id == 93 {
            eprintln!("FLUSH_NORMAL: group={} wrote {} packs, segments_written now={}",
                group_id, pack_count, metadata.segments_written);
        }

        metadata.is_flushed = true;
        self.total_groups_flushed += 1;

        Ok(())
    }

    /// Flush all currently buffered groups to disk
    fn flush_all_groups(&mut self) -> Result<()> {
        if self.config.verbosity > 1 {
            println!("  Periodic flush: flushing {} groups", self.segment_groups.len());
        }

        let keys: Vec<_> = self.segment_groups.keys().cloned().collect();
        for key in keys {
            self.flush_group(&key)?;
        }

        self.contigs_since_flush = 0;
        Ok(())
    }

    /// Flush any remaining buffered segments for a group (called at finalize)
    /// This writes partial packs for groups with buffered segments
    fn flush_group_final(&mut self, key: &SegmentGroupKey) -> Result<()> {
        const PACK_CARDINALITY: usize = 50;

        // Get metadata and check if there are pending segments
        let metadata = match self.group_metadata.get_mut(key) {
            Some(m) if !m.pending_segments.is_empty() => m,
            _ => return Ok(()), // No pending segments
        };

        let group_id = metadata.group_id;
        let stream_id = metadata.stream_id;
        let use_lz_encoding = group_id >= 16; // NO_RAW_GROUPS

        eprintln!("Final flush for group {}: {} pending segments, segments_written={}",
            group_id, metadata.pending_segments.len(), metadata.segments_written);

        // Take all pending segments (even if less than PACK_CARDINALITY)
        let segments_to_pack: Vec<SegmentInfo> = metadata.pending_segments.drain(..).collect();

        eprintln!("Group {} will write {} packs in final flush",
            group_id, (segments_to_pack.len() + PACK_CARDINALITY - 1) / PACK_CARDINALITY);

        // Pack and write segments (may be partial pack)
        // Track the base in_group_id for this flush operation
        let mut current_segments_written = metadata.segments_written;

        for pack_segments in segments_to_pack.chunks(PACK_CARDINALITY) {
            let mut packed_data = Vec::new();

            for (idx_in_pack, seg_info) in pack_segments.iter().enumerate() {
                let global_in_group_id = if use_lz_encoding {
                    current_segments_written + idx_in_pack + 1
                } else {
                    current_segments_written + idx_in_pack
                };

                let contig_data = if use_lz_encoding {
                    let mut lz_diff = LZDiff::new(self.config.min_match_len);
                    if let Some(ref reference) = metadata.reference {
                        lz_diff.prepare(reference);
                        lz_diff.encode(&seg_info.data)
                    } else {
                        seg_info.data.clone()
                    }
                } else {
                    seg_info.data.clone()
                };

                packed_data.extend_from_slice(&contig_data);
                packed_data.push(CONTIG_SEPARATOR);

                if group_id < 50 || group_id == 93 {
                    eprintln!("REGISTER_FINAL: group={}, in_group_id={}, contig={}, seg_part={}, idx_in_pack={}, current_segments_written={}",
                        group_id, global_in_group_id, seg_info.contig_name, seg_info.seg_part_no, idx_in_pack, current_segments_written);
                }

                self.collection.add_segment_placed(
                    &seg_info.sample_name,
                    &seg_info.contig_name,
                    seg_info.seg_part_no,
                    group_id,
                    global_in_group_id as u32,
                    seg_info.is_rev_comp,
                    seg_info.data.len() as u32,
                )?;
            }

            // Compress and write
            let compressed = compress_segment(&packed_data)?;

            if compressed.len() + 1 < packed_data.len() {
                let mut compressed_with_marker = compressed;
                compressed_with_marker.push(0);
                self.archive.add_part(stream_id, &compressed_with_marker, packed_data.len() as u64)?;
            } else {
                self.archive.add_part(stream_id, &packed_data, 0)?;
            }

            eprintln!("Group {} wrote pack to archive (final flush), pack has {} segments", group_id, pack_segments.len());
            current_segments_written += pack_segments.len();
        }

        eprintln!("Group {} final flush complete: wrote {} total segments in this flush, new segments_written={}",
            group_id, current_segments_written - metadata.segments_written, current_segments_written);

        // Update the metadata's segments_written counter
        metadata.segments_written = current_segments_written;

        Ok(())
    }

    /// Store params stream (C++ compatibility)
    fn store_params_stream(&mut self) -> Result<()> {
        let mut params_data = Vec::new();
        let append_u32 = |data: &mut Vec<u8>, value: u32| {
            data.extend_from_slice(&value.to_le_bytes());
        };

        append_u32(&mut params_data, self.config.kmer_length);
        append_u32(&mut params_data, self.config.min_match_len);
        append_u32(&mut params_data, 50); // pack_cardinality
        append_u32(&mut params_data, self.config.segment_size);
        append_u32(&mut params_data, 16); // no_raw_groups (groups 0-15 are raw-only, 16+ use LZ)

        let stream_id = self.archive.register_stream("params");
        self.archive
            .add_part(stream_id, &params_data, params_data.len() as u64)?;

        Ok(())
    }

    /// Store empty splitters stream (C++ compatibility)
    fn store_splitters_stream(&mut self) -> Result<()> {
        let splitters_data = Vec::new();
        let stream_id = self.archive.register_stream("splitters");
        self.archive.add_part(stream_id, &splitters_data, 0)?;
        Ok(())
    }

    /// Store empty segment-splitters stream (C++ compatibility)
    fn store_segment_splitters_stream(&mut self) -> Result<()> {
        let seg_splitters_data = Vec::new();
        let stream_id = self.archive.register_stream("segment-splitters");
        self.archive.add_part(stream_id, &seg_splitters_data, 0)?;
        Ok(())
    }

    /// Store file_type_info stream (C++ compatibility)
    fn store_file_type_info(&mut self) -> Result<()> {
        let mut data = Vec::new();
        let append_str = |data: &mut Vec<u8>, s: &str| {
            data.extend_from_slice(s.as_bytes());
            data.push(0);
        };

        append_str(&mut data, "producer");
        append_str(&mut data, "ragc-streaming");

        append_str(&mut data, "producer_version_major");
        append_str(&mut data, &AGC_FILE_MAJOR.to_string());

        append_str(&mut data, "producer_version_minor");
        append_str(&mut data, &AGC_FILE_MINOR.to_string());

        append_str(&mut data, "producer_version_build");
        append_str(&mut data, "0");

        append_str(&mut data, "file_version_major");
        append_str(&mut data, &AGC_FILE_MAJOR.to_string());

        append_str(&mut data, "file_version_minor");
        append_str(&mut data, &AGC_FILE_MINOR.to_string());

        append_str(&mut data, "comment");
        append_str(
            &mut data,
            &format!("AGC (Rust streaming) v.{AGC_FILE_MAJOR}.{AGC_FILE_MINOR}"),
        );

        let stream_id = self.archive.register_stream("file_type_info");
        self.archive.add_part(stream_id, &data, 7)?;

        Ok(())
    }

    /// Finalize compression and write any remaining segments
    pub fn finalize(&mut self) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Finalizing compression...");
            println!("Total segments: {}", self.total_segments);
            println!("Groups flushed: {}", self.total_groups_flushed);
            println!("Remaining groups to flush: {}", self.segment_groups.len());

            // Count segments per group
            let mut group_sizes: Vec<usize> = self.segment_groups.values()
                .map(|v| v.len())
                .collect();
            group_sizes.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending

            let total_groups = self.segment_groups.len() + self.total_groups_flushed;
            println!("\nGroup statistics:");
            println!("  Total unique groups: {}", total_groups);
            println!("  Groups with 1 segment: {}", group_sizes.iter().filter(|&&s| s == 1).count());
            println!("  Groups with 2+ segments: {}", group_sizes.iter().filter(|&&s| s > 1).count());
            if !group_sizes.is_empty() {
                println!("  Largest group: {} segments", group_sizes[0]);
                println!("  Top 10 group sizes: {:?}", &group_sizes[..group_sizes.len().min(10)]);
            }
        }

        // Flush any remaining groups
        let remaining_keys: Vec<_> = self.segment_groups.keys().cloned().collect();
        if self.config.verbosity > 1 {
            println!("Flushing {} remaining groups", remaining_keys.len());
        }

        let mut groups_actually_flushed = 0;
        for (idx, key) in remaining_keys.iter().enumerate() {
            if self.config.verbosity > 1 && idx % 1000 == 0 {
                println!("  Flushing group {}/{}", idx, remaining_keys.len());
            }
            self.flush_group(key)?;
            groups_actually_flushed += 1;
        }

        if self.config.verbosity > 0 {
            println!("Actually flushed {} groups", groups_actually_flushed);
        }

        // Flush any remaining buffered segments (partial packs)
        let all_group_keys: Vec<_> = self.group_metadata.keys().cloned().collect();
        if self.config.verbosity > 1 {
            println!("Flushing buffered segments for {} groups", all_group_keys.len());
        }
        for key in all_group_keys {
            self.flush_group_final(&key)?;
        }

        // Store metadata streams
        self.store_params_stream()?;
        self.store_splitters_stream()?;
        self.store_segment_splitters_stream()?;

        // Store collection metadata
        self.collection
            .store_batch_sample_names(&mut self.archive)?;

        let num_samples = self.collection.get_no_samples();
        if num_samples > 0 {
            self.collection
                .store_contig_batch(&mut self.archive, 0, num_samples)?;
        }

        // Store file_type_info stream (must be last)
        self.store_file_type_info()?;

        // Close archive
        self.archive.close()?;

        if self.config.verbosity > 0 {
            println!("Compression complete!");
            println!("Total bases: {}", self.total_bases_processed);
            println!("Total groups: {}", self.group_metadata.len());
        }

        Ok(())
    }
}
