// Streaming AGC Compressor
// Memory-efficient implementation that flushes groups incrementally

use crate::{
    genome_io::GenomeIO,
    kmer::{Kmer, KmerMode},
    lz_diff::LZDiff,
    segment::split_at_splitters_with_size,
    segment_compression::compress_segment,
    splitters::determine_splitters,
};
use anyhow::{Context, Result};
use ragc_common::{
    stream_delta_name, Archive, CollectionV3, Contig, AGC_FILE_MAJOR, AGC_FILE_MINOR,
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
    stream_id: usize,
    reference: Option<Contig>, // First segment (for LZ encoding)
    segments_written: usize,
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
            let stream_name = stream_delta_name(archive_version, group_id);
            let stream_id = self.archive.register_stream(&stream_name);

            GroupMetadata {
                group_id,
                stream_id,
                reference: None,
                segments_written: 0,
                is_flushed: false,
            }
        });

        let group_id = metadata.group_id;
        let stream_id = metadata.stream_id;
        // CRITICAL: Don't use LZ encoding since we don't write no_raw_groups metadata
        // C++ AGC defaults to no_raw_groups=0, meaning all groups use `get` (LZ path)
        // But we're writing raw segments, so disable LZ to use raw-only path
        let use_lz_encoding = false; // All groups are raw-only

        // Set reference if not set yet
        if metadata.reference.is_none() && !segments.is_empty() {
            metadata.reference = Some(segments[0].data.clone());
        }

        const PACK_CARDINALITY: usize = 50;

        // Pack and write segments
        for pack_segments in segments.chunks(PACK_CARDINALITY) {
            let mut packed_data = Vec::new();

            for (idx_in_pack, seg_info) in pack_segments.iter().enumerate() {
                let global_in_group_id = metadata.segments_written + idx_in_pack;

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} adding segment: contig={}, seg_part={}, data_len={}",
                        group_id, seg_info.contig_name, seg_info.seg_part_no, seg_info.data.len());
                }

                let contig_data = if !use_lz_encoding || global_in_group_id == 0 {
                    // Raw segment
                    seg_info.data.clone()
                } else {
                    // LZ-encoded segment
                    let mut lz_diff = LZDiff::new(self.config.min_match_len);
                    if let Some(ref reference) = metadata.reference {
                        lz_diff.prepare(reference);
                        lz_diff.encode(&seg_info.data)
                    } else {
                        seg_info.data.clone()
                    }
                };

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} contig_data len after encoding: {}", group_id, contig_data.len());
                }

                packed_data.extend_from_slice(&contig_data);
                packed_data.push(CONTIG_SEPARATOR);

                // Register in collection
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

            metadata.segments_written += pack_segments.len();
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
