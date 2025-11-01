// Streaming AGC Compressor
// Memory-efficient implementation that flushes groups incrementally

#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::needless_range_loop)]

/// Helper to estimate memory usage of a data structure
#[allow(dead_code)]
fn estimate_memory_mb<T>(items: &[T]) -> f64 {
    std::mem::size_of_val(items) as f64 / (1024.0 * 1024.0)
}

use crate::{
    contig_iterator::ContigIterator,
    genome_io::GenomeIO,
    kmer::{reverse_complement, Kmer, KmerMode},
    lz_diff::LZDiff,
    priority_queue::{BoundedPriorityQueue, PopResult},
    segment::{split_at_splitters_with_size, MISSING_KMER},
    segment_compression::{compress_reference_segment, compress_segment_configured},
    splitters::determine_splitters,
};
use anyhow::{Context, Result};
use crossbeam::channel::bounded;
use dashmap::DashMap;
use ragc_common::{
    stream_delta_name, stream_ref_name, Archive, CollectionV3, Contig, AGC_FILE_MAJOR,
    AGC_FILE_MINOR, CONTIG_SEPARATOR,
};
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicU32, AtomicUsize, Ordering},
    Arc, Barrier, Mutex, RwLock,
};
use std::thread;

/// Configuration for the streaming compressor
#[derive(Debug, Clone)]
pub struct StreamingCompressorConfig {
    pub kmer_length: u32,
    pub segment_size: u32,
    pub min_match_len: u32,
    pub compression_level: i32,
    pub verbosity: u32,
    /// Flush groups to disk when they reach this size (0 = only flush at finalize)
    pub group_flush_threshold: usize,
    /// Flush all groups after processing this many contigs (0 = never auto-flush)
    pub periodic_flush_interval: usize,
    /// Number of worker threads (matching C++ AGC: no_threads or no_threads-1 if >8)
    pub num_threads: usize,
    /// Adaptive mode: find new splitters for samples that can't be segmented well (matches C++ AGC -a flag)
    pub adaptive_mode: bool,
    /// Concatenated genomes mode: treat all contigs as one continuous sample (disables split detection)
    pub concatenated_genomes: bool,
}

impl Default for StreamingCompressorConfig {
    fn default() -> Self {
        // Matching C++ AGC line 2160: no_workers = (no_threads < 8) ? no_threads : no_threads - 1
        let num_cpus = num_cpus::get();
        let num_threads = if num_cpus < 8 { num_cpus } else { num_cpus - 1 };

        StreamingCompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            compression_level: 17,
            verbosity: 1,
            // MEMORY OPTIMIZATION: Flush segments to disk immediately
            // The middle splitter optimization only needs group_terminators (k-mer pairings),
            // not the actual segment data. C++ AGC writes to disk immediately via v_segments[group_id]->add()
            group_flush_threshold: 1, // Flush immediately (flush_group batches internally)
            periodic_flush_interval: 0, // Not needed with immediate flushing
            num_threads,
            adaptive_mode: false, // Default matches C++ AGC (adaptive mode off)
            concatenated_genomes: false, // Default: normal mode (split detection enabled)
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

/// A buffered segment awaiting two-phase processing (analogous to C++ seg_part_t)
#[derive(Debug, Clone)]
struct BufferedSegment {
    key: SegmentGroupKey,
    segment: SegmentInfo,
    is_new: bool, // true = creates new group (NEW), false = joins existing (KNOWN)
}

/// An uncompressed pack ready for compression (sent from workers to compression pool)
#[derive(Debug)]
struct UncompressedPack {
    group_id: u32,
    stream_id: usize,
    uncompressed_data: Vec<u8>, // LZ-encoded segments concatenated (NOT yet ZSTD compressed)
    segments: Vec<SegmentMetadata>, // Metadata for collection registration
}

/// A compressed pack ready to write to archive (sent from compression pool to writer)
#[derive(Debug)]
struct CompressedPack {
    group_id: u32,
    stream_id: usize,
    compressed_data: Vec<u8>, // Compressed pack (with marker byte if compressed)
    uncompressed_size: u64,   // Original packed size (0 if writing uncompressed)
    segments: Vec<SegmentMetadata>, // Metadata for collection registration
}

/// Either an uncompressed or compressed pack (returned from add_segment)
#[derive(Debug)]
enum PackToWrite {
    Compressed(CompressedPack),     // Reference packs (compressed immediately)
    Uncompressed(UncompressedPack), // Buffered packs (will be compressed by thread pool)
}

/// Segment metadata for collection registration
#[derive(Debug, Clone)]
struct SegmentMetadata {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    in_group_id: u32,
    is_rev_comp: bool,
    data_len: u32,
}

/// Prepared segment with normalized key, ready to add to group buffer
struct PreparedSegment {
    key: SegmentGroupKey,
    segment: SegmentInfo,
}

/// Contig task for C++ AGC-style contig-level parallelism
/// Workers pull contigs from queue, segment them, and add to shared groups
#[derive(Clone, PartialEq, Eq)]
enum ContigTask {
    /// Normal contig processing
    Normal {
        sample_name: String,
        contig_name: String,
        sequence: Contig,
        seq_num: u64,  // Sequence number to preserve file order
    },
    /// Synchronization point: workers wait, thread 0 registers groups
    Registration,
}

/// Pending segment awaiting registration
/// Matches C++ AGC's seg_part_t structure (agc_compressor.h:29-120)
/// Workers create these during segmentation, thread 0 registers them at barriers
#[derive(Clone)]
struct PendingSegment {
    key: SegmentGroupKey,
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    segment_data: Vec<u8>,
    kmer_front: u64,
    kmer_back: u64,
    seq_num: u64,  // Sequence number from parent ContigTask (preserves file order)
}

/// Segment group identified by flanking k-mers
/// IMPORTANT: Keys are normalized so kmer_front <= kmer_back (matching C++ minmax logic)
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct SegmentGroupKey {
    kmer_front: u64,
    kmer_back: u64,
}

impl SegmentGroupKey {
    /// Create a normalized key (matching C++ pk = minmax(kmer1, kmer2) logic)
    /// Returns (key, needs_rc) where needs_rc indicates if segment should be reverse complemented
    ///
    /// CRITICAL: C++ only normalizes when BOTH k-mers are present (lines 1313-1324)
    /// When one k-mer is MISSING, it uses a different code path (lines 1325-1372)
    fn new_normalized(kmer_front: u64, kmer_back: u64) -> (Self, bool) {
        // Only normalize if BOTH k-mers are valid (not MISSING_KMER)
        if kmer_front != MISSING_KMER && kmer_back != MISSING_KMER {
            // Both present: normalize (C++ lines 1317-1323)
            if kmer_front < kmer_back {
                (
                    SegmentGroupKey {
                        kmer_front,
                        kmer_back,
                    },
                    false,
                )
            } else {
                // Swap and set needs_rc=true (matching C++ lines 1319-1323)
                (
                    SegmentGroupKey {
                        kmer_front: kmer_back,
                        kmer_back: kmer_front,
                    },
                    true,
                )
            }
        } else {
            // One or both are MISSING: don't normalize (C++ lines 1297-1372)
            (
                SegmentGroupKey {
                    kmer_front,
                    kmer_back,
                },
                false,
            )
        }
    }
}

/// Group metadata for tracking what's been written
struct GroupMetadata {
    group_id: u32,
    stream_id: usize,                   // Delta stream ID
    ref_stream_id: Option<usize>,       // Reference stream ID (for LZ groups)
    reference: Option<Contig>,          // First segment (for LZ encoding)
    ref_written: bool,                  // Whether reference has been written to archive
    segments_written: usize, // Number of segments written to archive (not including buffered)
    pending_segments: Vec<SegmentInfo>, // Buffered segments waiting for complete pack
    is_flushed: bool,
}

/// Per-group writer state for parallel segment writing
/// Each group has its own buffer and writes packs when buffer reaches PACK_CARDINALITY
struct GroupWriter {
    group_id: u32,
    stream_id: usize,
    ref_stream_id: Option<usize>,
    reference: Option<Contig>,
    ref_written: bool,
    segments_written: usize,
    pending_segments: Vec<SegmentInfo>,
}

impl GroupWriter {
    fn new(group_id: u32, stream_id: usize, ref_stream_id: Option<usize>) -> Self {
        GroupWriter {
            group_id,
            stream_id,
            ref_stream_id,
            reference: None,
            ref_written: false,
            segments_written: 0,
            pending_segments: Vec::new(),
        }
    }

    /// Get reference segment data if available (for split cost calculation)
    fn get_reference(&self) -> Option<&Contig> {
        self.reference.as_ref()
    }

    /// Add a segment to this group's buffer
    /// Returns Some(CompressedPack) if buffer is full and pack should be written
    /// Returns Some(pack) if this is the first segment of an LZ group or if buffer is full
    fn add_segment(
        &mut self,
        segment: SegmentInfo,
        config: &StreamingCompressorConfig,
    ) -> Result<Option<PackToWrite>> {
        const PACK_CARDINALITY: usize = 50;
        const NO_RAW_GROUPS: u32 = 16;

        let use_lz_encoding = self.group_id >= NO_RAW_GROUPS;

        // For LZ groups: prepare reference pack on first segment
        if use_lz_encoding && !self.ref_written {
            self.reference = Some(segment.data.clone());
            self.ref_written = true;

            if let Some(ref_stream_id) = self.ref_stream_id {
                // Match C++ AGC's compression strategy (segment.h:218-256)
                // Check repetitiveness and use tuple packing for non-repetitive data
                let (compressed, marker_byte) = compress_reference_segment(&segment.data)?;

                let (compressed_data, uncompressed_size) =
                    if compressed.len() + 1 < segment.data.len() {
                        let mut compressed_with_marker = compressed;
                        compressed_with_marker.push(marker_byte);
                        (compressed_with_marker, segment.data.len() as u64)
                    } else {
                        (segment.data.clone(), 0)
                    };

                // Create metadata for reference segment
                let segments = vec![SegmentMetadata {
                    sample_name: segment.sample_name.clone(),
                    contig_name: segment.contig_name.clone(),
                    seg_part_no: segment.seg_part_no,
                    in_group_id: 0,
                    is_rev_comp: segment.is_rev_comp,
                    data_len: segment.data.len() as u32,
                }];

                return Ok(Some(PackToWrite::Compressed(CompressedPack {
                    group_id: self.group_id,
                    stream_id: ref_stream_id,
                    compressed_data,
                    uncompressed_size,
                    segments,
                })));
            }
        }

        // Buffer the segment
        self.pending_segments.push(segment);

        // Prepare pack if buffer is full
        if self.pending_segments.len() >= PACK_CARDINALITY {
            return Ok(self.prepare_pack(config)?.map(PackToWrite::Uncompressed));
        }

        Ok(None)
    }

    /// Prepare an UNCOMPRESSED pack from buffered segments (LZ encoding only, NO ZSTD compression)
    /// Compression will happen in separate compression thread pool
    fn prepare_pack(
        &mut self,
        config: &StreamingCompressorConfig,
    ) -> Result<Option<UncompressedPack>> {
        if self.pending_segments.is_empty() {
            return Ok(None);
        }

        const NO_RAW_GROUPS: u32 = 16;
        let use_lz_encoding = self.group_id >= NO_RAW_GROUPS;

        // Build pack (LZ encoding + concatenation)
        let mut packed_data = Vec::new();

        if config.verbosity > 2 {
            eprintln!(
                "[PACK] Creating pack for group {} with {} segments",
                self.group_id,
                self.pending_segments.len()
            );
        }

        // CRITICAL FIX: Reuse ONE LZDiff instance for all segments in this pack!
        // Creating a new LZDiff per segment causes massive allocations (116 GB for yeast235)
        // because each prepare() clones reference and rebuilds the entire HashMap
        let mut lz_diff = if use_lz_encoding {
            let mut lzd = LZDiff::new(config.min_match_len);
            if let Some(ref reference) = self.reference {
                lzd.prepare(reference);
            }
            Some(lzd)
        } else {
            None
        };

        for (idx, seg_info) in self.pending_segments.iter().enumerate() {
            let contig_data = if let Some(ref mut lzd) = lz_diff {
                lzd.encode(&seg_info.data)
            } else {
                seg_info.data.clone()
            };

            if config.verbosity > 2 {
                eprintln!(
                    "[PACK]   Segment {}: {} bytes ({}:{}), adding separator",
                    idx,
                    contig_data.len(),
                    seg_info.sample_name,
                    seg_info.contig_name
                );
            }

            packed_data.extend_from_slice(&contig_data);
            packed_data.push(CONTIG_SEPARATOR);
        }

        if config.verbosity > 2 {
            eprintln!("[PACK] Total pack size: {} bytes", packed_data.len());
        }

        // Build segment metadata (NO COMPRESSION HERE!)
        let mut segments = Vec::new();
        for (idx_in_pack, seg_info) in self.pending_segments.iter().enumerate() {
            let global_in_group_id = if use_lz_encoding {
                self.segments_written + idx_in_pack + 1 // +1 because reference is at 0
            } else {
                self.segments_written + idx_in_pack
            };

            segments.push(SegmentMetadata {
                sample_name: seg_info.sample_name.clone(),
                contig_name: seg_info.contig_name.clone(),
                seg_part_no: seg_info.seg_part_no,
                in_group_id: global_in_group_id as u32,
                is_rev_comp: seg_info.is_rev_comp,
                data_len: seg_info.data.len() as u32,
            });
        }

        self.segments_written += self.pending_segments.len();
        self.pending_segments.clear();

        // Return UNCOMPRESSED pack - compression will happen in compression thread pool!
        Ok(Some(UncompressedPack {
            group_id: self.group_id,
            stream_id: self.stream_id,
            uncompressed_data: packed_data,
            segments,
        }))
    }

    /// Flush any remaining buffered segments as a final pack
    /// Called at the end of compression to ensure all segments are written
    fn flush_final(
        &mut self,
        config: &StreamingCompressorConfig,
    ) -> Result<Option<PackToWrite>> {
        // Prepare pack from any remaining segments (even if less than PACK_CARDINALITY)
        Ok(self.prepare_pack(config)?.map(PackToWrite::Uncompressed))
    }
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

    // K-mer terminator tracking (for missing middle splitter optimization)
    // Maps each k-mer to the sorted vector of k-mers it pairs with in existing groups
    // This enables splitting segments to reuse existing groups
    // NOTE: Matches C++ AGC approach - single RwLock protecting entire map (seg_map_mtx)
    group_terminators: HashMap<u64, Vec<u64>>,

    // Adaptive mode: reference k-mers for exclusion (matching C++ AGC v_candidate_kmers and v_duplicated_kmers)
    reference_kmers_singletons: HashSet<u64>, // Singleton k-mers from reference (candidates)
    reference_kmers_duplicates: HashSet<u64>, // Duplicate k-mers from reference (non-singletons)
    current_splitters: HashSet<u64>,          // Mutable splitter set (grows in adaptive mode)

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
            next_group_id: 16, // Reserve 0-15 for raw groups (matching C++ AGC)
            group_terminators: HashMap::new(),
            reference_kmers_singletons: HashSet::new(),
            reference_kmers_duplicates: HashSet::new(),
            current_splitters: HashSet::new(),
            total_bases_processed: 0,
            total_segments: 0,
            total_groups_flushed: 0,
            contigs_since_flush: 0,
        })
    }

    /// Find new splitters for contigs that can't be segmented well (adaptive mode)
    ///
    /// This implements C++ AGC's find_new_splitters function (lines 2046-2077):
    /// 1. Extract all k-mers from contig
    /// 2. Find singletons
    /// 3. Exclude k-mers from reference genome (both singletons and duplicates)
    /// 4. Return remaining k-mers as new splitters
    ///
    /// # Arguments
    /// * `contig` - The contig sequence to find new splitters from
    ///
    /// # Returns
    /// HashSet of new splitter k-mers (not in reference)
    fn find_new_splitters(&self, contig: &Contig) -> HashSet<u64> {
        use crate::kmer_extract::enumerate_kmers;
        use rdst::RadixSort;

        // 1. Extract all k-mers from this contig
        let mut contig_kmers = enumerate_kmers(contig, self.config.kmer_length as usize);

        // 2. Sort and find singletons
        contig_kmers.radix_sort_unstable();
        crate::kmer_extract::remove_non_singletons(&mut contig_kmers, 0);

        // 3. Exclude k-mers from reference genome
        let mut new_splitters = HashSet::new();

        for kmer in contig_kmers {
            // Exclude if in reference singletons OR reference duplicates
            if !self.reference_kmers_singletons.contains(&kmer)
                && !self.reference_kmers_duplicates.contains(&kmer)
            {
                new_splitters.insert(kmer);
            }
        }

        new_splitters
    }

    /// Add a FASTA file to the archive (streaming mode)

    /// Detect if a FASTA file contains multiple samples in headers
    /// (format: >sample#haplotype#chromosome)
    pub fn detect_multi_sample_fasta(fasta_path: &Path) -> Result<bool> {
        let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
            .context("Failed to open FASTA for detection")?;

        // Read first few headers and collect unique sample names
        let mut sample_names = std::collections::HashSet::new();
        let mut has_multi_part_headers = false;

        for _ in 0..10 {
            // Check more headers to be sure
            if let Some((header, sample_name, _, _)) = reader.read_contig_with_sample()? {
                let parts: Vec<&str> = header.split('#').collect();
                if parts.len() >= 3 {
                    has_multi_part_headers = true;
                    sample_names.insert(sample_name);

                    // If we found 2+ unique samples, it's definitely multi-sample
                    if sample_names.len() >= 2 {
                        return Ok(true);
                    }
                }
            } else {
                break;
            }
        }

        // Only multi-sample if we found the right header format AND multiple samples
        Ok(has_multi_part_headers && sample_names.len() >= 2)
    }

    /// Add a multi-sample FASTA file with splitter-based segmentation
    /// This handles files where samples are encoded in headers (>sample#hap#chr format)
    pub fn add_multi_sample_fasta_with_splitters(&mut self, fasta_path: &Path) -> Result<()> {
        eprintln!("DEBUG: Entered add_multi_sample_fasta_with_splitters");
        if self.config.verbosity > 0 {
            println!("=== Processing multi-sample FASTA (grouping by sample names in headers) ===");
            println!("Input: {fasta_path:?}");
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
            while let Some((_full_header, sample_name, _contig_name, sequence)) =
                reader.read_contig_with_sample()?
            {
                if !sequence.is_empty() {
                    if reference_sample.is_empty() {
                        reference_sample = sample_name.clone();
                        if self.config.verbosity > 0 {
                            println!("Using first sample ({reference_sample}) as reference");
                        }
                    }

                    // CRITICAL: Only collect contigs from FIRST sample (matching C++ AGC behavior)
                    // C++ AGC's determine_splitters() only processes the reference file/sample
                    if sample_name == reference_sample {
                        reference_contigs.push(sequence);
                        contig_count += 1;
                    }
                }
            }

            if self.config.verbosity > 0 {
                println!("Collected {contig_count} reference contigs from first sample ({reference_sample})");
            }
        }

        // Find splitters (also get singleton and duplicate k-mers for adaptive mode)
        let (splitters, singletons, duplicates) = determine_splitters(
            &reference_contigs,
            self.config.kmer_length as usize,
            self.config.segment_size as usize,
        );

        // Store reference k-mers for adaptive mode
        if self.config.adaptive_mode {
            self.reference_kmers_singletons = singletons;
            self.reference_kmers_duplicates = duplicates;
            self.current_splitters = splitters.clone();

            if self.config.verbosity > 0 {
                println!("Adaptive mode: Stored {} singleton k-mers and {} duplicate k-mers from reference",
                    self.reference_kmers_singletons.len(), self.reference_kmers_duplicates.len());
            }
        }

        if self.config.verbosity > 0 {
            println!("Found {} actually-used splitter k-mers", splitters.len());

            // DEBUG: Output first 20 splitters sorted
            let mut splitter_vec: Vec<u64> = splitters.iter().copied().collect();
            splitter_vec.sort_unstable();
            eprintln!("RUST_SPLITTERS: First 20 splitters:");
            for (i, &spl) in splitter_vec.iter().take(20).enumerate() {
                eprintln!("  {i}: {spl}");
            }

            println!();
            println!("=== Pass 2: Segmenting and compressing all samples ===");
        }

        drop(reference_contigs); // Free memory

        // Pass 2: THREE-PHASE parallelization (matching C++ AGC architecture!)
        // Phase 1: Single-threaded segmentation (collect all segments)
        // Phase 2: Single-threaded grouping (group segments by key)
        // Phase 3: Parallel group processing (each worker gets exclusive groups - NO MUTEX!)
        //
        // This eliminates the mutex contention from the producer-consumer pattern where
        // workers competed for the same ~770 groups.

        if self.config.verbosity > 0 {
            println!(
                "=== Pass 2: Group-level parallelization with {} worker threads ===",
                self.config.num_threads
            );
        }

        // ==================================================================
        // PHASE 1: Single-threaded segmentation (collect all segments)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 1: Segmenting all contigs (single-threaded)...");
        }

        let mut all_segments: Vec<PreparedSegment> = Vec::new();
        let mut total_segments_count = 0;
        let mut total_bases_count = 0;
        let mut contig_count = 0;

        {
            let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                .context("Failed to open multi-sample FASTA for segmentation")?;

            while let Some((full_header, sample_name, _contig_name, sequence)) =
                reader.read_contig_with_sample()?
            {
                if !sequence.is_empty() {
                    contig_count += 1;

                    // PRE-REGISTER CONTIG: Ensures contigs are registered in file order
                    // before parallel Phase 3 processing (fixes data corruption!)
                    self.collection
                        .register_sample_contig(&sample_name, &full_header)?;

                    // Segment at splitters
                    let segments = split_at_splitters_with_size(
                        &sequence,
                        &splitters,
                        self.config.kmer_length as usize,
                        self.config.segment_size as usize,
                    );

                    // Process each segment: normalize and prepare for grouping
                    for (seg_idx, segment) in segments.into_iter().enumerate() {
                        total_bases_count += segment.data.len();
                        total_segments_count += 1;

                        // Prepare segment info (k-mer normalization, RC if needed)
                        // IMPORTANT: Use full_header as contig name to match C++ AGC

                        // DEBUG: Print k-mers from first few segments
                        if total_segments_count < 10 {
                            eprintln!(
                                "DEBUG segment {}: front_kmer={}, back_kmer={}, len={}",
                                total_segments_count,
                                segment.front_kmer,
                                segment.back_kmer,
                                segment.data.len()
                            );
                        }

                        let seg_info = Self::prepare_segment_info(
                            &self.config,
                            &sample_name,
                            &full_header, // Use full header instead of parsed contig_name
                            seg_idx,
                            segment.data,
                            Some(segment.front_kmer),
                            Some(segment.back_kmer),
                        )?;

                        // DEBUG: Print key from first few prepared segments
                        if total_segments_count < 10 {
                            eprintln!(
                                "DEBUG prepared {}: key.front={}, key.back={}",
                                total_segments_count,
                                seg_info.key.kmer_front,
                                seg_info.key.kmer_back
                            );
                        }

                        all_segments.push(seg_info);
                    }

                    if self.config.verbosity > 1 && contig_count % 1000 == 0 {
                        eprintln!("Phase 1: Segmented {contig_count} contigs ({total_segments_count} segments)");
                    }
                }
            }
        }

        if self.config.verbosity > 0 {
            println!(
                "Phase 1 complete: {contig_count} contigs → {total_segments_count} segments ({total_bases_count} bases)"
            );
        }

        // ==================================================================
        // PHASE 1.5: Fix one-kmer segment keys using candidate search
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 1.5: Finding candidate groups for one-kmer segments...");
        }

        // Build terminators map from segments with both k-mers
        let mut group_terminators: HashMap<u64, Vec<u64>> = HashMap::new();
        for segment in &all_segments {
            if segment.key.kmer_front != MISSING_KMER && segment.key.kmer_back != MISSING_KMER {
                // Track both directions
                group_terminators
                    .entry(segment.key.kmer_front)
                    .or_default()
                    .push(segment.key.kmer_back);
                if segment.key.kmer_front != segment.key.kmer_back {
                    group_terminators
                        .entry(segment.key.kmer_back)
                        .or_default()
                        .push(segment.key.kmer_front);
                }
            }
        }

        // Fix one-kmer segments to use first candidate
        let mut fixed_count = 0;
        for segment in &mut all_segments {
            let has_front = segment.key.kmer_front != MISSING_KMER;
            let has_back = segment.key.kmer_back != MISSING_KMER;

            if (has_front && !has_back) || (!has_front && has_back) {
                // One k-mer present - search for candidates
                let present_kmer = if has_front {
                    segment.key.kmer_front
                } else {
                    segment.key.kmer_back
                };

                if let Some(candidates) = group_terminators.get(&present_kmer) {
                    if !candidates.is_empty() {
                        // Use FIRST candidate to avoid creating new MISSING_KMER groups
                        let first_cand = candidates[0];
                        let old_key = segment.key.clone();
                        let new_key = if first_cand < present_kmer {
                            SegmentGroupKey {
                                kmer_front: first_cand,
                                kmer_back: present_kmer,
                            }
                        } else {
                            SegmentGroupKey {
                                kmer_front: present_kmer,
                                kmer_back: first_cand,
                            }
                        };

                        segment.key = new_key;
                        fixed_count += 1;
                    }
                }
            }
        }

        if self.config.verbosity > 0 {
            println!("Phase 1.5 complete: Fixed {fixed_count} one-kmer segments");
        }

        // ==================================================================
        // PHASE 2: Single-threaded grouping (group segments by key)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 2: Grouping segments by k-mer pairs (single-threaded)...");
        }

        let mut groups: HashMap<SegmentGroupKey, Vec<PreparedSegment>> = HashMap::new();
        for segment in all_segments {
            groups.entry(segment.key.clone()).or_default().push(segment);
        }

        let num_groups = groups.len();
        if self.config.verbosity > 0 {
            println!("Phase 2 complete: {num_groups} unique groups");
        }

        // Convert to Vec for partitioning across workers
        let groups_vec: Vec<(SegmentGroupKey, Vec<PreparedSegment>)> = groups.into_iter().collect();

        // ==================================================================
        // PHASE 3: Parallel group processing (matching C++ AGC architecture!)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!(
                "Phase 3: Processing groups in parallel ({} threads)...",
                self.config.num_threads
            );
        }

        // Group ID counter for sequential assignment
        let next_group_id = Arc::new(AtomicUsize::new(self.next_group_id as usize));

        // Clone config so we don't borrow self in the parallel closure!
        let config = self.config.clone();

        // Configure Rayon thread pool to respect num_threads setting
        // CRITICAL: Default par_iter() uses ALL cores, ignoring our config!
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.num_threads)
            .build()
            .expect("Failed to build thread pool");

        // Process groups in parallel using Rayon (like C++ AGC's parallel loop)
        // Each worker gets exclusive groups - NO CHANNELS, NO MUTEX!
        // MEMORY FIX: Use into_par_iter() to consume and move segments instead of cloning!
        let all_packs: Vec<CompressedPack> = pool.install(|| {
            groups_vec
                .into_par_iter()
                .flat_map(|(_key, segments)| {
                    // DEBUG: Check if we're actually running on multiple threads
                    eprintln!(
                        "Thread {:?} processing a group",
                        std::thread::current().id()
                    );

                    // Assign group ID sequentially
                    let gid = next_group_id.fetch_add(1, Ordering::SeqCst) as u32;
                    let stream_id = gid as usize;

                    const NO_RAW_GROUPS: u32 = 16;
                    let ref_stream_id = if gid >= NO_RAW_GROUPS {
                        Some(100000 + gid as usize) // Ref streams at offset 100000
                    } else {
                        None
                    };

                    // Create group writer for this group (thread-local, no contention!)
                    let mut group_writer = GroupWriter::new(gid, stream_id, ref_stream_id);
                    let mut packs = Vec::new();

                    // Process all segments in this group
                    for prepared_seg in segments {
                        // MEMORY FIX: Move segment instead of cloning (saves 220 MB!)
                        if let Ok(Some(pack)) =
                            group_writer.add_segment(prepared_seg.segment, &config)
                        {
                            match pack {
                                PackToWrite::Compressed(compressed_pack) => {
                                    packs.push(compressed_pack);
                                }
                                PackToWrite::Uncompressed(uncompressed_pack) => {
                                    // Compress immediately (parallel ZSTD!)
                                    if let Ok(compressed_data) = compress_segment_configured(
                                        &uncompressed_pack.uncompressed_data,
                                        config.compression_level,
                                    ) {
                                        let (final_data, uncompressed_size) = if compressed_data
                                            .len()
                                            + 1
                                            < uncompressed_pack.uncompressed_data.len()
                                        {
                                            let mut data_with_marker = compressed_data;
                                            data_with_marker.push(0);
                                            (
                                                data_with_marker,
                                                uncompressed_pack.uncompressed_data.len() as u64,
                                            )
                                        } else {
                                            (uncompressed_pack.uncompressed_data, 0)
                                        };

                                        packs.push(CompressedPack {
                                            group_id: uncompressed_pack.group_id,
                                            stream_id: uncompressed_pack.stream_id,
                                            compressed_data: final_data,
                                            uncompressed_size,
                                            segments: uncompressed_pack.segments,
                                        });
                                    }
                                }
                            }
                        }
                    }

                    // Flush remaining segments for this group
                    if let Ok(Some(uncompressed_pack)) = group_writer.prepare_pack(&config) {
                        if let Ok(compressed_data) = compress_segment_configured(
                            &uncompressed_pack.uncompressed_data,
                            config.compression_level,
                        ) {
                            let (final_data, uncompressed_size) = if compressed_data.len() + 1
                                < uncompressed_pack.uncompressed_data.len()
                            {
                                let mut data_with_marker = compressed_data;
                                data_with_marker.push(0);
                                (
                                    data_with_marker,
                                    uncompressed_pack.uncompressed_data.len() as u64,
                                )
                            } else {
                                (uncompressed_pack.uncompressed_data, 0)
                            };

                            packs.push(CompressedPack {
                                group_id: uncompressed_pack.group_id,
                                stream_id: uncompressed_pack.stream_id,
                                compressed_data: final_data,
                                uncompressed_size,
                                segments: uncompressed_pack.segments,
                            });
                        }
                    }

                    packs
                })
                .collect()
        });

        if self.config.verbosity > 0 {
            println!(
                "Phase 3 complete: processed {} groups → {} packs",
                num_groups,
                all_packs.len()
            );
        }

        // ==================================================================
        // PHASE 4: Sequential write (like C++ AGC's single-threaded registration)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 4: Writing to archive (sequential)...");
        }

        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
        let mut stream_id_map = HashMap::new();
        let mut packs_written = 0;
        let mut segments_written = 0;

        for pack in all_packs {
            // Register stream if first time seeing it and get actual stream_id
            let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                id
            } else {
                let stream_name = if pack.stream_id >= 100000 {
                    // Reference stream
                    let group_id = (pack.stream_id - 100000) as u32;
                    stream_ref_name(archive_version, group_id)
                } else {
                    // Delta stream
                    stream_delta_name(archive_version, pack.group_id)
                };
                let actual_id = self.archive.register_stream(&stream_name);
                stream_id_map.insert(pack.stream_id, actual_id);
                actual_id
            };

            // Register segments in collection
            for seg_meta in &pack.segments {
                // Register contig if first time
                self.collection
                    .register_sample_contig(&seg_meta.sample_name, &seg_meta.contig_name)?;

                // Register segment placement
                self.collection.add_segment_placed(
                    &seg_meta.sample_name,
                    &seg_meta.contig_name,
                    seg_meta.seg_part_no,
                    pack.group_id,
                    seg_meta.in_group_id,
                    seg_meta.is_rev_comp,
                    seg_meta.data_len,
                )?;
            }

            // Write pack to archive using actual stream_id
            self.archive.add_part(
                actual_stream_id,
                &pack.compressed_data,
                pack.uncompressed_size,
            )?;

            packs_written += 1;
            segments_written += pack.segments.len();

            if self.config.verbosity > 1 && packs_written % 100 == 0 {
                eprintln!("Wrote {packs_written} packs ({segments_written} segments)");
            }
        }

        if self.config.verbosity > 0 {
            println!(
                "Processing complete: wrote {packs_written} packs ({segments_written} segments)"
            );
        }

        // Update state with counts from Phase 1
        self.total_segments = total_segments_count;
        self.total_bases_processed = total_bases_count;

        // Update next_group_id for subsequent operations
        self.next_group_id = Arc::try_unwrap(next_group_id)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap next_group_id Arc"))?
            .into_inner() as u32;

        if self.config.verbosity > 0 {
            println!("Processed {contig_count} contigs total");
        }

        Ok(())
    }

    /// Unified method: Add contigs using an iterator (single file or multiple files)
    ///
    /// This is the unified input handling approach that works identically for both:
    /// - Single pansn-format FASTA file (via PansnFileIterator)
    /// - Multiple per-sample FASTA files (via MultiFileIterator)
    ///
    /// Uses the two-pass algorithm with parallel streaming:
    /// Pass 1: Collect reference contigs → find splitter k-mers
    /// Pass 2: Stream through all contigs → segment → parallel compress
    /// Process contigs using a ContigIterator (for multi-sample FASTAs)
    ///
    /// This is the entry point for single-file multi-sample FASTAs (e.g., yeast235.fa.bgz)
    /// It finds splitters from the first sample, then processes all samples using
    /// the unified inline split logic.
    pub fn add_contigs_with_splitters(
        &mut self,
        mut iterator: Box<dyn ContigIterator>,
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("=== Pass 1: Finding splitter k-mers from first sample ===");
        }

        // Pass 1: Collect reference contigs for splitter finding
        let mut reference_contigs = Vec::new();
        let mut reference_sample = String::new();
        let mut contig_count = 0;

        while let Some((sample_name, _contig_name, sequence)) = iterator.next_contig()? {
            if !sequence.is_empty() {
                if reference_sample.is_empty() {
                    reference_sample = sample_name.clone();
                    if self.config.verbosity > 0 {
                        println!("Using first sample ({reference_sample}) as reference");
                    }
                }

                // Only collect contigs from FIRST sample (matching C++ AGC)
                if sample_name == reference_sample {
                    reference_contigs.push(sequence);
                    contig_count += 1;
                }
            }
        }

        if self.config.verbosity > 0 {
            println!("Collected {contig_count} reference contigs");
        }

        // Find splitters (also get singleton and duplicate k-mers for adaptive mode)
        let (splitters, singletons, duplicates) = determine_splitters(
            &reference_contigs,
            self.config.kmer_length as usize,
            self.config.segment_size as usize,
        );

        // Store reference k-mers for adaptive mode
        if self.config.adaptive_mode {
            self.reference_kmers_singletons = singletons;
            self.reference_kmers_duplicates = duplicates;
            self.current_splitters = splitters.clone();
        }

        if self.config.verbosity > 0 {
            println!("Found {} splitter k-mers", splitters.len());
            println!();
            println!("=== Pass 2: Processing all samples with unified logic ===");
        }

        drop(reference_contigs); // Free memory

        // Reset iterator for second pass
        iterator.reset()?;

        // Use the unified implementation
        self.add_segments_with_inline_splits(iterator, &splitters)?;

        Ok(())
    }

    /// Add multiple FASTA files with splitter-based segmentation (two-pass)
    ///
    /// This is the real AGC algorithm:
    /// Pass 1: Read all genomes → find singleton k-mers (splitters)
    /// Pass 2: Re-read genomes → segment at splitters → compress
    ///
    /// This achieves much better compression than treating each contig as one segment!
    pub fn add_fasta_files_with_splitters(
        &mut self,
        fasta_paths: &[(String, &Path)],
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("=== Pass 1: Finding splitter k-mers across all genomes ===");
        }

        // Pass 1: Stream through FIRST FILE ONLY to find splitters
        // (matching C++ AGC: only reference file used for splitter determination)
        // IMPORTANT: We read ALL contigs from the first file, finding singletons
        // within that file's contigs (not across all genomes)
        //
        // MEMORY-EFFICIENT: Streams through file twice instead of loading all contigs
        // For yeast (12MB genome): ~100MB Vec vs ~2.8GB loading all contigs!
        let (splitters, singletons, duplicates) = if let Some((ref_sample_name, ref_fasta_path)) =
            fasta_paths.first()
        {
            if self.config.verbosity > 0 {
                println!("Using first file ({ref_sample_name}) as reference to find splitters (streaming)...");
            }

            // Read ALL contigs from FIRST FILE ONLY (matching C++ AGC behavior)
            crate::splitters::determine_splitters_streaming(
                ref_fasta_path,
                self.config.kmer_length as usize,
                self.config.segment_size as usize,
            )?
        } else {
            anyhow::bail!("No input files provided");
        };

        // Store reference k-mers for adaptive mode
        if self.config.adaptive_mode {
            self.reference_kmers_singletons = singletons;
            self.reference_kmers_duplicates = duplicates;
            self.current_splitters = splitters.clone();

            if self.config.verbosity > 0 {
                println!("Adaptive mode: Stored {} singleton k-mers and {} duplicate k-mers from reference",
                    self.reference_kmers_singletons.len(), self.reference_kmers_duplicates.len());
            }
        }

        if self.config.verbosity > 0 {
            println!("Found {} actually-used splitter k-mers", splitters.len());
            println!();
            println!("=== Pass 2: Segmenting and compressing genomes ===");
        }

        // Pass 2: Use C++ AGC worker-based architecture (parallel contig processing)
        // Create MultiFileIterator from file paths
        let file_path_bufs: Vec<PathBuf> = fasta_paths
            .iter()
            .map(|(_, path)| path.to_path_buf())
            .collect();
        let iterator = Box::new(crate::contig_iterator::MultiFileIterator::new(
            file_path_bufs,
        )?);

        // Use worker-based processing (parallel segmentation!)
        self.add_segments_with_workers(iterator, &splitters)?;

        Ok(())
    }

    /// Find shared k-mer between two k-mers (C++ AGC find_cand_segment_with_missing_middle_splitter)
    /// Returns the first k-mer that appears in both connection lists
    fn find_shared_kmer(
        &self,
        kmer_front: u64,
        kmer_back: u64,
        kmer_connections: &HashMap<u64, Vec<u64>>,
    ) -> Option<u64> {
        let front_connections = kmer_connections.get(&kmer_front)?;
        let back_connections = kmer_connections.get(&kmer_back)?;

        // Find intersection (C++ AGC lines 1554-1569)
        // Both lists should be sorted from how we build them
        for &kmer in front_connections {
            if kmer != u64::MAX && back_connections.contains(&kmer) {
                return Some(kmer);
            }
        }

        None
    }

    /// C++ AGC exact architecture: SEQUENTIAL SAMPLE PROCESSING with inline split checking
    /// Process samples ONE AT A TIME (like C++ AGC) so kmer_connections builds naturally
    /// For each segment: check if we can split into existing groups BEFORE creating new group
    ///
    /// UNIFIED IMPLEMENTATION: Works with any ContigIterator (MultiFileIterator, BufferedPansnFileIterator, etc.)
    fn add_segments_with_inline_splits(
        &mut self,
        mut iterator: Box<dyn ContigIterator>,
        splitters: &HashSet<u64>,
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Starting C++ AGC-style compression (sequential sample processing)...");
        }

        // Track which segment keys we've seen (persistent across samples)
        let mut known_groups = HashMap::<SegmentGroupKey, u32>::new();

        // Track contigs and bases
        let mut total_contigs = 0;
        let mut total_bases = 0;
        let mut total_segments = 0;
        let mut current_sample = String::new();

        // Group storage and management
        let mut groups = HashMap::<SegmentGroupKey, (u32, GroupWriter)>::new();
        let mut group_terminators = HashMap::<u64, Vec<u64>>::new();
        let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));

        if self.config.verbosity > 0 {
            println!("=== Sequential Processing: Process segments one-at-a-time ===");
        }

        // ================================================================
        // SET UP WRITER THREAD (for Phases 2/3/4)
        // ================================================================
        let (pack_tx, pack_rx) = bounded::<CompressedPack>(10);
        let archive = Arc::new(Mutex::new(std::mem::replace(
            &mut self.archive,
            Archive::new_writer(),
        )));
        let collection = Arc::new(Mutex::new(std::mem::replace(
            &mut self.collection,
            CollectionV3::new(),
        )));

        let writer_handle = {
            let archive = archive.clone();
            let collection = collection.clone();
            let config = self.config.clone();

            thread::spawn(move || -> Result<usize> {
                let mut stream_id_map: HashMap<usize, usize> = HashMap::new();
                let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
                let mut packs_written = 0;
                let mut segments_registered = 0;
                let mut seen_segments: HashMap<(String, String, usize), (u32, u32)> =
                    HashMap::new();

                while let Ok(pack) = pack_rx.recv() {
                    let mut archive_lock = archive.lock().unwrap();
                    let mut collection_lock = collection.lock().unwrap();

                    if config.verbosity > 2 && (pack.group_id == 11813 || pack.group_id < 100) {
                        println!(
                            "  Writer received pack: group_id={}, stream_id={}, {} segments",
                            pack.group_id,
                            pack.stream_id,
                            pack.segments.len()
                        );
                    }

                    // Register stream on-demand
                    let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                        if config.verbosity > 2 && (pack.group_id == 11813 || pack.group_id < 100) {
                            println!(
                                "    Stream already registered: stream_id={} -> archive_id={}",
                                pack.stream_id, id
                            );
                        }
                        id
                    } else {
                        let stream_name = if pack.stream_id >= 100000 {
                            stream_ref_name(archive_version, (pack.stream_id - 100000) as u32)
                        } else {
                            stream_delta_name(archive_version, pack.group_id)
                        };
                        if config.verbosity > 1 {
                            println!(
                                "  Registering stream: {} for group_id={}",
                                stream_name, pack.group_id
                            );
                        }
                        let id = archive_lock.register_stream(&stream_name);

                        // Check for collision before inserting
                        if let Some(&existing_id) = stream_id_map.get(&pack.stream_id) {
                            eprintln!("ERROR: stream_id collision detected!");
                            eprintln!(
                                "  stream_id={} already maps to archive_id={}",
                                pack.stream_id, existing_id
                            );
                            eprintln!(
                                "  Trying to register stream '{}' for group_id={} as archive_id={}",
                                stream_name, pack.group_id, id
                            );
                            eprintln!(
                                "  This means two different groups are using the same stream_id!"
                            );
                        }

                        stream_id_map.insert(pack.stream_id, id);

                        if config.verbosity > 2 {
                            println!(
                                "    Mapped: stream_id={} -> archive_id={} ({})",
                                pack.stream_id, id, stream_name
                            );
                        }
                        id
                    };

                    // Register segments in collection (seg_part_no already correct from sequential processing!)
                    for seg_meta in &pack.segments {
                        // Check for duplicate seg_part_no
                        let key = (
                            seg_meta.sample_name.clone(),
                            seg_meta.contig_name.clone(),
                            seg_meta.seg_part_no,
                        );
                        if let Some((prev_group_id, prev_in_group_id)) = seen_segments.get(&key) {
                            println!("WARNING: Duplicate segment detected!");
                            println!(
                                "  Sample: {}, Contig: {}, seg_part_no: {}",
                                seg_meta.sample_name, seg_meta.contig_name, seg_meta.seg_part_no
                            );
                            println!(
                                "  Previous: group_id={}, in_group_id={}",
                                prev_group_id, prev_in_group_id
                            );
                            println!(
                                "  Current:  group_id={}, in_group_id={}",
                                pack.group_id, seg_meta.in_group_id
                            );
                        }

                        if pack.group_id > 50000 {
                            println!("WARNING: Registering segment with huge group_id={}! Sample={}, Contig={}, seg_part_no={}",
                                pack.group_id, seg_meta.sample_name, seg_meta.contig_name, seg_meta.seg_part_no);
                        }

                        seen_segments.insert(key, (pack.group_id, seg_meta.in_group_id));

                        collection_lock
                            .register_sample_contig(&seg_meta.sample_name, &seg_meta.contig_name)?;
                        collection_lock.add_segment_placed(
                            &seg_meta.sample_name,
                            &seg_meta.contig_name,
                            seg_meta.seg_part_no,
                            pack.group_id,
                            seg_meta.in_group_id,
                            seg_meta.is_rev_comp,
                            seg_meta.data_len,
                        )?;
                        segments_registered += 1;
                    }

                    // Write pack
                    archive_lock.add_part(
                        actual_stream_id,
                        &pack.compressed_data,
                        pack.uncompressed_size,
                    )?;
                    packs_written += 1;
                }

                if config.verbosity > 0 {
                    println!(
                        "Writer thread complete: wrote {} packs, registered {} segments",
                        packs_written, segments_registered
                    );
                }
                Ok(packs_written)
            })
        };

        // ================================================================
        // SEQUENTIAL PROCESSING: Process each contig segment-by-segment
        // ================================================================
        // This matches C++ AGC: assign seg_part_no DURING processing, accounting for splits immediately
        // NO pre-assigned part numbers = NO duplicate conflicts!

        while let Some((sample_name, contig_name, sequence)) = iterator.next_contig()? {
            if sequence.is_empty() {
                continue;
            }

            // Log sample changes
            if current_sample != sample_name {
                if !current_sample.is_empty() && self.config.verbosity > 0 {
                    println!("Finished processing sample: {}", current_sample);
                }
                current_sample = sample_name.clone();
                if self.config.verbosity > 0 {
                    println!("Processing sample: {}", sample_name);
                }
            }

            total_contigs += 1;
            total_bases += sequence.len();

            // Segment contig at splitters
            let segments = split_at_splitters_with_size(
                &sequence,
                splitters,
                self.config.kmer_length as usize,
                self.config.segment_size as usize,
            );

            // Process each segment IMMEDIATELY (sequential!)
            let mut seg_part_no = 0;

            for segment in segments {
                let kmer_front = segment.front_kmer;
                let kmer_back = segment.back_kmer;
                let segment_data = segment.data;

                // Create normalized key
                let (key, _) = SegmentGroupKey::new_normalized(kmer_front, kmer_back);

                // **KEY FIX**: Try to SPLIT FIRST (before checking if new/known)
                // C++ AGC checks for splits on ALL segments, not just known ones
                let can_split = kmer_front != MISSING_KMER
                    && kmer_back != MISSING_KMER
                    && !self.config.concatenated_genomes;

                let mut did_split = false;

                if can_split {
                    // Try to find shared k-mer
                    if let Some(middle_kmer) =
                        self.find_shared_kmer(kmer_front, kmer_back, &group_terminators)
                    {
                        // Check if both target groups exist
                        let key1 = SegmentGroupKey::new_normalized(kmer_front, middle_kmer).0;
                        let key2 = SegmentGroupKey::new_normalized(middle_kmer, kmer_back).0;

                        if std::env::var("RAGC_DEBUG_SPLIT").is_ok() {
                            eprintln!("SPLIT_ATTEMPT: front={} back={} middle={} key1_exists={} key2_exists={}",
                                kmer_front, kmer_back, middle_kmer,
                                groups.contains_key(&key1), groups.contains_key(&key2));
                        }

                        if groups.contains_key(&key1) && groups.contains_key(&key2) {
                            // Calculate split position
                            if let Some(split_pos) = Self::calculate_split_position(
                                kmer_front,
                                kmer_back,
                                middle_kmer,
                                &segment_data,
                                &groups,
                                &self.config,
                            ) {
                                // SPLIT THE SEGMENT!
                                let kmer_len = self.config.kmer_length as usize;
                                let seg2_start_pos = split_pos.saturating_sub(kmer_len / 2);

                                // Ensure seg1_end doesn't exceed segment bounds
                                let seg1_end =
                                    std::cmp::min(seg2_start_pos + kmer_len, segment_data.len());
                                let seg1_data = segment_data[..seg1_end].to_vec();
                                let seg2_data = segment_data[seg2_start_pos..].to_vec();

                                // **C++ AGC compatibility**: Only split if BOTH segments are >= kmer_length
                                // If either segment would be too short, skip the split
                                if seg1_data.len() < kmer_len || seg2_data.len() < kmer_len {
                                    if std::env::var("RAGC_DEBUG_SPLIT").is_ok() {
                                        eprintln!(
                                            "SPLIT_REJECTED: seg1_len={} seg2_len={} < kmer_len={} - using full segment instead",
                                            seg1_data.len(),
                                            seg2_data.len(),
                                            kmer_len
                                        );
                                    }
                                    // Don't split - continue to normal segment processing below
                                } else {
                                // Valid split - both segments are long enough

                                // Create first split segment (seg_part_no)
                                let seg1_info = Self::prepare_segment_info(
                                    &self.config,
                                    &sample_name,
                                    &contig_name,
                                    seg_part_no,
                                    seg1_data,
                                    Some(kmer_front),
                                    Some(middle_kmer),
                                )?;

                                // Create second split segment (seg_part_no + 1)
                                let seg2_info = Self::prepare_segment_info(
                                    &self.config,
                                    &sample_name,
                                    &contig_name,
                                    seg_part_no + 1,
                                    seg2_data,
                                    Some(middle_kmer),
                                    Some(kmer_back),
                                )?;

                                // Add both segments to their respective groups
                                for (seg_key, seg_info_prepared) in [
                                    (seg1_info.key.clone(), seg1_info),
                                    (seg2_info.key.clone(), seg2_info),
                                ] {
                                    let group_entry =
                                        groups.entry(seg_key.clone()).or_insert_with(|| {
                                            let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                                            known_groups.insert(seg_key.clone(), gid);

                                            // Add terminators
                                            if seg_key.kmer_front != MISSING_KMER
                                                && seg_key.kmer_back != MISSING_KMER
                                            {
                                                group_terminators
                                                    .entry(seg_key.kmer_front)
                                                    .or_default()
                                                    .push(seg_key.kmer_back);
                                                if seg_key.kmer_front != seg_key.kmer_back {
                                                    group_terminators
                                                        .entry(seg_key.kmer_back)
                                                        .or_default()
                                                        .push(seg_key.kmer_front);
                                                }
                                            }

                                            let stream_id = gid as usize;
                                            let ref_stream_id = if gid >= 16 {
                                                Some(100000 + gid as usize)
                                            } else {
                                                None
                                            };
                                            (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
                                        });

                                    if let Some(pack) = group_entry
                                        .1
                                        .add_segment(seg_info_prepared.segment, &self.config)?
                                    {
                                        let compressed_pack = match pack {
                                            PackToWrite::Compressed(cp) => cp,
                                            PackToWrite::Uncompressed(up) => {
                                                let compressed_data = compress_segment_configured(
                                                    &up.uncompressed_data,
                                                    self.config.compression_level,
                                                )?;
                                                let (final_data, uncompressed_size) =
                                                    if compressed_data.len() + 1
                                                        < up.uncompressed_data.len()
                                                    {
                                                        let mut data_with_marker = compressed_data;
                                                        data_with_marker.push(0);
                                                        (
                                                            data_with_marker,
                                                            up.uncompressed_data.len() as u64,
                                                        )
                                                    } else {
                                                        (up.uncompressed_data, 0)
                                                    };
                                                CompressedPack {
                                                    group_id: up.group_id,
                                                    stream_id: up.stream_id,
                                                    compressed_data: final_data,
                                                    uncompressed_size,
                                                    segments: up.segments,
                                                }
                                            }
                                        };
                                        pack_tx.send(compressed_pack)?;
                                    }
                                }

                                seg_part_no += 2; // Increment by 2 for split
                                total_segments += 2;
                                did_split = true;
                                } // End of valid split (both segments >= kmer_length)
                            }
                        }
                    }
                }

                // Only process as normal segment if we didn't split
                if !did_split {
                    // Classify: NEW or KNOWN?
                    let group_id_opt = known_groups.get(&key).copied();
                    let is_new = group_id_opt.is_none();

                    if is_new {
                        // NEW SEGMENT: Create group (but should we split first?)
                        if std::env::var("RAGC_DEBUG_SPLIT").is_ok()
                            && kmer_front != MISSING_KMER
                            && kmer_back != MISSING_KMER
                        {
                            eprintln!(
                                "NEW_SEGMENT: front={} back={} size={}",
                                kmer_front,
                                kmer_back,
                                segment_data.len()
                            );
                        }

                        let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                        known_groups.insert(key.clone(), gid);

                        // Optional debug logging for compatibility testing
                        if std::env::var("RAGC_DEBUG_KMER_PAIRS").is_ok() {
                            eprintln!(
                                "GROUP_KMER: group_id={} front={} back={}",
                                gid, key.kmer_front, key.kmer_back
                            );
                        }

                        // Add terminators (use NORMALIZED key values to match C++ AGC!)
                        if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                            group_terminators
                                .entry(key.kmer_front)
                                .or_default()
                                .push(key.kmer_back);
                            if key.kmer_front != key.kmer_back {
                                group_terminators
                                    .entry(key.kmer_back)
                                    .or_default()
                                    .push(key.kmer_front);
                            }
                        }

                        // Prepare segment
                        let seg_info = Self::prepare_segment_info(
                            &self.config,
                            &sample_name,
                            &contig_name,
                            seg_part_no,
                            segment_data,
                            Some(kmer_front),
                            Some(kmer_back),
                        )?;

                        // Create group and add segment
                        let stream_id = gid as usize;
                        let ref_stream_id = if gid >= 16 {
                            Some(100000 + gid as usize)
                        } else {
                            None
                        };
                        let mut group_writer = GroupWriter::new(gid, stream_id, ref_stream_id);

                        if let Some(pack) =
                            group_writer.add_segment(seg_info.segment, &self.config)?
                        {
                            let compressed_pack = match pack {
                                PackToWrite::Compressed(cp) => cp,
                                PackToWrite::Uncompressed(up) => {
                                    let compressed_data = compress_segment_configured(
                                        &up.uncompressed_data,
                                        self.config.compression_level,
                                    )?;
                                    let (final_data, uncompressed_size) =
                                        if compressed_data.len() + 1 < up.uncompressed_data.len() {
                                            let mut data_with_marker = compressed_data;
                                            data_with_marker.push(0);
                                            (data_with_marker, up.uncompressed_data.len() as u64)
                                        } else {
                                            (up.uncompressed_data, 0)
                                        };
                                    CompressedPack {
                                        group_id: up.group_id,
                                        stream_id: up.stream_id,
                                        compressed_data: final_data,
                                        uncompressed_size,
                                        segments: up.segments,
                                    }
                                }
                            };
                            pack_tx.send(compressed_pack)?;
                        }

                        groups.insert(key, (gid, group_writer));

                        seg_part_no += 1; // Increment by 1 for normal segment
                        total_segments += 1;
                    } else {
                        // KNOWN SEGMENT: Add to existing group (split was already attempted above)
                        let seg_info = Self::prepare_segment_info(
                            &self.config,
                            &sample_name,
                            &contig_name,
                            seg_part_no,
                            segment_data,
                            Some(kmer_front),
                            Some(kmer_back),
                        )?;

                        let group_entry = groups.entry(key.clone()).or_insert_with(|| {
                            let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                            known_groups.insert(key.clone(), gid);

                            // Add terminators
                            if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                                group_terminators
                                    .entry(key.kmer_front)
                                    .or_default()
                                    .push(key.kmer_back);
                                if key.kmer_front != key.kmer_back {
                                    group_terminators
                                        .entry(key.kmer_back)
                                        .or_default()
                                        .push(key.kmer_front);
                                }
                            }

                            let stream_id = gid as usize;
                            let ref_stream_id = if gid >= 16 {
                                Some(100000 + gid as usize)
                            } else {
                                None
                            };
                            (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
                        });

                        if let Some(pack) =
                            group_entry.1.add_segment(seg_info.segment, &self.config)?
                        {
                            let compressed_pack = match pack {
                                PackToWrite::Compressed(cp) => cp,
                                PackToWrite::Uncompressed(up) => {
                                    let compressed_data = compress_segment_configured(
                                        &up.uncompressed_data,
                                        self.config.compression_level,
                                    )?;
                                    let (final_data, uncompressed_size) =
                                        if compressed_data.len() + 1 < up.uncompressed_data.len() {
                                            let mut data_with_marker = compressed_data;
                                            data_with_marker.push(0);
                                            (data_with_marker, up.uncompressed_data.len() as u64)
                                        } else {
                                            (up.uncompressed_data, 0)
                                        };
                                    CompressedPack {
                                        group_id: up.group_id,
                                        stream_id: up.stream_id,
                                        compressed_data: final_data,
                                        uncompressed_size,
                                        segments: up.segments,
                                    }
                                }
                            };
                            pack_tx.send(compressed_pack)?;
                        }

                        seg_part_no += 1; // Increment by 1
                        total_segments += 1;
                    }
                }
            }
        }

        // ================================================================
        // FLUSH REMAINING GROUPS
        // ================================================================
        if self.config.verbosity > 0 {
            println!("Flushing {} remaining groups...", groups.len());
        }

        for (_, (gid, mut group_writer)) in groups {
            if let Some(uncompressed_pack) = group_writer.prepare_pack(&self.config)? {
                if self.config.verbosity > 1 && (gid == 11813 || gid < 100) {
                    println!(
                        "  Flushing delta pack for group_id={}, {} segments buffered",
                        gid,
                        uncompressed_pack.segments.len()
                    );
                }
                let compressed_data = compress_segment_configured(
                    &uncompressed_pack.uncompressed_data,
                    self.config.compression_level,
                )?;

                let (final_data, uncompressed_size) =
                    if compressed_data.len() + 1 < uncompressed_pack.uncompressed_data.len() {
                        let mut data_with_marker = compressed_data;
                        data_with_marker.push(0);
                        (
                            data_with_marker,
                            uncompressed_pack.uncompressed_data.len() as u64,
                        )
                    } else {
                        (uncompressed_pack.uncompressed_data, 0)
                    };

                let compressed_pack = CompressedPack {
                    group_id: uncompressed_pack.group_id,
                    stream_id: uncompressed_pack.stream_id,
                    compressed_data: final_data,
                    uncompressed_size,
                    segments: uncompressed_pack.segments,
                };

                pack_tx.send(compressed_pack)?;
            }
        }

        // ================================================================
        // SIGNAL WRITER AND WAIT
        // ================================================================
        drop(pack_tx);

        let packs_written = writer_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

        if self.config.verbosity > 0 {
            println!("Writer complete: {} packs written", packs_written);
        }

        // Update next_group_id for future calls
        self.next_group_id = next_group_id.load(Ordering::SeqCst);

        // Update total_segments for finalize() (critical for collection metadata)
        self.total_segments = total_segments;

        // Restore archive and collection
        self.archive = Arc::try_unwrap(archive)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc for archive"))?
            .into_inner()
            .unwrap();
        self.collection = Arc::try_unwrap(collection)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc for collection"))?
            .into_inner()
            .unwrap();

        if self.config.verbosity > 0 {
            println!(
                "Compression complete: {} contigs, {} bases, {} segments",
                total_contigs, total_bases, total_segments
            );
        }

        Ok(())
    }

    /// Worker-based contig processing (matching C++ AGC architecture)
    ///
    /// This replaces single-threaded segmentation with parallel workers:
    /// - Main thread: Streams contigs into BoundedPriorityQueue
    /// - Worker threads: Pull contigs, segment, compress, send to writer
    /// - Writer thread: Receives compressed packs, writes to archive
    ///
    /// C++ AGC equivalent: start_compressing_threads() + AddSampleFiles()
    fn add_segments_with_workers(
        &mut self,
        mut iterator: Box<dyn ContigIterator>,
        splitters: &HashSet<u64>,
    ) -> Result<()> {
        let num_workers = self.config.num_threads;

        if self.config.verbosity > 0 {
            println!("Starting C++ AGC-style worker-based compression ({} workers)...", num_workers);
        }

        // ================================================================
        // SETUP: Shared state (thread-safe)
        // ================================================================

        // Queue capacity: 2GB or 192MB/thread (matching C++ AGC line 2150)
        let queue_capacity = std::cmp::max(
            2usize << 30,                          // 2GB
            num_workers * (192usize << 20)         // 192MB/thread
        );

        // Priority queue for contigs (matching C++ AGC CBoundedPQueue)
        let contig_queue = BoundedPriorityQueue::new(1, queue_capacity);

        // Shared state for workers (thread-safe with DashMap)
        let known_groups = Arc::new(DashMap::<SegmentGroupKey, u32>::new());
        let groups = Arc::new(DashMap::<SegmentGroupKey, (u32, Mutex<GroupWriter>)>::new());
        let group_terminators = Arc::new(DashMap::<u64, Vec<u64>>::new());
        let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));

        // Pending segments storage (C++ AGC-style batch registration)
        // Workers add segments here during Normal processing
        // Thread 0 registers them deterministically at Registration barriers
        let pending_segments = Arc::new(Mutex::new(Vec::<PendingSegment>::new()));

        // Segments ready for compression (after registration)
        // Map from group_id to Vec<PendingSegment>
        let segments_to_compress = Arc::new(DashMap::<u32, Vec<PendingSegment>>::new());

        // Tracking counters
        let total_segments = Arc::new(AtomicUsize::new(0));
        let total_bases = Arc::new(AtomicUsize::new(0));
        let total_contigs = Arc::new(AtomicUsize::new(0));

        // Channel for compressed packs (writer thread)
        let (pack_tx, pack_rx) = bounded::<CompressedPack>(100);

        // ================================================================
        // WRITER THREAD (same as existing code)
        // ================================================================
        let archive = Arc::new(Mutex::new(std::mem::replace(
            &mut self.archive,
            Archive::new_writer(),
        )));
        let collection = Arc::new(Mutex::new(std::mem::replace(
            &mut self.collection,
            CollectionV3::new(),
        )));

        let writer_handle = {
            let archive = archive.clone();
            let collection = collection.clone();
            let config = self.config.clone();

            thread::spawn(move || -> Result<usize> {
                let mut stream_id_map: HashMap<usize, usize> = HashMap::new();
                let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
                let mut packs_written = 0;
                let mut segments_registered = 0;

                while let Ok(pack) = pack_rx.recv() {
                    // Wrap processing in error handler to capture exact failure point
                    if let Err(e) = (|| -> Result<()> {
                        let mut archive_lock = archive.lock().unwrap();
                        let mut collection_lock = collection.lock().unwrap();

                        // Register stream on-demand
                        let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                            id
                        } else {
                            let stream_name = if pack.stream_id >= 100000 {
                                stream_ref_name(archive_version, (pack.stream_id - 100000) as u32)
                            } else {
                                stream_delta_name(archive_version, pack.group_id)
                            };
                            let id = archive_lock.register_stream(&stream_name);
                            stream_id_map.insert(pack.stream_id, id);
                            id
                        };

                        // Register segments in collection
                        for seg_meta in &pack.segments {
                            // Register sample/contig first (required before add_segment_placed)
                            if let Err(e) = collection_lock.register_sample_contig(
                                &seg_meta.sample_name,
                                &seg_meta.contig_name,
                            ) {
                                eprintln!("ERROR in writer thread: register_sample_contig failed");
                                eprintln!("  Sample: {}", seg_meta.sample_name);
                                eprintln!("  Contig: {}", seg_meta.contig_name);
                                eprintln!("  Error: {}", e);
                                return Err(e);
                            }

                            if let Err(e) = collection_lock.add_segment_placed(
                                &seg_meta.sample_name,
                                &seg_meta.contig_name,
                                seg_meta.seg_part_no,
                                pack.group_id,
                                seg_meta.in_group_id,
                                seg_meta.is_rev_comp,
                                seg_meta.data_len,
                            ) {
                                eprintln!("ERROR in writer thread: add_segment_placed failed");
                                eprintln!("  Sample: {}", seg_meta.sample_name);
                                eprintln!("  Contig: {}", seg_meta.contig_name);
                                eprintln!("  seg_part_no: {}", seg_meta.seg_part_no);
                                eprintln!("  group_id: {}", pack.group_id);
                                eprintln!("  in_group_id: {}", seg_meta.in_group_id);
                                eprintln!("  is_rev_comp: {}", seg_meta.is_rev_comp);
                                eprintln!("  data_len: {}", seg_meta.data_len);
                                eprintln!("  Error: {}", e);
                                return Err(e);
                            }
                            segments_registered += 1;
                        }

                        // Write compressed data to archive
                        if let Err(e) = archive_lock.add_part(
                            actual_stream_id,
                            &pack.compressed_data,
                            pack.uncompressed_size,
                        ) {
                            eprintln!("ERROR in writer thread: add_part failed");
                            eprintln!("  stream_id: {}", actual_stream_id);
                            eprintln!("  group_id: {}", pack.group_id);
                            eprintln!("  compressed_size: {}", pack.compressed_data.len());
                            eprintln!("  uncompressed_size: {}", pack.uncompressed_size);
                            eprintln!("  Error: {}", e);
                            return Err(e);
                        }

                        Ok(())
                    })() {
                        eprintln!("Writer thread exiting due to error processing pack for group {}", pack.group_id);
                        return Err(e);
                    }

                    packs_written += 1;
                }

                if config.verbosity > 0 {
                    println!("Writer thread complete: wrote {} packs, registered {} segments",
                        packs_written, segments_registered);
                }
                Ok(packs_written)
            })
        };

        // ================================================================
        // WORKER THREADS (matching C++ AGC architecture)
        // ================================================================

        // Create barrier for synchronization (C++ AGC-style batch processing)
        let barrier = Arc::new(Barrier::new(num_workers));

        let mut worker_handles = Vec::new();

        for worker_id in 0..num_workers {
            let queue = contig_queue.clone();
            let barrier = barrier.clone();
            let known_groups = known_groups.clone();
            let groups = groups.clone();
            let group_terminators = group_terminators.clone();
            let next_group_id = next_group_id.clone();
            let pending_segments = pending_segments.clone();
            let segments_to_compress = segments_to_compress.clone();
            let pack_tx = pack_tx.clone();
            let splitters = splitters.clone();
            let config = self.config.clone();
            let total_segments_counter = total_segments.clone();
            let total_bases_counter = total_bases.clone();

            let handle = thread::spawn(move || -> Result<()> {
                // Worker loop (matching C++ AGC lines 1108-1275)
                loop {
                    let (result, task_opt) = queue.pop_large();

                    match result {
                        PopResult::Empty => continue,
                        PopResult::Completed => break,
                        PopResult::Normal => {
                            let task: ContigTask = task_opt.unwrap();

                            // Handle task based on type
                            match task {
                                ContigTask::Normal { sample_name, contig_name, sequence, seq_num } => {
                                    // Segment contig at splitters
                                    let segments = split_at_splitters_with_size(
                                        &sequence,
                                        &splitters,
                                        config.kmer_length as usize,
                                        config.segment_size as usize,
                                    );

                                    if config.verbosity > 1 {
                                        println!("Worker {} processing contig {} (seq_num={}, len={}, {} segments)",
                                            worker_id, contig_name, seq_num, sequence.len(), segments.len());
                                    }

                                    total_bases_counter.fetch_add(sequence.len(), Ordering::Relaxed);

                                    // Process each segment
                                    let mut seg_part_no = 0usize;

                                    for segment in segments {
                                        let kmer_front = segment.front_kmer;
                                        let kmer_back = segment.back_kmer;
                                        let segment_data = segment.data;

                                        // Create normalized key
                                        let (key, needs_rc) = SegmentGroupKey::new_normalized(kmer_front, kmer_back);

                                        // Store pending segment instead of creating group
                                        // Thread 0 will register groups deterministically at Registration barrier
                                        pending_segments.lock().unwrap().push(PendingSegment {
                                            key,
                                            sample_name: sample_name.clone(),
                                            contig_name: contig_name.clone(),
                                            seg_part_no,
                                            segment_data,
                                            kmer_front,
                                            kmer_back,
                                            seq_num,  // Preserve file order through sequence number
                                        });

                                        seg_part_no += 1;
                                        total_segments_counter.fetch_add(1, Ordering::Relaxed);
                                    }
                                }

                                ContigTask::Registration => {
                                    // C++ AGC-style batch synchronization (lines 123-130 in THREADING_REFACTOR.md)
                                    // All workers wait at barrier for deterministic group registration
                                    if config.verbosity > 1 {
                                        println!("Worker {} reached registration barrier", worker_id);
                                    }

                                    // Barrier 1: All workers stop processing
                                    barrier.wait();

                                    // Phase 3.1: Thread 0 registers groups deterministically (C++ AGC: agc_compressor.cpp:1123)
                                    if worker_id == 0 {
                                        // Take all pending segments (swap with empty vec to avoid long lock)
                                        let mut all_pending = Vec::new();
                                        std::mem::swap(&mut all_pending, &mut *pending_segments.lock().unwrap());

                                        // CRITICAL: Sort by seq_num to restore FILE ORDER!
                                        // Each ContigTask gets a seq_num when queued (preserves file order).
                                        // Workers process tasks in parallel → segments added in non-deterministic order.
                                        // Sorting by seq_num restores file order for deterministic group ID assignment.
                                        //
                                        // Note: Group IDs MUST match file order, not alphabetical order!
                                        // (Alphabetical: "chrIX" < "chrV", but file order may have chrV before chrIX)
                                        all_pending.sort_by_key(|seg| seg.seq_num);

                                        // Group segments by k-mer key in FILE ORDER (restored via seq_num sort!)
                                        // Groups created in order k-mer pairs first appear in sorted (= file order) list
                                        use std::collections::HashMap;
                                        let mut key_to_segments: HashMap<(u64, u64), Vec<PendingSegment>> = HashMap::new();
                                        let mut key_order: Vec<(u64, u64)> = Vec::new(); // Track first occurrence in sorted order

                                        for seg in all_pending {
                                            let key_tuple = (seg.key.kmer_front, seg.key.kmer_back);
                                            if !key_to_segments.contains_key(&key_tuple) {
                                                key_order.push(key_tuple); // First time seeing this k-mer pair in sorted list
                                            }
                                            key_to_segments.entry(key_tuple).or_insert_with(Vec::new).push(seg);
                                        }

                                        // Assign group IDs in SORTED ORDER (matching C++ AGC std::set iteration!)
                                        for (kmer_front, kmer_back) in key_order {
                                            let mut segments = key_to_segments.remove(&(kmer_front, kmer_back)).unwrap();

                                            // Sort segments within group by (sample, contig, seg_part_no) for correct ordering
                                            segments.sort_by(|a, b| {
                                                (&a.sample_name, &a.contig_name, a.seg_part_no)
                                                    .cmp(&(&b.sample_name, &b.contig_name, b.seg_part_no))
                                            });

                                            let key = SegmentGroupKey { kmer_front, kmer_back };

                                            // Check if group already exists, or create new one
                                            let gid = match known_groups.entry(key.clone()) {
                                                dashmap::mapref::entry::Entry::Occupied(e) => {
                                                    // Group already exists from previous registration
                                                    *e.get()
                                                }
                                                dashmap::mapref::entry::Entry::Vacant(e) => {
                                                    // NEW GROUP: Assign sequential ID
                                                    let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                                                    e.insert(gid);

                                                    // Add terminators (matches C++ AGC logic)
                                                    if kmer_front != MISSING_KMER && kmer_back != MISSING_KMER {
                                                        group_terminators
                                                            .entry(kmer_front)
                                                            .or_insert_with(Vec::new)
                                                            .push(kmer_back);
                                                        if kmer_front != kmer_back {
                                                            group_terminators
                                                                .entry(kmer_back)
                                                                .or_insert_with(Vec::new)
                                                                .push(kmer_front);
                                                        }
                                                    }

                                                    // Create GroupWriter for this group
                                                    let stream_id = gid as usize;
                                                    let ref_stream_id = if gid >= 16 {
                                                        Some(100000 + gid as usize)
                                                    } else {
                                                        None
                                                    };
                                                    groups.insert(
                                                        key.clone(),
                                                        (gid, Mutex::new(GroupWriter::new(gid, stream_id, ref_stream_id))),
                                                    );

                                                    if config.verbosity > 1 {
                                                        println!("  Created new group {} for k-mers ({}, {})", gid, kmer_front, kmer_back);
                                                    }

                                                    gid
                                                }
                                            };

                                            // Store segments for compression phase (Phase 4)
                                            segments_to_compress.insert(gid, segments);
                                        }

                                        if config.verbosity > 1 {
                                            println!("Thread 0 registration complete. Total groups: {}", known_groups.len());
                                        }
                                    }

                                    // Barrier 2: All workers can continue (after registration)
                                    barrier.wait();

                                    // Phase 4.1: All workers compress registered segments in parallel
                                    // (matches C++ AGC store_segments - agc_compressor.cpp:974-1054)
                                    if config.verbosity > 1 {
                                        println!("Worker {} starting compression phase", worker_id);
                                    }

                                    // Collect all group IDs and sort for deterministic processing
                                    let mut all_group_ids: Vec<u32> = segments_to_compress.iter()
                                        .map(|entry| *entry.key())
                                        .collect();
                                    all_group_ids.sort_unstable();

                                    // Process groups in sorted order (work-stealing with sorted list)
                                    while !all_group_ids.is_empty() {
                                        // Try to claim the next group (remove from front for FIFO order)
                                        let group_opt = if !all_group_ids.is_empty() {
                                            let gid = all_group_ids.remove(0);
                                            segments_to_compress.remove(&gid).map(|(k, v)| (k, v))
                                        } else {
                                            None
                                        };

                                        if let Some((gid, segments)) = group_opt {

                                            if config.verbosity > 1 {
                                                println!("Worker {} processing group {} ({} segments)",
                                                    worker_id, gid, segments.len());
                                            }

                                            // Process all segments in this group
                                            for segment in segments {
                                                // Prepare segment info
                                                let seg_info = Self::prepare_segment_info(
                                                    &config,
                                                    &segment.sample_name,
                                                    &segment.contig_name,
                                                    segment.seg_part_no,
                                                    segment.segment_data,
                                                    Some(segment.kmer_front),
                                                    Some(segment.kmer_back),
                                                )?;

                                                // Add to group (group exists from registration)
                                                if let Some(group_entry) = groups.get(&segment.key) {
                                                    let (_, group_writer_mutex) = &*group_entry;
                                                    let mut group_writer = group_writer_mutex.lock().unwrap();

                                                    if let Some(pack) = group_writer.add_segment(seg_info.segment, &config)? {
                                                        // Compress pack if needed
                                                        let compressed_pack = match pack {
                                                            PackToWrite::Compressed(cp) => cp,
                                                            PackToWrite::Uncompressed(up) => {
                                                                let compressed_data = compress_segment_configured(
                                                                    &up.uncompressed_data,
                                                                    config.compression_level,
                                                                )?;
                                                                let (final_data, uncompressed_size) =
                                                                    if compressed_data.len() + 1 < up.uncompressed_data.len() {
                                                                        let mut data_with_marker = compressed_data;
                                                                        data_with_marker.push(0);
                                                                        (data_with_marker, up.uncompressed_data.len() as u64)
                                                                    } else {
                                                                        (up.uncompressed_data, 0)
                                                                    };
                                                                CompressedPack {
                                                                    group_id: up.group_id,
                                                                    stream_id: up.stream_id,
                                                                    compressed_data: final_data,
                                                                    uncompressed_size,
                                                                    segments: up.segments,
                                                                }
                                                            }
                                                        };
                                                        pack_tx.send(compressed_pack)?;
                                                    }
                                                }
                                            }
                                        } else {
                                            // No more groups available
                                            break;
                                        }
                                    }

                                    // Flush any remaining partial packs in all groups (only worker 0)
                                    if worker_id == 0 {
                                        if config.verbosity > 1 {
                                            println!("Worker 0 flushing remaining partial packs from {} groups", groups.len());
                                        }

                                        // Sort by group_id for deterministic flushing (same order as registration)
                                        let mut groups_with_ids: Vec<(u32, SegmentGroupKey)> = groups.iter()
                                            .map(|entry| {
                                                let key = entry.key().clone();
                                                let (gid, _) = entry.value();
                                                (*gid, key)
                                            })
                                            .collect();
                                        groups_with_ids.sort_by_key(|(gid, _)| *gid);

                                        for (_gid, key) in groups_with_ids {
                                            if let Some(group_entry) = groups.get(&key) {
                                                let (_, group_writer_mutex) = &*group_entry;
                                                let mut group_writer = group_writer_mutex.lock().unwrap();

                                                // Flush any remaining segments
                                                if let Some(pack) = group_writer.flush_final(&config)? {
                                                    let compressed_pack = match pack {
                                                        PackToWrite::Compressed(cp) => cp,
                                                        PackToWrite::Uncompressed(up) => {
                                                            let compressed_data = compress_segment_configured(
                                                                &up.uncompressed_data,
                                                                config.compression_level,
                                                            )?;
                                                            let (final_data, uncompressed_size) =
                                                                if compressed_data.len() + 1 < up.uncompressed_data.len() {
                                                                    let mut data_with_marker = compressed_data;
                                                                    data_with_marker.push(0);
                                                                    (data_with_marker, up.uncompressed_data.len() as u64)
                                                                } else {
                                                                    (up.uncompressed_data, 0)
                                                                };
                                                            CompressedPack {
                                                                group_id: up.group_id,
                                                                stream_id: up.stream_id,
                                                                compressed_data: final_data,
                                                                uncompressed_size,
                                                                segments: up.segments,
                                                            }
                                                        }
                                                    };
                                                    pack_tx.send(compressed_pack)?;
                                                }
                                            }
                                        }
                                    }

                                    if config.verbosity > 1 {
                                        println!("Worker {} completed compression phase", worker_id);
                                    }
                                }
                            }
                        }
                    }
                }

                if config.verbosity > 1 {
                    println!("Worker {} complete", worker_id);
                }
                Ok(())
            });

            worker_handles.push(handle);
        }

        // ================================================================
        // MAIN THREAD: Stream contigs into queue
        // ================================================================
        let mut contig_priority = usize::MAX; // Higher priority = processed first

        while let Some((sample_name, contig_name, sequence)) = iterator.next_contig()? {
            if sequence.is_empty() {
                continue;
            }

            let contig_count = total_contigs.fetch_add(1, Ordering::Relaxed) + 1;
            let seq_num = contig_count as u64;  // Use contig count as sequence number (preserves file order)

            if self.config.verbosity > 1 {
                println!("Main thread: Feeding contig {} (seq_num={}, priority={})", contig_name, seq_num, contig_priority);
            }

            let task = ContigTask::Normal {
                sample_name,
                contig_name,
                sequence: sequence.clone(),
                seq_num,
            };

            let cost = sequence.len();
            contig_queue.emplace(task, contig_priority, cost);
            contig_priority = contig_priority.saturating_sub(1);

            // TODO: Periodically inject Registration task for batch processing (C++ AGC-style)
            // Disabled for now to simplify debugging
            // const BATCH_SIZE: usize = 1000;
            // if contig_count % BATCH_SIZE == 0 {
            //     contig_queue.emplace_many_no_cost(ContigTask::Registration, usize::MAX - 1, num_workers);
            // }
        }

        // Final registration phase for remaining segments
        // Use priority 0 so Registration is processed AFTER all contigs (contigs have usize::MAX down)
        let registration_priority = 0;
        if self.config.verbosity > 1 {
            println!("Main thread: All contigs fed. Injecting final Registration task (priority {})", registration_priority);
        }
        contig_queue.emplace_many_no_cost(ContigTask::Registration, registration_priority, num_workers);

        // Signal completion
        if self.config.verbosity > 1 {
            println!("Main thread: Marking queue as completed");
        }
        contig_queue.mark_completed();

        // Wait for workers
        for (idx, handle) in worker_handles.into_iter().enumerate() {
            match handle.join() {
                Ok(Ok(())) => {
                    if self.config.verbosity > 1 {
                        println!("Worker {} joined successfully", idx);
                    }
                }
                Ok(Err(e)) => {
                    eprintln!("ERROR: Worker {} failed: {}", idx, e);
                    return Err(e);
                }
                Err(_) => {
                    return Err(anyhow::anyhow!("Worker {} panicked", idx));
                }
            }
        }

        // Flush remaining groups (sorted for determinism - matches C++ AGC)
        let mut group_keys: Vec<_> = groups.iter().map(|e| e.key().clone()).collect();
        group_keys.sort_unstable();

        if self.config.verbosity > 0 {
            println!("Flushing {} remaining groups...", group_keys.len());
        }

        for key in group_keys {
            // Remove from DashMap to take ownership (no more ref holding)
            if let Some((_, (gid, group_writer_mutex))) = groups.remove(&key) {
                let mut group_writer = group_writer_mutex.lock().unwrap();

                if let Some(uncompressed_pack) = group_writer.prepare_pack(&self.config)? {
                let compressed_data = compress_segment_configured(
                    &uncompressed_pack.uncompressed_data,
                    self.config.compression_level,
                )?;
                let (final_data, uncompressed_size) =
                    if compressed_data.len() + 1 < uncompressed_pack.uncompressed_data.len() {
                        let mut data_with_marker = compressed_data;
                        data_with_marker.push(0);
                        (data_with_marker, uncompressed_pack.uncompressed_data.len() as u64)
                    } else {
                        (uncompressed_pack.uncompressed_data, 0)
                    };

                let compressed_pack = CompressedPack {
                    group_id: uncompressed_pack.group_id,
                    stream_id: uncompressed_pack.stream_id,
                    compressed_data: final_data,
                    uncompressed_size,
                    segments: uncompressed_pack.segments,
                };

                if let Err(e) = pack_tx.send(compressed_pack) {
                    eprintln!("ERROR: Failed to send pack for group {}: {}", gid, e);
                    eprintln!("This usually means the writer thread has exited early.");
                    eprintln!("Checking writer thread status...");

                    // Check if writer is still running
                    drop(pack_tx);
                    match writer_handle.join() {
                        Ok(Ok(count)) => {
                            eprintln!("Writer thread exited successfully with {} packs written", count);
                            eprintln!("But we still have groups to flush - this is unexpected!");
                        }
                        Ok(Err(writer_err)) => {
                            eprintln!("Writer thread failed with error: {}", writer_err);
                            return Err(writer_err);
                        }
                        Err(_) => {
                            eprintln!("Writer thread panicked");
                            return Err(anyhow::anyhow!("Writer thread panicked"));
                        }
                    }

                    return Err(anyhow::anyhow!("Channel send failed: {}", e));
                }
            }  // end if let Some(uncompressed_pack)
            }  // end if let Some(groups.remove)
        }  // end for key

        // Signal writer and wait
        drop(pack_tx);

        let packs_written = match writer_handle.join() {
            Ok(Ok(count)) => count,
            Ok(Err(e)) => {
                eprintln!("ERROR: Writer thread failed: {}", e);
                return Err(e);
            }
            Err(_) => {
                return Err(anyhow::anyhow!("Writer thread panicked"));
            }
        };

        if self.config.verbosity > 0 {
            println!("Writer complete: {} packs written", packs_written);
        }

        // Update state
        self.next_group_id = next_group_id.load(Ordering::SeqCst);
        self.total_segments = total_segments.load(Ordering::Relaxed);
        self.total_bases_processed = total_bases.load(Ordering::Relaxed);

        // Restore archive and collection
        self.archive = Arc::try_unwrap(archive)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc for archive"))?
            .into_inner()
            .unwrap();
        self.collection = Arc::try_unwrap(collection)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc for collection"))?
            .into_inner()
            .unwrap();

        if self.config.verbosity > 0 {
            println!(
                "Compression complete: {} contigs, {} bases, {} segments",
                total_contigs.load(Ordering::Relaxed),
                total_bases.load(Ordering::Relaxed),
                total_segments.load(Ordering::Relaxed)
            );
        }

        Ok(())
    }

    fn calculate_split_position(
        kmer_front: u64,
        kmer_back: u64,
        middle_kmer: u64,
        segment: &Contig,
        groups: &HashMap<SegmentGroupKey, (u32, GroupWriter)>,
        config: &StreamingCompressorConfig,
    ) -> Option<usize> {
        // Look up existing groups using normalized keys
        let (key1, _) = SegmentGroupKey::new_normalized(kmer_front, middle_kmer);
        let (key2, _) = SegmentGroupKey::new_normalized(middle_kmer, kmer_back);

        // Get reference segments from both groups
        let ref_seg1 = {
            let group_entry = groups.get(&key1)?;
            group_entry.1.reference.as_ref()?.clone()
        };

        let ref_seg2 = {
            let group_entry = groups.get(&key2)?;
            group_entry.1.reference.as_ref()?.clone()
        };

        // Compute reverse complement of segment (matching C++ exactly)
        let segment_rc: Vec<u8> = segment
            .iter()
            .rev()
            .map(|&b| match b {
                0 => 3, // A -> T
                1 => 2, // C -> G
                2 => 1, // G -> C
                3 => 0, // T -> A
                _ => b, // N stays N
            })
            .collect();

        // Prepare LZ encoders with reference segments
        let mut lz_diff1 = LZDiff::new(config.min_match_len);
        let mut lz_diff2 = LZDiff::new(config.min_match_len);

        lz_diff1.prepare(&ref_seg1);
        lz_diff2.prepare(&ref_seg2);

        // Compute cost vector for first part (C++ lines 1541-1551)
        let mut v_costs1 = if kmer_front < middle_kmer {
            lz_diff1.get_coding_cost_vector(segment, true)
        } else {
            let mut costs = lz_diff1.get_coding_cost_vector(&segment_rc, false);
            costs.reverse();
            costs
        };

        if v_costs1.is_empty() {
            return None;
        }

        // Apply forward partial_sum to v_costs1
        for i in 1..v_costs1.len() {
            v_costs1[i] += v_costs1[i - 1];
        }

        // Compute cost vector for second part (C++ lines 1553-1578)
        let mut v_costs2 = if middle_kmer < kmer_back {
            lz_diff2.get_coding_cost_vector(segment, false)
        } else {
            lz_diff2.get_coding_cost_vector(&segment_rc, true)
        };

        if v_costs2.is_empty() || v_costs1.len() != v_costs2.len() {
            return None;
        }

        // Apply partial_sum based on orientation
        if middle_kmer < kmer_back {
            // Reverse partial_sum
            for i in (0..v_costs2.len() - 1).rev() {
                v_costs2[i] += v_costs2[i + 1];
            }
        } else {
            // Forward partial_sum, then reverse
            for i in 1..v_costs2.len() {
                v_costs2[i] += v_costs2[i - 1];
            }
            v_costs2.reverse();
        }

        // Find position with minimum total cost (C++ lines 1606-1617)
        let mut best_sum = u32::MAX;
        let mut best_pos = 0;

        for i in 0..v_costs1.len() {
            let cs = v_costs1[i] + v_costs2[i];
            if cs < best_sum {
                best_sum = cs;
                best_pos = i;
            }
        }

        // Boundary validation (C++ lines 1619-1622)
        let kmer_length = config.kmer_length as usize;
        if best_pos < kmer_length + 1 {
            best_pos = 0;
        }
        if best_pos + kmer_length + 1 > v_costs1.len() {
            best_pos = v_costs1.len();
        }

        Some(best_pos)
    }

    fn prepare_segment_info(
        config: &StreamingCompressorConfig,
        sample_name: &str,
        contig_name: &str,
        seg_idx: usize,
        mut data: Contig,
        front_kmer: Option<u64>,
        back_kmer: Option<u64>,
    ) -> Result<PreparedSegment> {
        // Create normalized key
        let (key, needs_rc) = SegmentGroupKey::new_normalized(
            front_kmer.unwrap_or(MISSING_KMER),
            back_kmer.unwrap_or(MISSING_KMER),
        );

        // Apply RC if needed
        if needs_rc {
            data = data
                .iter()
                .rev()
                .map(|&b| match b {
                    0 => 3, // A -> T
                    1 => 2, // C -> G
                    2 => 1, // G -> C
                    3 => 0, // T -> A
                    _ => b, // N stays N
                })
                .collect();
        }

        Ok(PreparedSegment {
            key,
            segment: SegmentInfo {
                data,
                sample_name: sample_name.to_string(),
                contig_name: contig_name.to_string(),
                seg_part_no: seg_idx,
                is_rev_comp: needs_rc,
            },
        })
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
            let front_kmer_val = if front.is_full() {
                // C++ lines 1648-1654: When one k-mer is missing, orientation matters
                // If only this k-mer is present, check is_dir_oriented to determine position
                front.data()
            } else {
                MISSING_KMER
            };
            let front_is_dir = front.is_dir_oriented();

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
            let back_kmer_val = if back.is_full() {
                back.data()
            } else {
                MISSING_KMER
            };
            let back_is_dir = back.is_dir_oriented();

            // Match C++ logic for one-k-mer cases (lines 1326-1372)
            if front_kmer_val != MISSING_KMER && back_kmer_val == MISSING_KMER {
                // Only front k-mer present (C++ lines 1326-1346)
                // C++ checks is_dir_oriented() directly
                let result = if front_is_dir {
                    (front_kmer_val, MISSING_KMER) // dir-oriented: (kmer, MISSING)
                } else {
                    (MISSING_KMER, front_kmer_val) // RC-oriented: (MISSING, kmer)
                };

                // DEBUG: Log a few cases
                if self.total_segments < 20 {
                    eprintln!(
                        "RUST_ONE_KMER: FRONT kmer={}, is_dir={}, result=({}, {})",
                        front_kmer_val, front_is_dir, result.0, result.1
                    );
                }

                result
            } else if front_kmer_val == MISSING_KMER && back_kmer_val != MISSING_KMER {
                // Only back k-mer present (C++ lines 1348-1372)
                // CRITICAL: C++ calls kmer.swap_dir_rc() which swaps dir and rc
                // This INVERTS the is_dir_oriented() check!
                // So we need to use !back_is_dir
                let result = if !back_is_dir {
                    (back_kmer_val, MISSING_KMER) // After swap, dir-oriented: (kmer, MISSING)
                } else {
                    (MISSING_KMER, back_kmer_val) // After swap, RC-oriented: (MISSING, kmer)
                };

                // DEBUG: Log a few cases
                if self.total_segments < 20 {
                    eprintln!(
                        "RUST_ONE_KMER: BACK kmer={}, is_dir={}, !is_dir={}, result=({}, {})",
                        back_kmer_val, back_is_dir, !back_is_dir, result.0, result.1
                    );
                }

                result
            } else {
                // Both present or both missing - keep as is
                (front_kmer_val, back_kmer_val)
            }
        } else {
            (MISSING_KMER, MISSING_KMER)
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

    /// Find middle splitter and optimal split position (C++ AGC optimization)
    /// Returns (middle_kmer, split_position) if a better split is found
    /// This implements the algorithm from agc_compressor.cpp lines 1505-1623
    fn find_middle_splitter(
        &self,
        kmer_front: u64,
        kmer_back: u64,
        segment: &Contig,
    ) -> Option<(u64, usize)> {
        // Both k-mers must be present (not MISSING_KMER)
        if kmer_front == MISSING_KMER || kmer_back == MISSING_KMER {
            if self.total_segments < 20 {
                eprintln!("RUST_MS_DEBUG: MISSING_KMER check failed");
            }
            return None;
        }

        // Get terminators for both k-mers (use read lock)
        let front_terminators = self.group_terminators.get(&kmer_front).cloned();
        let back_terminators = self.group_terminators.get(&kmer_back).cloned();

        if self.total_segments < 20 {
            eprintln!(
                "RUST_MS_DEBUG: front_term={:?}, back_term={:?}",
                front_terminators.as_ref().map(|v| v.len()),
                back_terminators.as_ref().map(|v| v.len())
            );
        }

        let front_terminators = front_terminators?;
        let back_terminators = back_terminators?;

        if self.config.verbosity > 2 {
            eprintln!("    Looking for middle splitter: front={} has {} terminators, back={} has {} terminators",
                kmer_front, front_terminators.len(), kmer_back, back_terminators.len());
        }

        // Find intersection of terminators (shared k-mers)
        // C++ uses std::set_intersection on sorted vectors
        // Both vectors are sorted, so we can use binary search
        let shared_kmers: Vec<u64> = front_terminators
            .iter()
            .filter(|&&k| back_terminators.binary_search(&k).is_ok())
            .copied()
            .filter(|&k| k != MISSING_KMER)
            .collect();

        if self.config.verbosity > 2 && !shared_kmers.is_empty() {
            eprintln!(
                "    Found {} shared k-mers between terminators!",
                shared_kmers.len()
            );
        }

        if shared_kmers.is_empty() {
            return None;
        }

        // Take first shared k-mer as middle (C++ does this)
        let middle_kmer = shared_kmers[0];

        // Look up existing groups using normalized keys (matching C++ lines 1535-1536)
        // C++: auto segment_id1 = map_segments[minmax(kmer_front.data(), middle)];
        let (key1, _) = SegmentGroupKey::new_normalized(kmer_front, middle_kmer);
        let (key2, _) = SegmentGroupKey::new_normalized(middle_kmer, kmer_back);

        let group1 = self.segment_groups.get(&key1)?;
        let group2 = self.segment_groups.get(&key2)?;

        if group1.is_empty() || group2.is_empty() {
            if self.config.verbosity > 2 {
                eprintln!("    Groups exist but are empty");
            }
            return None;
        }

        // Get reference segments (first segment in each group)
        let ref_seg1 = &group1[0].data;
        let ref_seg2 = &group2[0].data;

        if self.config.verbosity > 2 {
            eprintln!(
                "    Found reference groups: ({}, {}) with {} segs, ({}, {}) with {} segs",
                kmer_front,
                middle_kmer,
                group1.len(),
                middle_kmer,
                kmer_back,
                group2.len()
            );
            eprintln!(
                "    Reference segment sizes: {} and {} bytes",
                ref_seg1.len(),
                ref_seg2.len()
            );
        }

        // Compute reverse complement of segment (matching C++ exactly)
        let segment_rc: Vec<u8> = segment
            .iter()
            .rev()
            .map(|&b| {
                match b {
                    0 => 3, // A -> T
                    1 => 2, // C -> G
                    2 => 1, // G -> C
                    3 => 0, // T -> A
                    _ => b, // N stays N
                }
            })
            .collect();

        // Prepare LZ encoders with reference segments
        let mut lz_diff1 = LZDiff::new(self.config.min_match_len);
        let mut lz_diff2 = LZDiff::new(self.config.min_match_len);

        lz_diff1.prepare(ref_seg1);
        lz_diff2.prepare(ref_seg2);

        // Compute cost vector for first part (matching C++ seg1_run logic)
        // C++ lines 1541-1551
        let mut v_costs1 = if kmer_front < middle_kmer {
            // Use forward segment with prefix_costs=true
            lz_diff1.get_coding_cost_vector(segment, true)
        } else {
            // Use RC segment with prefix_costs=false, then reverse
            let mut costs = lz_diff1.get_coding_cost_vector(&segment_rc, false);
            costs.reverse();
            costs
        };

        if v_costs1.is_empty() {
            if self.config.verbosity > 2 {
                eprintln!("    v_costs1 is empty after computation");
            }
            return None;
        }

        // Apply forward partial_sum to v_costs1
        for i in 1..v_costs1.len() {
            v_costs1[i] += v_costs1[i - 1];
        }

        // Compute cost vector for second part (matching C++ seg2_run logic)
        // C++ lines 1553-1578
        let mut v_costs2 = if middle_kmer < kmer_back {
            // Use forward segment with prefix_costs=false

            lz_diff2.get_coding_cost_vector(segment, false)
        } else {
            // Use RC segment with prefix_costs=true

            lz_diff2.get_coding_cost_vector(&segment_rc, true)
        };

        if v_costs2.is_empty() || v_costs1.len() != v_costs2.len() {
            if self.config.verbosity > 2 {
                eprintln!(
                    "    Cost vector issue: v_costs2 len={}, v_costs1 len={}",
                    v_costs2.len(),
                    v_costs1.len()
                );
            }
            return None;
        }

        // Apply partial_sum based on orientation
        if middle_kmer < kmer_back {
            // Reverse partial_sum: partial_sum(v_costs2.rbegin(), v_costs2.rend(), v_costs2.rbegin())
            // This means: v_costs2[i] = sum from i to end
            for i in (0..v_costs2.len() - 1).rev() {
                v_costs2[i] += v_costs2[i + 1];
            }
        } else {
            // Forward partial_sum, then reverse
            for i in 1..v_costs2.len() {
                v_costs2[i] += v_costs2[i - 1];
            }
            v_costs2.reverse();
        }

        // Find position with minimum total cost (C++ lines 1606-1617)
        let mut best_sum = u32::MAX;
        let mut best_pos = 0;

        for i in 0..v_costs1.len() {
            let cs = v_costs1[i] + v_costs2[i];
            if cs < best_sum {
                best_sum = cs;
                best_pos = i;
            }
        }

        // Boundary validation (C++ lines 1619-1622)
        // Don't split too close to edges - need at least k+1 bases
        let kmer_length = self.config.kmer_length as usize;
        if best_pos < kmer_length + 1 {
            best_pos = 0;
        }
        if best_pos + kmer_length + 1 > v_costs1.len() {
            best_pos = v_costs1.len();
        }

        if self.config.verbosity > 2 {
            eprintln!(
                "    Best split position: {best_pos} (after boundary check) with cost {best_sum}"
            );
        }

        Some((middle_kmer, best_pos))
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
        // CRITICAL FIX: When one k-mer is MISSING, we need to check orientation
        // to determine which position to put the k-mer in (matching C++ lines 1648-1654)
        let (final_front, final_back) = if kmer_front != MISSING_KMER && kmer_back == MISSING_KMER {
            // Only front k-mer present - extract k-mer object to check orientation
            let k = self.config.kmer_length as usize;
            if segment.len() >= k {
                let mut front = Kmer::new(self.config.kmer_length, KmerMode::Canonical);
                for i in 0..k {
                    if segment[i] > 3 {
                        front.reset();
                        break;
                    }
                    front.insert(segment[i] as u64);
                }

                if front.is_full() && front.data() == kmer_front {
                    // C++ logic (lines 1648-1654): check is_dir_oriented()
                    if front.is_dir_oriented() {
                        (kmer_front, MISSING_KMER) // dir-oriented: (kmer, MISSING)
                    } else {
                        (MISSING_KMER, kmer_front) // RC-oriented: (MISSING, kmer)
                    }
                } else {
                    (kmer_front, kmer_back) // Fallback if extraction failed
                }
            } else {
                (kmer_front, kmer_back)
            }
        } else if kmer_front == MISSING_KMER && kmer_back != MISSING_KMER {
            // Only back k-mer present - extract k-mer object to check orientation
            let k = self.config.kmer_length as usize;
            if segment.len() >= k {
                let mut back = Kmer::new(self.config.kmer_length, KmerMode::Canonical);
                let start = segment.len() - k;
                for i in 0..k {
                    if segment[start + i] > 3 {
                        back.reset();
                        break;
                    }
                    back.insert(segment[start + i] as u64);
                }

                if back.is_full() && back.data() == kmer_back {
                    // CRITICAL: C++ calls swap_dir_rc() for back k-mers (line 1363)
                    // which inverts is_dir_oriented(). So use !is_dir_oriented()
                    if !back.is_dir_oriented() {
                        (kmer_back, MISSING_KMER) // After swap, dir-oriented: (kmer, MISSING)
                    } else {
                        (MISSING_KMER, kmer_back) // After swap, RC-oriented: (MISSING, kmer)
                    }
                } else {
                    (kmer_front, kmer_back) // Fallback if extraction failed
                }
            } else {
                (kmer_front, kmer_back)
            }
        } else {
            // Both present or both missing - keep as-is
            (kmer_front, kmer_back)
        };

        // Normalize key (matching C++ pk = minmax(kmer1, kmer2) at line 1310)
        let (key, needs_rc) = SegmentGroupKey::new_normalized(final_front, final_back);

        // If normalization flipped the k-mers, we need to store the RC
        let final_is_rev_comp = is_rev_comp ^ needs_rc;

        // Reverse complement the segment if needed to match normalized orientation
        let final_segment = if needs_rc {
            segment
                .iter()
                .rev()
                .map(|&b| match b {
                    0 => 3, // A -> T
                    1 => 2, // C -> G
                    2 => 1, // G -> C
                    3 => 0, // T -> A
                    _ => b, // N stays N
                })
                .collect()
        } else {
            segment.clone()
        };

        // C++ AGC optimization: missing middle splitter
        // If this group doesn't exist, try to split the segment using existing groups
        // CRITICAL: Must pass NORMALIZED k-mers since terminators are stored with normalized keys!
        if !self.segment_groups.contains_key(&key) {
            // DEBUG: Log when we try middle splitter
            if self.total_segments < 10 && kmer_front != MISSING_KMER && kmer_back != MISSING_KMER {
                eprintln!("RUST_MS_TRY: segment={}, orig_kmers=({}, {}), normalized=({}, {}), group_exists={}",
                    self.total_segments, kmer_front, kmer_back, key.kmer_front, key.kmer_back, self.segment_groups.contains_key(&key));
            }

            // Use NORMALIZED k-mers for lookup since that's how terminators are stored
            if let Some((middle_kmer, split_pos)) =
                self.find_middle_splitter(key.kmer_front, key.kmer_back, &final_segment)
            {
                // DEBUG: Log when middle splitter succeeds
                if self.total_segments < 10 {
                    eprintln!("RUST_MS_SUCCESS: middle_kmer={middle_kmer}, split_pos={split_pos}");
                }

                // Matching C++ lines 1400-1444: handle zero-size cases
                let left_size = split_pos;
                let right_size = final_segment.len() - split_pos;

                if left_size == 0 {
                    // C++ lines 1403-1410: left side is empty, just update the key
                    // Don't split - just add with the new key (middle_kmer, kmer_back)
                    if self.config.verbosity > 1 {
                        eprintln!("  Middle splitter: left_size=0, updating key from ({kmer_front}, {kmer_back}) to ({middle_kmer}, {kmer_back})");
                    }

                    // Add with updated key
                    return self.add_segment_with_kmers(
                        sample_name,
                        contig_name,
                        seg_part_no,
                        final_segment,
                        middle_kmer,
                        kmer_back,
                        final_is_rev_comp,
                    );
                } else if right_size == 0 {
                    // C++ lines 1411-1418: right side is empty, just update the key
                    // Don't split - just add with the new key (kmer_front, middle_kmer)
                    if self.config.verbosity > 1 {
                        eprintln!("  Middle splitter: right_size=0, updating key from ({kmer_front}, {kmer_back}) to ({kmer_front}, {middle_kmer})");
                    }

                    // Add with updated key
                    return self.add_segment_with_kmers(
                        sample_name,
                        contig_name,
                        seg_part_no,
                        final_segment,
                        kmer_front,
                        middle_kmer,
                        final_is_rev_comp,
                    );
                } else {
                    // Both sizes > 0: Split with overlap (C++ lines 1419-1444)
                    if self.config.verbosity > 1 {
                        eprintln!(
                            "  Middle splitter optimization: splitting segment at position {split_pos}"
                        );
                        eprintln!("    Original group: ({kmer_front}, {kmer_back})");
                        eprintln!(
                            "    Split into: ({kmer_front}, {middle_kmer}) and ({middle_kmer}, {kmer_back})"
                        );
                    }

                    // C++ line 1425: Create overlap of kmer_length/2 between segments
                    let kmer_len = self.config.kmer_length as usize;
                    let seg2_start_pos = split_pos.saturating_sub(kmer_len / 2);

                    // C++ lines 1426-1428:
                    // segment2 = [seg2_start_pos, end)
                    // segment1 = [0, seg2_start_pos + kmer_length)
                    // Ensure seg1_end doesn't exceed segment bounds
                    let seg1_end = std::cmp::min(seg2_start_pos + kmer_len, final_segment.len());
                    let seg1 = final_segment[..seg1_end].to_vec();
                    let seg2 = final_segment[seg2_start_pos..].to_vec();

                    if self.config.verbosity > 1 {
                        eprintln!(
                            "    Overlap: seg1_len={}, seg2_len={}, overlap_size={}",
                            seg1.len(),
                            seg2.len(),
                            seg1.len() + seg2.len() - final_segment.len()
                        );
                    }

                    // Recursively add both segments (they will be normalized inside)
                    // Note: This could trigger more middle splitter searches
                    self.add_segment_with_kmers(
                        sample_name,
                        contig_name,
                        seg_part_no,
                        seg1,
                        kmer_front,
                        middle_kmer,
                        final_is_rev_comp,
                    )?;

                    self.add_segment_with_kmers(
                        sample_name,
                        contig_name,
                        seg_part_no,
                        seg2,
                        middle_kmer,
                        kmer_back,
                        final_is_rev_comp,
                    )?;

                    return Ok(());
                }
            }
        }

        // Check if this is a new group (before adding)
        let is_new_group = !self.segment_groups.contains_key(&key);

        if is_new_group {
            eprintln!(
                "RUST_NEW_GROUP: group_id={}, front_kmer={}, back_kmer={} (normalized)",
                self.segment_groups.len(),
                key.kmer_front,
                key.kmer_back
            );

            // DEBUG: Print first few bases of first segment
            if self.segment_groups.is_empty() {
                let bases_to_print = final_segment.len().min(50);
                let base_str: String = final_segment[..bases_to_print]
                    .iter()
                    .map(|&b| match b {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        3 => 'T',
                        _ => 'N',
                    })
                    .collect();
                eprintln!(
                    "RUST_FIRST_SEG: sample={}, contig={}, len={}, first_bases={}",
                    sample_name,
                    contig_name,
                    final_segment.len(),
                    base_str
                );
            }
        }

        let seg_info = SegmentInfo {
            sample_name: sample_name.to_string(),
            contig_name: contig_name.to_string(),
            seg_part_no,
            data: final_segment,
            is_rev_comp: final_is_rev_comp,
        };

        // Add to group
        self.segment_groups
            .entry(key.clone())
            .or_default()
            .push(seg_info);
        self.total_segments += 1;

        // Record k-mer pairing for middle splitter optimization (if new group)
        // CRITICAL: Must happen here, not in flush_group, so terminators are available
        // for subsequent segments that need splitting
        // Use NORMALIZED k-mers (matching C++ lines 1018-1025 which use pk.first/pk.second)
        if is_new_group && key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
            // Add back to front's terminator vector and sort
            let vec = self.group_terminators.entry(key.kmer_front).or_default();
            vec.push(key.kmer_back);
            vec.sort_unstable();

            // Add front to back's terminator vector and sort (if different)
            if key.kmer_front != key.kmer_back {
                let vec = self.group_terminators.entry(key.kmer_back).or_default();
                vec.push(key.kmer_front);
                vec.sort_unstable();
            }

            if self.config.verbosity > 2 {
                println!(
                    "    Recorded terminators: {} <-> {}",
                    key.kmer_front, key.kmer_back
                );
            }
        }

        // MEMORY OPTIMIZATION: Flush immediately to write segments to disk
        // This prevents buffering thousands of segments in memory
        if self.config.group_flush_threshold > 0 {
            self.flush_group(&key)?;
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
                eprintln!(
                    "DEBUG: Creating group {} for key front_kmer={}, back_kmer={}",
                    group_id, key.kmer_front, key.kmer_back
                );
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
                // CRITICAL FIX: Use compress_reference_segment which checks repetitiveness
                // and uses tuple packing for low-repetitiveness data (matching C++ AGC)
                let (compressed, marker) = compress_reference_segment(&ref_segment.data)?;

                if compressed.len() + 1 < ref_segment.data.len() {
                    let mut compressed_with_marker = compressed;
                    compressed_with_marker.push(marker); // Marker byte: 0 = plain ZSTD, 1 = tuple-packed
                    self.archive.add_part(
                        ref_stream_id,
                        &compressed_with_marker,
                        ref_segment.data.len() as u64,
                    )?;
                } else {
                    // Uncompressed
                    self.archive.add_part(ref_stream_id, &ref_segment.data, 0)?;
                }

                // Register reference in collection (in_group_id = 0)
                // NOTE: Contig MUST already be registered during file read (line 1733)!
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
                    eprintln!(
                        "DEBUG: Group {} wrote reference segment to ref stream, len={}",
                        group_id,
                        ref_segment.data.len()
                    );
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

        // MEMORY OPTIMIZATION: Write all segments immediately instead of waiting for packs of 50
        // This matches C++ AGC's immediate write behavior
        metadata.pending_segments.append(&mut new_delta_segments);

        if metadata.pending_segments.is_empty() {
            return Ok(());
        }

        // Write ALL pending segments (drain completely for memory efficiency)
        let segments_to_pack: Vec<SegmentInfo> = metadata.pending_segments.drain(..).collect();
        let segments_to_write = segments_to_pack.len();

        if self.config.verbosity > 2 {
            eprintln!(
                "Group {group_id} writing {segments_to_write} segments immediately (memory-efficient mode)"
            );
        }

        // Step 1: Build all packs (LZ encoding + concatenation)
        let compression_level = self.config.compression_level;
        let mut packs: Vec<(Vec<SegmentInfo>, Vec<u8>)> = Vec::new();

        for pack_segments in segments_to_pack.chunks(PACK_CARDINALITY) {
            let mut packed_data = Vec::new();

            for seg_info in pack_segments.iter() {
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

                packed_data.extend_from_slice(&contig_data);
                packed_data.push(CONTIG_SEPARATOR);
            }

            packs.push((pack_segments.to_vec(), packed_data));
        }

        // Step 2: Compress all packs in parallel
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.num_threads)
            .build()
            .expect("Failed to build thread pool");
        let compressed_packs: Vec<(Vec<SegmentInfo>, Vec<u8>, Vec<u8>)> = pool.install(|| {
            packs
                .par_iter()
                .map(|(seg_infos, packed_data)| {
                    let compressed = compress_segment_configured(packed_data, compression_level)
                        .expect("Compression failed");
                    (seg_infos.clone(), packed_data.clone(), compressed)
                })
                .collect()
        });

        // Step 3: Write compressed packs and register segments sequentially
        for (pack_segments, packed_data, compressed) in compressed_packs {
            for (idx_in_pack, seg_info) in pack_segments.iter().enumerate() {
                // For LZ groups: in_group_id starts at 1 (reference is 0)
                // For raw groups: in_group_id starts at 0
                let global_in_group_id = if use_lz_encoding {
                    metadata.segments_written + idx_in_pack + 1
                } else {
                    metadata.segments_written + idx_in_pack
                };

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!("DEBUG: Group {} adding segment: contig={}, seg_part={}, in_group_id={}, data_len={}",
                        group_id, seg_info.contig_name, seg_info.seg_part_no, global_in_group_id, seg_info.data.len());
                }

                // Register in collection
                if self.config.verbosity > 2 && group_id == 801 {
                    eprintln!("REGISTER: group={}, in_group_id={}, pack_id would be {}, contig={}, seg_part={}",
                        group_id, global_in_group_id, (global_in_group_id - 1) / 50, seg_info.contig_name, seg_info.seg_part_no);
                }

                // NOTE: Contig MUST already be registered during file read (line 1733)!
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

                self.archive.add_part(
                    stream_id,
                    &compressed_with_marker,
                    packed_data.len() as u64,
                )?;
            } else {
                // Compression not beneficial, write uncompressed
                // CRITICAL: Do NOT add marker byte for uncompressed data!
                // C++ AGC expects raw data when uncompressed_size = 0

                if self.config.verbosity > 2 && (group_id == 16 || group_id == 17) {
                    eprintln!(
                        "DEBUG: Group {} pack: packed_len={}, writing uncompressed (no marker)",
                        group_id,
                        packed_data.len()
                    );
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
            println!(
                "  Periodic flush: flushing {} groups",
                self.segment_groups.len()
            );
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

        // Take all pending segments (even if less than PACK_CARDINALITY)
        let segments_to_pack: Vec<SegmentInfo> = metadata.pending_segments.drain(..).collect();

        // Pack and write segments (may be partial pack)
        // Track the base in_group_id for this flush operation
        let mut current_segments_written = metadata.segments_written;

        // Step 1: Build all packs (LZ encoding + concatenation)
        let compression_level = self.config.compression_level;
        let mut packs: Vec<(Vec<SegmentInfo>, Vec<u8>)> = Vec::new();

        for pack_segments in segments_to_pack.chunks(PACK_CARDINALITY) {
            let mut packed_data = Vec::new();

            for seg_info in pack_segments.iter() {
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
            }

            packs.push((pack_segments.to_vec(), packed_data));
        }

        // Step 2: Compress all packs in parallel
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.num_threads)
            .build()
            .expect("Failed to build thread pool");
        let compressed_packs: Vec<(Vec<SegmentInfo>, Vec<u8>, Vec<u8>)> = pool.install(|| {
            packs
                .par_iter()
                .map(|(seg_infos, packed_data)| {
                    let compressed = compress_segment_configured(packed_data, compression_level)
                        .expect("Compression failed");
                    (seg_infos.clone(), packed_data.clone(), compressed)
                })
                .collect()
        });

        // Step 3: Write compressed packs and register segments sequentially
        for (pack_segments, packed_data, compressed) in compressed_packs {
            for (idx_in_pack, seg_info) in pack_segments.iter().enumerate() {
                let global_in_group_id = if use_lz_encoding {
                    current_segments_written + idx_in_pack + 1
                } else {
                    current_segments_written + idx_in_pack
                };

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
            if compressed.len() + 1 < packed_data.len() {
                let mut compressed_with_marker = compressed;
                compressed_with_marker.push(0);
                self.archive.add_part(
                    stream_id,
                    &compressed_with_marker,
                    packed_data.len() as u64,
                )?;
            } else {
                self.archive.add_part(stream_id, &packed_data, 0)?;
            }

            current_segments_written += pack_segments.len();
        }

        // Update the metadata's segments_written counter
        metadata.segments_written = current_segments_written;

        Ok(())
    }

    /// Store params stream (C++ compatibility)
    /// Format: kmer_length (u32) + min_match_len (u32) + pack_cardinality (u32) + segment_size (u32)
    fn store_params_stream(&mut self) -> Result<()> {
        let mut params_data = Vec::new();
        let append_u32 = |data: &mut Vec<u8>, value: u32| {
            data.extend_from_slice(&value.to_le_bytes());
        };

        append_u32(&mut params_data, self.config.kmer_length);
        append_u32(&mut params_data, self.config.min_match_len);
        append_u32(&mut params_data, 50); // pack_cardinality
        append_u32(&mut params_data, self.config.segment_size);
        // Note: no_raw_groups is not included (C++ AGC doesn't store this)

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
            let mut group_sizes: Vec<usize> =
                self.segment_groups.values().map(|v| v.len()).collect();
            group_sizes.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending

            let total_groups = self.segment_groups.len() + self.total_groups_flushed;
            println!("\nGroup statistics:");
            println!("  Total unique groups: {total_groups}");
            println!(
                "  Groups with 1 segment: {}",
                group_sizes.iter().filter(|&&s| s == 1).count()
            );
            println!(
                "  Groups with 2+ segments: {}",
                group_sizes.iter().filter(|&&s| s > 1).count()
            );
            if !group_sizes.is_empty() {
                println!("  Largest group: {} segments", group_sizes[0]);
                println!(
                    "  Top 10 group sizes: {:?}",
                    &group_sizes[..group_sizes.len().min(10)]
                );
            }
        }

        // Flush any remaining groups (sorted for determinism)
        let mut remaining_keys: Vec<_> = self.segment_groups.keys().cloned().collect();
        remaining_keys.sort_unstable();
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
            println!("Actually flushed {groups_actually_flushed} groups");
        }

        // Flush any remaining buffered segments (partial packs - sorted for determinism)
        let mut all_group_keys: Vec<_> = self.group_metadata.keys().cloned().collect();
        all_group_keys.sort_unstable();
        if self.config.verbosity > 0 {
            println!(
                "Flushing buffered segments for {} groups",
                all_group_keys.len()
            );

            // Count how many have pending segments
            let pending_count: usize = all_group_keys
                .iter()
                .filter(|k| {
                    self.group_metadata
                        .get(k)
                        .map(|m| !m.pending_segments.is_empty())
                        .unwrap_or(false)
                })
                .count();
            println!("  Groups with pending segments: {pending_count}");
        }
        for key in all_group_keys {
            self.flush_group_final(&key)?;
        }

        // DEBUG: Print first 50 groups sorted by group_id to compare with C++
        if self.config.verbosity > 0 {
            println!("\n=== First 50 unique groups (sorted by group_id) ===");
            let mut groups_with_ids: Vec<(u32, u64, u64)> = self
                .group_metadata
                .iter()
                .map(|(key, meta)| (meta.group_id, key.kmer_front, key.kmer_back))
                .collect();
            groups_with_ids.sort_by_key(|(id, _, _)| *id);

            for (group_id, front_kmer, back_kmer) in groups_with_ids.iter().take(50) {
                println!("  Group {group_id}: front_kmer={front_kmer}, back_kmer={back_kmer}");
            }
        }

        // Store metadata streams
        self.store_params_stream()?;
        self.store_splitters_stream()?;
        self.store_segment_splitters_stream()?;

        // Store collection metadata
        if self.config.verbosity > 0 {
            println!("DEBUG: Storing collection metadata...");
            println!("  Number of samples: {}", self.collection.get_no_samples());
        }

        self.collection
            .store_batch_sample_names(&mut self.archive)?;

        let num_samples = self.collection.get_no_samples();
        if num_samples > 0 {
            if self.config.verbosity > 0 {
                println!("  Storing contig batch for {} samples", num_samples);
            }
            self.collection
                .store_contig_batch(&mut self.archive, 0, num_samples)?;
        } else {
            if self.config.verbosity > 0 {
                println!("  WARNING: No samples to store!");
            }
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
