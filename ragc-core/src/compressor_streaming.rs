// Streaming AGC Compressor
// Memory-efficient implementation that flushes groups incrementally

#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::needless_range_loop)]

/// Helper to estimate memory usage of a data structure
fn estimate_memory_mb<T>(items: &[T]) -> f64 {
    std::mem::size_of_val(items) as f64 / (1024.0 * 1024.0)
}

use crate::{
    contig_iterator::ContigIterator,
    genome_io::GenomeIO,
    kmer::{Kmer, KmerMode},
    lz_diff::LZDiff,
    segment::{split_at_splitters_with_size, MISSING_KMER},
    segment_compression::compress_segment_configured,
    splitters::determine_splitters,
};
use anyhow::{Context, Result};
use dashmap::DashMap;
use ragc_common::{
    stream_delta_name, stream_ref_name, Archive, CollectionV3, Contig, AGC_FILE_MAJOR,
    AGC_FILE_MINOR, CONTIG_SEPARATOR,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::Path;
use std::sync::{
    atomic::{AtomicU32, AtomicUsize, Ordering},
    Arc, Mutex, RwLock,
};
use std::thread;
use crossbeam::channel::bounded;

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
struct ContigTask {
    sample_name: String,
    contig_name: String,
    sequence: Contig,
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
                let compressed =
                    compress_segment_configured(&segment.data, config.compression_level)?;

                let (compressed_data, uncompressed_size) =
                    if compressed.len() + 1 < segment.data.len() {
                        let mut compressed_with_marker = compressed;
                        compressed_with_marker.push(0);
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

        for (idx, seg_info) in self.pending_segments.iter().enumerate() {
            let contig_data = if use_lz_encoding {
                let mut lz_diff = LZDiff::new(config.min_match_len);
                if let Some(ref reference) = self.reference {
                    lz_diff.prepare(reference);
                    lz_diff.encode(&seg_info.data)
                } else {
                    seg_info.data.clone()
                }
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
            next_group_id: 0,
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
    pub fn add_fasta_file(&mut self, sample_name: &str, fasta_path: &Path) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Processing sample: {sample_name} from {fasta_path:?}");
        }

        let mut reader =
            GenomeIO::<Box<dyn Read>>::open(fasta_path).context("Failed to open FASTA file")?;

        let mut contigs_processed = 0;

        // Read contigs with conversion (ASCII -> numeric)
        while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
            self.add_contig(sample_name, &contig_name, sequence)?;
            contigs_processed += 1;

            // Periodic flush of all groups
            if self.config.periodic_flush_interval > 0
                && self.contigs_since_flush >= self.config.periodic_flush_interval
            {
                self.flush_all_groups()?;
            }

            if self.config.verbosity > 1 && contigs_processed % 100 == 0 {
                println!(
                    "  Processed {contigs_processed} contigs, {} groups in memory, {} flushed",
                    self.segment_groups.len(),
                    self.total_groups_flushed
                );
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
                            println!("Using {reference_sample} as reference to find splitters...");
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
                println!("Collected {contig_count} reference contigs from {reference_sample}");
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

            while let Some((_full_header, sample_name, contig_name, sequence)) =
                reader.read_contig_with_sample()?
            {
                if !sequence.is_empty() {
                    contig_count += 1;

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
                        let seg_info = Self::prepare_segment_info(
                            &self.config,
                            &sample_name,
                            &contig_name,
                            seg_idx,
                            segment.data,
                            segment.front_kmer,
                            segment.back_kmer,
                        )?;

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
                    Some(10000 + gid as usize) // Ref streams at offset 10000
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
                                    let (final_data, uncompressed_size) = if compressed_data.len()
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
                let stream_name = if pack.stream_id >= 10000 {
                    // Reference stream
                    let group_id = (pack.stream_id - 10000) as u32;
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
    pub fn add_contigs_with_splitters(
        &mut self,
        mut iterator: Box<dyn ContigIterator>,
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("=== Processing contigs with unified iterator ===");
            println!();
        }

        // Pass 1: Collect reference contigs for splitter finding
        if self.config.verbosity > 0 {
            println!("=== Pass 1: Finding splitter k-mers from reference ===");
        }

        let mut reference_contigs = Vec::new();
        let mut reference_sample = String::new();

        let mut contig_count = 0;
        while let Some((sample_name, _contig_name, sequence)) = iterator.next_contig()? {
            if !sequence.is_empty() {
                if reference_sample.is_empty() {
                    reference_sample = sample_name.clone();
                    if self.config.verbosity > 0 {
                        println!("Using {reference_sample} as reference to find splitters...");
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
            println!("Collected {contig_count} reference contigs from {reference_sample}");
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
            println!();
            println!("=== Pass 2: Segmenting and compressing all samples ===");
        }

        drop(reference_contigs); // Free memory

        // Reset iterator for second pass
        iterator.reset()?;

        // Pass 2: Collect all contigs in memory (simple approach for now)
        // TODO: Optimize to streaming later if memory becomes an issue
        if self.config.verbosity > 0 {
            println!("=== Pass 2: Collecting all contigs ===");
        }

        let mut all_contigs = Vec::new();
        while let Some((sample_name, contig_name, sequence)) = iterator.next_contig()? {
            if !sequence.is_empty() {
                all_contigs.push((sample_name, contig_name, sequence));
            }
        }

        if self.config.verbosity > 0 {
            println!(
                "Collected {} contigs, starting parallel processing...",
                all_contigs.len()
            );
        }

        // Now use parallel pipeline
        use crossbeam::channel::{bounded, Receiver, Sender};

        type ContigTask = (String, String, Contig);

        // Reduced buffer sizes from 100 → 10 to lower memory pressure
        // Each buffered item contains Vec<u8> with segment data (~10KB+)
        // 100 slots × 3 channels = significant memory overhead
        let (contig_tx, contig_rx): (Sender<ContigTask>, Receiver<ContigTask>) = bounded(10);

        // THREE-STAGE PIPELINE for parallel compression
        let (uncompressed_tx, uncompressed_rx): (
            Sender<UncompressedPack>,
            Receiver<UncompressedPack>,
        ) = bounded(10);

        let (compressed_tx, compressed_rx): (Sender<CompressedPack>, Receiver<CompressedPack>) =
            bounded(10);

        let verbosity = self.config.verbosity;

        // Producer thread: feeds contigs from memory to channel
        let producer_handle = std::thread::spawn(move || -> Result<usize> {
            let contig_count = all_contigs.len();
            for (sample_name, contig_name, sequence) in all_contigs {
                contig_tx
                    .send((sample_name, contig_name, sequence))
                    .context("Failed to send contig to workers")?;
            }

            // Close channel to signal workers we're done
            drop(contig_tx);

            if verbosity > 0 {
                println!("Producer: finished sending {contig_count} contigs");
            }

            Ok(contig_count)
        });

        // Shared state for middle splitter optimization
        // Matches C++ AGC: single RwLock protecting entire map (seg_map_mtx)
        let group_terminators = Arc::new(RwLock::new(std::mem::take(&mut self.group_terminators)));

        // Per-group buffers for segment writing
        let group_writers =
            Arc::new(DashMap::<SegmentGroupKey, std::sync::Mutex<GroupWriter>>::new());

        // Group ID counter (atomic, shared across workers)
        let next_group_id = Arc::new(AtomicUsize::new(self.next_group_id as usize));

        // Collect contig registrations for sequential processing
        let contig_registrations = Arc::new(DashMap::new());

        let total_segments = Arc::new(AtomicUsize::new(self.total_segments));
        let total_bases = Arc::new(AtomicUsize::new(self.total_bases_processed));
        let config = Arc::new(self.config.clone());
        let splitters_arc = Arc::new(splitters);

        let contig_registrations_ref = contig_registrations.clone();

        // Worker threads: consume from channel and process
        let mut worker_handles = Vec::new();
        for _worker_id in 0..self.config.num_threads {
            let rx = contig_rx.clone();
            let uncompressed_tx = uncompressed_tx.clone();
            let compressed_tx = compressed_tx.clone();
            let writers = group_writers.clone();
            let _group_terms = group_terminators.clone();
            let registrations = contig_registrations.clone();
            let total_segs = total_segments.clone();
            let total_b = total_bases.clone();
            let cfg = config.clone();
            let splitters = splitters_arc.clone();
            let group_id_counter = next_group_id.clone();

            let worker_handle = std::thread::spawn(move || -> Result<usize> {
                let mut worker_contigs = 0;

                // Consume from channel until closed
                while let Ok((sample_name, contig_name, sequence)) = rx.recv() {
                    worker_contigs += 1;

                    // Record registration for later
                    registrations.insert((sample_name.clone(), contig_name.clone()), ());

                    // Segment at splitters
                    let segments = split_at_splitters_with_size(
                        &sequence,
                        &splitters,
                        cfg.kmer_length as usize,
                        cfg.segment_size as usize,
                    );

                    // Process each segment: buffer in per-group writers
                    for (seg_idx, segment) in segments.into_iter().enumerate() {
                        total_b.fetch_add(segment.data.len(), Ordering::Relaxed);
                        total_segs.fetch_add(1, Ordering::Relaxed);

                        // Prepare segment info
                        let seg_info = Self::prepare_segment_info(
                            &cfg,
                            &sample_name,
                            &contig_name,
                            seg_idx,
                            segment.data,
                            segment.front_kmer,
                            segment.back_kmer,
                        )?;

                        // Get or create group writer
                        let group_writer =
                            writers.entry(seg_info.key.clone()).or_insert_with(|| {
                                let gid = group_id_counter.fetch_add(1, Ordering::SeqCst) as u32;

                                // TEMPORARILY DISABLED: Update terminators to test if this is the bottleneck
                                // if seg_info.key.kmer_front != MISSING_KMER && seg_info.key.kmer_back != MISSING_KMER {
                                //     let mut map = group_terms.write().unwrap();
                                //     map.entry(seg_info.key.kmer_front)
                                //         .or_default()
                                //         .push(seg_info.key.kmer_back);
                                //     map.get_mut(&seg_info.key.kmer_front).unwrap().sort_unstable();
                                //
                                //     if seg_info.key.kmer_front != seg_info.key.kmer_back {
                                //         map.entry(seg_info.key.kmer_back)
                                //             .or_default()
                                //             .push(seg_info.key.kmer_front);
                                //         map.get_mut(&seg_info.key.kmer_back).unwrap().sort_unstable();
                                //     }
                                // }

                                let stream_id = gid;

                                const NO_RAW_GROUPS: u32 = 16;
                                let ref_stream_id = if gid >= NO_RAW_GROUPS {
                                    Some(10000 + gid as usize)
                                } else {
                                    None
                                };

                                std::sync::Mutex::new(GroupWriter::new(
                                    gid,
                                    stream_id as usize,
                                    ref_stream_id,
                                ))
                            });

                        // Add segment to group buffer
                        let mut writer_guard = group_writer.lock().unwrap();
                        if let Some(pack) = writer_guard.add_segment(seg_info.segment, &cfg)? {
                            drop(writer_guard); // Release lock before sending

                            match pack {
                                PackToWrite::Compressed(compressed_pack) => {
                                    // Reference pack - already compressed, send to writer
                                    compressed_tx.send(compressed_pack).context(
                                        "Failed to send compressed reference pack to writer",
                                    )?;
                                }
                                PackToWrite::Uncompressed(uncompressed_pack) => {
                                    // Buffer is full - send to compression pool
                                    uncompressed_tx.send(uncompressed_pack).context(
                                        "Failed to send uncompressed pack to compression pool",
                                    )?;
                                }
                            }
                        }
                    }
                }

                Ok(worker_contigs)
            });

            worker_handles.push(worker_handle);
        }

        // Compression thread pool: receives uncompressed packs and compresses in parallel
        let compression_level = config.compression_level;
        let num_compressors = self.config.num_threads;
        for _compressor_id in 0..num_compressors {
            let uncompressed_rx = uncompressed_rx.clone();
            let compressed_tx = compressed_tx.clone();

            std::thread::spawn(move || {
                while let Ok(uncompressed_pack) = uncompressed_rx.recv() {
                    // Compress the pack
                    let compressed = compress_segment_configured(
                        &uncompressed_pack.uncompressed_data,
                        compression_level,
                    );

                    if let Ok(compressed_data) = compressed {
                        // Decide whether to use compressed or raw data
                        let (final_data, uncompressed_size) = if compressed_data.len() + 1
                            < uncompressed_pack.uncompressed_data.len()
                        {
                            let mut data_with_marker = compressed_data;
                            data_with_marker.push(0); // Marker byte
                            (
                                data_with_marker,
                                uncompressed_pack.uncompressed_data.len() as u64,
                            )
                        } else {
                            (uncompressed_pack.uncompressed_data, 0)
                        };

                        // Send compressed pack to writer
                        let compressed_pack = CompressedPack {
                            group_id: uncompressed_pack.group_id,
                            stream_id: uncompressed_pack.stream_id,
                            compressed_data: final_data,
                            uncompressed_size,
                            segments: uncompressed_pack.segments,
                        };

                        let _ = compressed_tx.send(compressed_pack);
                    }
                }
            });
        }

        // Spawn a thread to signal when workers finish and flush remaining packs
        let writers_for_flush = group_writers.clone();
        let cfg_for_flush = config.clone();
        let (workers_done_tx, workers_done_rx) = bounded(1);
        std::thread::spawn(move || {
            let _ = producer_handle.join();
            for handle in worker_handles {
                let _ = handle.join();
            }

            // Workers done - flush all remaining buffers
            for entry in writers_for_flush.iter() {
                let mut writer = entry.value().lock().unwrap();
                if let Ok(Some(uncompressed_pack)) = writer.prepare_pack(&cfg_for_flush) {
                    let _ = uncompressed_tx.send(uncompressed_pack);
                }
            }

            drop(uncompressed_tx); // Close uncompressed channel
            let _ = workers_done_tx.send(());
        });

        // Drop original compressed_tx so compression pool has the only remaining senders
        drop(compressed_tx);

        // Register streams in archive
        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
        let mut stream_id_map = HashMap::new();

        // Write packs as they arrive from compression pool
        let mut packs_written = 0;
        let mut segments_written = 0;

        loop {
            match compressed_rx.try_recv() {
                Ok(pack) => {
                    // Register stream if first time seeing it
                    let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                        id
                    } else {
                        let stream_name = if pack.stream_id >= 10000 {
                            // Reference stream
                            let group_id = (pack.stream_id - 10000) as u32;
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
                        self.collection
                            .register_sample_contig(&seg_meta.sample_name, &seg_meta.contig_name)?;

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

                    // Write pack to archive
                    self.archive.add_part(
                        actual_stream_id,
                        &pack.compressed_data,
                        pack.uncompressed_size,
                    )?;
                    packs_written += 1;
                    segments_written += pack.segments.len();
                }
                Err(crossbeam::channel::TryRecvError::Empty) => {
                    // No pack available yet - check if workers are done
                    if workers_done_rx.try_recv().is_ok() {
                        // Workers done and channel empty - drain remaining
                        while let Ok(pack) = compressed_rx.recv() {
                            let actual_stream_id =
                                if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                                    id
                                } else {
                                    let stream_name = if pack.stream_id >= 10000 {
                                        let group_id = (pack.stream_id - 10000) as u32;
                                        stream_ref_name(archive_version, group_id)
                                    } else {
                                        stream_delta_name(archive_version, pack.group_id)
                                    };
                                    let actual_id = self.archive.register_stream(&stream_name);
                                    stream_id_map.insert(pack.stream_id, actual_id);
                                    actual_id
                                };

                            for seg_meta in &pack.segments {
                                self.collection.register_sample_contig(
                                    &seg_meta.sample_name,
                                    &seg_meta.contig_name,
                                )?;
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

                            self.archive.add_part(
                                actual_stream_id,
                                &pack.compressed_data,
                                pack.uncompressed_size,
                            )?;
                            packs_written += 1;
                            segments_written += pack.segments.len();
                        }
                        break;
                    }
                    // Workers not done yet - sleep briefly and retry
                    std::thread::sleep(std::time::Duration::from_millis(10));
                }
                Err(crossbeam::channel::TryRecvError::Disconnected) => break,
            }
        }

        if self.config.verbosity > 0 {
            println!(
                "Streaming complete: wrote {packs_written} packs ({segments_written} segments)"
            );
        }

        // Move shared state back from RwLock to HashMap
        {
            let rwlock = Arc::try_unwrap(group_terminators)
                .map_err(|_| anyhow::anyhow!("Failed to unwrap group_terminators Arc"))?;
            self.group_terminators = rwlock.into_inner().unwrap();
        }
        self.total_segments = Arc::try_unwrap(total_segments)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap total_segments Arc"))?
            .into_inner();
        self.total_bases_processed = Arc::try_unwrap(total_bases)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap total_bases Arc"))?
            .into_inner();

        if self.config.verbosity > 0 {
            println!("Processed {} contigs total", contig_registrations_ref.len());
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
    pub fn add_fasta_files_with_splitters(
        &mut self,
        fasta_paths: &[(String, &Path)],
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("=== Pass 1: Finding splitter k-mers across all genomes ===");
        }

        // Pass 1: Stream through FIRST genome to find splitters
        // This is the reference-based approach: use splitters from first genome
        // to segment all genomes at the same positions
        //
        // MEMORY-EFFICIENT: Streams through file twice instead of loading all contigs
        // For yeast (12MB genome): ~100MB Vec vs ~2.8GB loading all contigs!
        let (splitters, singletons, duplicates) = if let Some((ref_sample_name, ref_fasta_path)) =
            fasta_paths.first()
        {
            if self.config.verbosity > 0 {
                println!("Using {ref_sample_name} as reference to find splitters (streaming)...");
            }

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

        // Pass 2: Choose processing strategy based on adaptive mode
        if self.config.adaptive_mode {
            // Adaptive mode: must process samples sequentially (each sample may add splitters)
            self.add_fasta_files_sequential_adaptive(fasta_paths, &splitters)?;
        } else {
            // Non-adaptive mode: use C++ AGC exact architecture
            self.add_fasta_files_cpp_agc_style(fasta_paths, &splitters)?;
        }

        Ok(())
    }

    /// Sequential adaptive mode processing (old approach, for adaptive mode only)
    fn add_fasta_files_sequential_adaptive(
        &mut self,
        fasta_paths: &[(String, &Path)],
        splitters: &HashSet<u64>,
    ) -> Result<()> {
        // Pass 2: Re-read files and segment using splitters
        for (sample_name, fasta_path) in fasta_paths {
            if self.config.verbosity > 0 {
                println!("Processing sample: {sample_name} from {fasta_path:?}");
            }

            let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                .context(format!("Failed to re-open {sample_name}"))?;

            let mut contigs_processed = 0;

            // Adaptive mode: track contigs that need new splitters for this sample
            let mut contigs_needing_splitters: Vec<Contig> = Vec::new();

            // Clone the active splitters to avoid holding an immutable borrow of self
            // (needed so we can mutably borrow self later in the loop)
            let active_splitters = if self.config.adaptive_mode {
                self.current_splitters.clone()
            } else {
                splitters.clone()
            };

            while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
                if sequence.is_empty() {
                    continue;
                }

                // Register in collection
                self.collection
                    .register_sample_contig(sample_name, &contig_name)?;

                // Split at splitter positions with minimum segment size
                let segments = split_at_splitters_with_size(
                    &sequence,
                    &active_splitters,
                    self.config.kmer_length as usize,
                    self.config.segment_size as usize,
                );

                // Adaptive mode: check if this contig needs new splitters
                // Matching C++ AGC lines 2033-2039: find_new_splitters called for contigs >= segment_size
                if self.config.adaptive_mode && sequence.len() >= self.config.segment_size as usize
                {
                    // Check if segmentation is poor (max segment > 10x segment_size)
                    let max_segment_len = segments.iter().map(|s| s.data.len()).max().unwrap_or(0);
                    if max_segment_len > (self.config.segment_size as usize) * 10 {
                        // This contig needs new splitters
                        contigs_needing_splitters.push(sequence.clone());

                        if self.config.verbosity > 1 {
                            println!("  Contig {} has poor segmentation (max segment {} > {}), will find new splitters",
                                contig_name, max_segment_len, self.config.segment_size * 10);
                        }
                    }
                }

                if self.config.verbosity > 1 {
                    println!(
                        "  Contig {} split into {} segments",
                        contig_name,
                        segments.len()
                    );
                }

                // Add each segment
                for (seg_idx, segment) in segments.iter().enumerate() {
                    self.total_bases_processed += segment.data.len();

                    // CRITICAL: Check for segments < k (will cause C++ AGC errors)
                    if seg_idx > 0 && segment.data.len() < self.config.kmer_length as usize {
                        eprintln!(
                            "WARNING: Segment {} of contig {} has size {} < k={} bytes!",
                            seg_idx,
                            contig_name,
                            segment.data.len(),
                            self.config.kmer_length
                        );
                        eprintln!(
                            "  Sample: {}, front_kmer={}, back_kmer={}",
                            sample_name, segment.front_kmer, segment.back_kmer
                        );
                        eprintln!("  This will cause 'Corrupted archive!' errors in C++ AGC");
                    }

                    if self.config.verbosity > 2 {
                        println!(
                            "    Contig {}, Segment {}: front_kmer={}, back_kmer={}, len={}",
                            contig_name,
                            seg_idx,
                            segment.front_kmer,
                            segment.back_kmer,
                            segment.data.len()
                        );
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
                if self.config.periodic_flush_interval > 0
                    && self.contigs_since_flush >= self.config.periodic_flush_interval
                {
                    self.flush_all_groups()?;
                }

                if self.config.verbosity > 1 && contigs_processed % 100 == 0 {
                    println!(
                        "  Processed {} contigs, {} groups in memory, {} flushed",
                        contigs_processed,
                        self.segment_groups.len(),
                        self.total_groups_flushed
                    );
                }
            }

            if self.config.verbosity > 0 {
                println!("  Processed {contigs_processed} contigs from {sample_name}");
            }

            // Adaptive mode: after sample completes, find and add new splitters
            // Matching C++ AGC synchronization (lines 2186-2191 and 1154-1177)
            if self.config.adaptive_mode && !contigs_needing_splitters.is_empty() {
                if self.config.verbosity > 0 {
                    println!(
                        "Adaptive mode: Finding new splitters from {} problematic contigs...",
                        contigs_needing_splitters.len()
                    );
                }

                let mut new_splitters_for_sample = HashSet::new();

                for contig in &contigs_needing_splitters {
                    let contig_new_splitters = self.find_new_splitters(contig);
                    new_splitters_for_sample.extend(contig_new_splitters);
                }

                if !new_splitters_for_sample.is_empty() {
                    if self.config.verbosity > 0 {
                        println!("Adaptive mode: Adding {} new splitters to global set (was {} splitters)",
                            new_splitters_for_sample.len(), self.current_splitters.len());
                    }

                    // Add new splitters to current set for use with next sample
                    self.current_splitters.extend(new_splitters_for_sample);

                    if self.config.verbosity > 0 {
                        println!(
                            "Adaptive mode: Now have {} total splitters",
                            self.current_splitters.len()
                        );
                    }
                }
            }
        }

        Ok(())
    }

    /// Parallel non-adaptive mode processing (3-phase architecture matching C++ AGC!)
    fn add_fasta_files_parallel_non_adaptive(
        &mut self,
        fasta_paths: &[(String, &Path)],
        splitters: &HashSet<u64>,
    ) -> Result<()> {
        // ==================================================================
        // PHASE 1: Single-threaded segmentation (collect ALL segments from ALL samples)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 1: Segmenting all genomes...");
        }

        let mut all_segments: Vec<PreparedSegment> = Vec::new();
        let mut total_segments_count = 0;
        let mut total_bases_count = 0;

        for (sample_name, fasta_path) in fasta_paths {
            if self.config.verbosity > 0 {
                println!("  Reading sample: {sample_name} from {fasta_path:?}");
            }

            let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                .context(format!("Failed to open {sample_name}"))?;

            while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
                if sequence.is_empty() {
                    continue;
                }

                // Register in collection
                self.collection
                    .register_sample_contig(sample_name, &contig_name)?;

                // Split at splitter positions with minimum segment size
                let segments = split_at_splitters_with_size(
                    &sequence,
                    splitters,
                    self.config.kmer_length as usize,
                    self.config.segment_size as usize,
                );

                // Collect all segments
                for (seg_idx, segment) in segments.into_iter().enumerate() {
                    total_bases_count += segment.data.len();
                    total_segments_count += 1;

                    // CRITICAL: Check for segments < k (will cause C++ AGC errors)
                    if seg_idx > 0 && segment.data.len() < self.config.kmer_length as usize {
                        eprintln!(
                            "WARNING: Segment {} of contig {} has size {} < k={} bytes!",
                            seg_idx,
                            contig_name,
                            segment.data.len(),
                            self.config.kmer_length
                        );
                        eprintln!("  This will cause 'Corrupted archive!' errors in C++ AGC");
                    }

                    // Prepare segment info for grouping
                    let seg_info = Self::prepare_segment_info(
                        &self.config,
                        sample_name,
                        &contig_name,
                        seg_idx,
                        segment.data,
                        segment.front_kmer,
                        segment.back_kmer,
                    )?;
                    all_segments.push(seg_info);
                }
            }
        }

        if self.config.verbosity > 0 {
            println!(
                "Phase 1 complete: collected {total_segments_count} segments ({total_bases_count} bases)"
            );

            // Memory profiling
            let segment_mem = estimate_memory_mb(&all_segments);
            eprintln!(
                "[MEMORY] all_segments: {:.2} MB ({} items × {} bytes)",
                segment_mem,
                all_segments.len(),
                std::mem::size_of::<PreparedSegment>()
            );

            // Estimate actual data size
            let total_data_mb: f64 = all_segments
                .iter()
                .map(|s| s.segment.data.len())
                .sum::<usize>() as f64
                / (1024.0 * 1024.0);
            eprintln!("[MEMORY] Segment data payload: {total_data_mb:.2} MB");
        }

        // ==================================================================
        // PHASE 2: Single-threaded grouping (group segments by k-mer keys)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 2: Grouping segments by k-mer keys...");
        }

        let mut groups: HashMap<SegmentGroupKey, Vec<PreparedSegment>> = HashMap::new();
        for segment in all_segments {
            groups.entry(segment.key.clone()).or_default().push(segment);
        }

        let num_groups = groups.len();
        let groups_vec: Vec<(SegmentGroupKey, Vec<PreparedSegment>)> = groups.into_iter().collect();

        if self.config.verbosity > 0 {
            println!("Phase 2 complete: created {num_groups} groups");

            // Memory profiling
            let groups_mem = estimate_memory_mb(&groups_vec);
            eprintln!(
                "[MEMORY] groups_vec: {:.2} MB ({} groups × {} bytes)",
                groups_mem,
                groups_vec.len(),
                std::mem::size_of::<(SegmentGroupKey, Vec<PreparedSegment>)>()
            );
        }

        // ==================================================================
        // PHASE 3: Parallel group processing (matching C++ AGC!)
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
                // Assign group ID sequentially
                let gid = next_group_id.fetch_add(1, Ordering::SeqCst) as u32;
                let stream_id = gid as usize;

                const NO_RAW_GROUPS: u32 = 16;
                let ref_stream_id = if gid >= NO_RAW_GROUPS {
                    Some(10000 + gid as usize) // Ref streams at offset 10000
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
                                    let (final_data, uncompressed_size) = if compressed_data.len()
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

            // Memory profiling
            let packs_mem = estimate_memory_mb(&all_packs);
            let total_compressed_data: usize = all_packs.iter().map(|p| p.compressed_data.len()).sum();
            eprintln!(
                "[MEMORY] all_packs Vec: {:.2} MB ({} packs × {} bytes struct)",
                packs_mem,
                all_packs.len(),
                std::mem::size_of::<CompressedPack>()
            );
            eprintln!("[MEMORY] Compressed data payload: {:.2} MB", total_compressed_data as f64 / (1024.0 * 1024.0));
        }

        // ==================================================================
        // PHASE 4: Sequential write (like C++ AGC's single-threaded registration)
        // ==================================================================
        if self.config.verbosity > 0 {
            println!("Phase 4: Writing packs to archive sequentially...");
        }

        // Map logical stream IDs (from Phase 3) → actual archive stream IDs
        let mut stream_id_map: HashMap<usize, usize> = HashMap::new();
        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;

        // First pass: Register all streams and build ID mapping
        for pack in &all_packs {
            if let std::collections::hash_map::Entry::Vacant(e) =
                stream_id_map.entry(pack.stream_id)
            {
                let stream_name = if pack.stream_id >= 10000 {
                    // Reference stream
                    let group_id = (pack.stream_id - 10000) as u32;
                    stream_ref_name(archive_version, group_id)
                } else {
                    // Delta stream
                    stream_delta_name(archive_version, pack.group_id)
                };
                let actual_id = self.archive.register_stream(&stream_name);
                e.insert(actual_id);

                // Also create group metadata entry for delta streams
                if pack.stream_id < 10000 {
                    let group_id = pack.group_id;
                    let key = SegmentGroupKey {
                        kmer_front: 0, // We don't have access to the key here, but it's not needed for finalization
                        kmer_back: 0,
                    };

                    const NO_RAW_GROUPS: u32 = 16;
                    let ref_stream_id = if group_id >= NO_RAW_GROUPS {
                        let _ref_name = stream_ref_name(archive_version, group_id);
                        // Check if ref stream is already registered
                        let ref_logical = 10000 + group_id as usize;
                        stream_id_map.get(&ref_logical).copied()
                    } else {
                        None
                    };

                    self.group_metadata
                        .entry(key)
                        .or_insert_with(|| GroupMetadata {
                            group_id,
                            stream_id: actual_id,
                            ref_stream_id,
                            reference: None,
                            ref_written: false,
                            segments_written: 0,
                            pending_segments: Vec::new(),
                            is_flushed: false,
                        });
                }
            }
        }

        // Second pass: Write packs and register segments
        let mut packs_written = 0;
        for pack in all_packs {
            let actual_stream_id = *stream_id_map
                .get(&pack.stream_id)
                .expect("Stream ID should be registered");

            // Register segments in collection (use logical group_id, not actual_stream_id!)
            for seg_meta in &pack.segments {
                if self.config.verbosity > 2 {
                    eprintln!("[DEBUG] Registering segment: sample={}, contig={}, part={}, group_id={}, in_group={}",
                        seg_meta.sample_name, seg_meta.contig_name, seg_meta.seg_part_no,
                        pack.group_id, seg_meta.in_group_id);
                }
                self.collection
                    .register_sample_contig(&seg_meta.sample_name, &seg_meta.contig_name)?;
                self.collection.add_segment_placed(
                    &seg_meta.sample_name,
                    &seg_meta.contig_name,
                    seg_meta.seg_part_no,
                    pack.group_id, // Use logical group_id, collection will compute stream name from this
                    seg_meta.in_group_id,
                    seg_meta.is_rev_comp,
                    seg_meta.data_len,
                )?;
            }

            // Write pack data to archive
            self.archive.add_part(
                actual_stream_id,
                &pack.compressed_data,
                pack.uncompressed_size,
            )?;

            packs_written += 1;
            if self.config.verbosity > 1 && packs_written % 100 == 0 {
                println!("  Written {packs_written} packs...");
            }
        }

        if self.config.verbosity > 0 {
            println!("Phase 4 complete: wrote {packs_written} packs");
        }

        // Update state
        self.total_segments = total_segments_count;
        self.total_bases_processed = total_bases_count;

        // Update next group ID for future calls
        self.next_group_id = Arc::try_unwrap(next_group_id)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc"))?
            .into_inner() as u32;

        if self.config.verbosity > 0 {
            println!(
                "Compression complete!\nTotal bases: {total_bases_count}\nTotal groups: {num_groups}"
            );
        }

        Ok(())
    }

    /// Streaming architecture: process segments as they're read, write immediately
    /// This replaces the batch loading approach to reduce memory from 731 MB → ~235 MB
    fn add_fasta_files_streaming(
        &mut self,
        fasta_paths: &[(String, &Path)],
        splitters: &HashSet<u64>,
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Starting streaming compression pipeline...");
        }

        // Create single shared segment channel - workers pull from it
        let (segment_tx, segment_rx) = bounded::<PreparedSegment>(100);

        // ================================================================
        // SEGMENTATION THREAD: Read FASTAs and stream segments
        // ================================================================
        let segmentation_handle = {
            // Convert &Path to PathBuf for 'static lifetime
            let fasta_paths: Vec<(String, std::path::PathBuf)> = fasta_paths
                .iter()
                .map(|(name, path)| (name.clone(), path.to_path_buf()))
                .collect();
            let splitters = splitters.clone();
            let config = self.config.clone();

            std::thread::spawn(move || -> Result<(usize, usize)> {
                let mut total_segments = 0;
                let mut total_bases = 0;

                for (sample_name, fasta_path) in &fasta_paths {
                    if config.verbosity > 0 {
                        println!("  Reading sample: {sample_name} from {fasta_path:?}");
                    }

                    let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                        .context(format!("Failed to open {sample_name}"))?;

                    while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
                        if sequence.is_empty() {
                            continue;
                        }

                        // Split at splitter positions with minimum segment size
                        let segments = split_at_splitters_with_size(
                            &sequence,
                            &splitters,
                            config.kmer_length as usize,
                            config.segment_size as usize,
                        );

                        for (seg_idx, segment) in segments.into_iter().enumerate() {
                            total_bases += segment.data.len();
                            total_segments += 1;

                            // CRITICAL: Check for segments < k (will cause C++ AGC errors)
                            if seg_idx > 0 && segment.data.len() < config.kmer_length as usize {
                                eprintln!(
                                    "WARNING: Segment {} of contig {} has size {} < k={} bytes!",
                                    seg_idx,
                                    contig_name,
                                    segment.data.len(),
                                    config.kmer_length
                                );
                                eprintln!("  This will cause 'Corrupted archive!' errors in C++ AGC");
                            }

                            // Prepare segment info
                            let seg_info = Self::prepare_segment_info(
                                &config,
                                &sample_name,
                                &contig_name,
                                seg_idx,
                                segment.data,
                                segment.front_kmer,
                                segment.back_kmer,
                            )?;

                            // Send to shared channel - workers pull as available
                            segment_tx
                                .send(seg_info)
                                .context("Failed to send segment")?;
                        }
                    }
                }

                // Signal completion to workers
                drop(segment_tx);
                if config.verbosity > 0 {
                    println!("Segmentation complete: {total_segments} segments ({total_bases} bases)");
                }
                Ok((total_segments, total_bases))
            })
        };

        // ================================================================
        // WRITE CHANNEL: Workers send packs here for immediate writing
        // Note: Writer runs in main thread because Archive/Collection aren't Send
        // ================================================================
        let (write_tx, write_rx) = bounded::<CompressedPack>(10);

        // ================================================================
        // WORKER THREADS: Process segments by group, send packs to writer
        // ================================================================
        let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));
        let config = self.config.clone();

        let worker_handles: Vec<_> = (0..self.config.num_threads)
            .map(|_worker_id| {
                let segment_rx = segment_rx.clone();
                let write_tx = write_tx.clone();
                let next_group_id = next_group_id.clone();
                let config = config.clone();

                std::thread::spawn(move || -> Result<()> {
                    // Each worker maintains its own groups (by key)
                    let mut my_groups: HashMap<SegmentGroupKey, (u32, GroupWriter)> = HashMap::new();

                    let mut segments_processed = 0;
                    while let Ok(prepared_seg) = segment_rx.recv() {
                        segments_processed += 1;

                        // Process segment - work stealing from shared channel

                        // Get or create group
                        let (_group_id, group_writer) = my_groups
                            .entry(prepared_seg.key.clone())
                            .or_insert_with(|| {
                                let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                                let stream_id = gid as usize;
                                let ref_stream_id = if gid >= 16 {
                                    Some(10000 + gid as usize)
                                } else {
                                    None
                                };
                                (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
                            });

                        // Add segment, may produce pack
                        if let Some(pack) = group_writer.add_segment(prepared_seg.segment.clone(), &config)? {
                            match pack {
                                PackToWrite::Compressed(compressed_pack) => {
                                    write_tx.send(compressed_pack).context("Failed to send pack to writer")?;
                                }
                                PackToWrite::Uncompressed(uncompressed_pack) => {
                                    // Compress immediately
                                    let compressed_data = compress_segment_configured(
                                        &uncompressed_pack.uncompressed_data,
                                        config.compression_level,
                                    )?;

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

                                    let compressed_pack = CompressedPack {
                                        group_id: uncompressed_pack.group_id,
                                        stream_id: uncompressed_pack.stream_id,
                                        compressed_data: final_data,
                                        uncompressed_size,
                                        segments: uncompressed_pack.segments,
                                    };
                                    write_tx.send(compressed_pack).context("Failed to send pack to writer")?;
                                }
                            }
                        }
                    }

                    // Flush all groups
                    for (_key, (_gid, mut group_writer)) in my_groups {
                        if let Some(uncompressed_pack) = group_writer.prepare_pack(&config)? {
                            // Compress final pack
                            let compressed_data = compress_segment_configured(
                                &uncompressed_pack.uncompressed_data,
                                config.compression_level,
                            )?;

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

                            let compressed_pack = CompressedPack {
                                group_id: uncompressed_pack.group_id,
                                stream_id: uncompressed_pack.stream_id,
                                compressed_data: final_data,
                                uncompressed_size,
                                segments: uncompressed_pack.segments,
                            };
                            write_tx.send(compressed_pack).context("Failed to send pack to writer")?;
                        }
                    }

                    Ok(())
                })
            })
            .collect();

        // Drop local references so channels close properly
        drop(segment_rx); // Workers have their clones
        drop(write_tx); // Workers have their clones

        // ================================================================
        // MAIN THREAD WRITING: Receive and write packs as they arrive
        // ================================================================
        let mut stream_id_map: HashMap<usize, usize> = HashMap::new();
        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
        let mut packs_written = 0;

        // Write packs as they arrive from workers
        while let Ok(pack) = write_rx.recv() {
            // Register stream on-demand
            let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                id
            } else {
                let stream_name = if pack.stream_id >= 10000 {
                    // Reference stream
                    let group_id = (pack.stream_id - 10000) as u32;
                    stream_ref_name(archive_version, group_id)
                } else {
                    // Delta stream
                    stream_delta_name(archive_version, pack.group_id)
                };
                let id = self.archive.register_stream(&stream_name);
                stream_id_map.insert(pack.stream_id, id);
                id
            };

            // Register segments in collection
            for seg_meta in &pack.segments {
                self.collection
                    .register_sample_contig(&seg_meta.sample_name, &seg_meta.contig_name)?;
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

            // Write pack data to archive immediately
            self.archive.add_part(
                actual_stream_id,
                &pack.compressed_data,
                pack.uncompressed_size,
            )?;

            packs_written += 1;
            if self.config.verbosity > 1 && packs_written % 100 == 0 {
                println!("  Written {packs_written} packs...");
            }
        }

        if self.config.verbosity > 0 {
            println!("Writing complete: wrote {packs_written} packs");
        }

        // Wait for segmentation to complete
        let (total_segments, total_bases) = segmentation_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Segmentation thread panicked"))??;

        // Wait for all workers
        for (idx, handle) in worker_handles.into_iter().enumerate() {
            handle
                .join()
                .map_err(|_| anyhow::anyhow!("Worker thread {} panicked", idx))??;
        }

        // Update state
        self.next_group_id = next_group_id.load(Ordering::SeqCst);
        self.total_segments = total_segments;
        self.total_bases_processed = total_bases;

        if self.config.verbosity > 0 {
            println!(
                "Streaming compression complete!\nTotal bases: {total_bases}\nTotal segments: {total_segments}\nPacks written: {packs_written}"
            );
        }

        Ok(())
    }

    /// C++ AGC exact architecture: contig-level parallelism with shared groups
    /// Matches C++ AGC's CBoundedPQueue + worker threads precisely
    /// Memory target: ~245 MB (vs C++ AGC's 205 MB)
    fn add_fasta_files_cpp_agc_style(
        &mut self,
        fasta_paths: &[(String, &Path)],
        splitters: &HashSet<u64>,
    ) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Starting C++ AGC-style compression (contig-level parallelism)...");
        }

        // Bounded queue for contigs (like C++ CBoundedPQueue)
        // C++ uses memory-based capacity, we use contig count for simplicity
        let queue_capacity = self.config.num_threads * 4;  // 4 contigs per thread in flight
        let (contig_tx, contig_rx) = bounded::<ContigTask>(queue_capacity);

        // Shared state: groups, archive, collection (like C++ AGC's shared groups with mutex)
        // Using DashMap for concurrent access (per-shard locking like C++ AGC's per-segment mutexes)
        let groups = Arc::new(DashMap::<SegmentGroupKey, (u32, GroupWriter)>::new());
        let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));
        let total_segments_counter = Arc::new(AtomicUsize::new(0));

        // Move archive and collection to Arc<Mutex> so writer thread can access them
        let archive = Arc::new(Mutex::new(std::mem::replace(&mut self.archive, Archive::new_writer())));
        let collection = Arc::new(Mutex::new(std::mem::replace(&mut self.collection, CollectionV3::new())));

        // ================================================================
        // READER THREAD: Push contigs to bounded queue
        // ================================================================
        let reader_handle = {
            let fasta_paths: Vec<(String, std::path::PathBuf)> = fasta_paths
                .iter()
                .map(|(name, path)| (name.clone(), path.to_path_buf()))
                .collect();
            let config = self.config.clone();

            thread::spawn(move || -> Result<(usize, usize)> {
                let mut total_contigs = 0;
                let mut total_bases = 0;

                for (sample_name, fasta_path) in &fasta_paths {
                    if config.verbosity > 0 {
                        println!("  Reading sample: {sample_name} from {fasta_path:?}");
                    }

                    let mut reader = GenomeIO::<Box<dyn Read>>::open(fasta_path)
                        .context(format!("Failed to open {sample_name}"))?;

                    while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
                        if sequence.is_empty() {
                            continue;
                        }

                        total_bases += sequence.len();
                        total_contigs += 1;

                        // Push contig to queue (blocks if queue is full - memory control!)
                        contig_tx.send(ContigTask {
                            sample_name: sample_name.clone(),
                            contig_name,
                            sequence,
                        }).context("Failed to send contig to queue")?;
                    }
                }

                drop(contig_tx);  // Signal workers: no more contigs
                if config.verbosity > 0 {
                    println!("Reader complete: {total_contigs} contigs ({total_bases} bases)");
                }
                Ok((total_contigs, total_bases))
            })
        };

        // ================================================================
        // WRITER THREAD: Receives packs and writes immediately (no buffering!)
        // ================================================================
        let (pack_tx, pack_rx) = bounded::<CompressedPack>(10);  // Small buffer for immediate writes

        let writer_handle = {
            let archive = archive.clone();
            let collection = collection.clone();
            let config = self.config.clone();

            thread::spawn(move || -> Result<usize> {
                let mut stream_id_map: HashMap<usize, usize> = HashMap::new();
                let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
                let mut packs_written = 0;

                while let Ok(pack) = pack_rx.recv() {
                    // Lock archive and collection for writing
                    let mut archive_lock = archive.lock()
                        .map_err(|e| anyhow::anyhow!("Failed to lock archive: {}", e))?;
                    let mut collection_lock = collection.lock()
                        .map_err(|e| anyhow::anyhow!("Failed to lock collection: {}", e))?;

                    // Register stream on-demand
                    let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                        id
                    } else {
                        let stream_name = if pack.stream_id >= 10000 {
                            stream_ref_name(archive_version, (pack.stream_id - 10000) as u32)
                        } else {
                            stream_delta_name(archive_version, pack.group_id)
                        };
                        let id = archive_lock.register_stream(&stream_name);
                        stream_id_map.insert(pack.stream_id, id);
                        id
                    };

                    // Register segments in collection
                    for seg_meta in &pack.segments {
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
                    }

                    // Write pack to archive immediately
                    archive_lock.add_part(
                        actual_stream_id,
                        &pack.compressed_data,
                        pack.uncompressed_size,
                    )?;

                    packs_written += 1;
                    if config.verbosity > 1 && packs_written % 100 == 0 {
                        println!("  Written {} packs...", packs_written);
                    }

                    // Locks released here automatically
                }

                if config.verbosity > 0 {
                    println!("Writer thread complete: wrote {packs_written} packs immediately");
                }
                Ok(packs_written)
            })
        };

        // ================================================================
        // WORKER THREADS: Pull contigs, segment, add to shared groups, send packs to writer
        // ================================================================
        let worker_handles: Vec<_> = (0..self.config.num_threads)
            .map(|worker_id| {
                let contig_rx = contig_rx.clone();
                let pack_tx = pack_tx.clone();
                let groups = groups.clone();
                let next_group_id = next_group_id.clone();
                let total_segments_counter = total_segments_counter.clone();
                let splitters = splitters.clone();
                let config = self.config.clone();

                thread::spawn(move || -> Result<()> {

                    // Per-thread reused buffers (like C++ AGC's per-thread ZSTD_CCtx)
                    // Currently not using buffer pooling, but this is where it would go

                    while let Ok(task) = contig_rx.recv() {
                        // Split contig into segments (contig-level parallelism!)
                        let segments = split_at_splitters_with_size(
                            &task.sequence,
                            &splitters,
                            config.kmer_length as usize,
                            config.segment_size as usize,
                        );

                        for (seg_idx, segment) in segments.into_iter().enumerate() {
                            total_segments_counter.fetch_add(1, Ordering::SeqCst);

                            // Check for segments < k
                            if seg_idx > 0 && segment.data.len() < config.kmer_length as usize {
                                eprintln!(
                                    "WARNING: Segment {} of contig {} has size {} < k={}!",
                                    seg_idx, task.contig_name,
                                    segment.data.len(), config.kmer_length
                                );
                            }

                            // Prepare segment info
                            let seg_info = Self::prepare_segment_info(
                                &config,
                                &task.sample_name,
                                &task.contig_name,
                                seg_idx,
                                segment.data,
                                segment.front_kmer,
                                segment.back_kmer,
                            )?;

                            // Access shared groups with DashMap (per-shard locking, not global!)
                            let pack_opt = {
                                let mut group_entry = groups.entry(seg_info.key.clone())
                                    .or_insert_with(|| {
                                        let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                                        let stream_id = gid as usize;
                                        let ref_stream_id = if gid >= 16 {
                                            Some(10000 + gid as usize)
                                        } else {
                                            None
                                        };
                                        (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
                                    });

                                // Add segment, may produce pack (only this entry is locked!)
                                group_entry.1.add_segment(seg_info.segment, &config)?
                            };  // Entry lock released here

                            // Process pack if ready (outside mutex!) - send to writer immediately
                            if let Some(pack) = pack_opt {
                                let compressed_pack = match pack {
                                    PackToWrite::Compressed(cp) => cp,
                                    PackToWrite::Uncompressed(uncompressed_pack) => {
                                        // Compress immediately
                                        let compressed_data = compress_segment_configured(
                                            &uncompressed_pack.uncompressed_data,
                                            config.compression_level,
                                        )?;

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

                                        CompressedPack {
                                            group_id: uncompressed_pack.group_id,
                                            stream_id: uncompressed_pack.stream_id,
                                            compressed_data: final_data,
                                            uncompressed_size,
                                            segments: uncompressed_pack.segments,
                                        }
                                    }
                                };
                                // Send to writer immediately - no buffering!
                                pack_tx.send(compressed_pack)
                                    .context("Failed to send pack to writer")?;
                            }
                        }
                    }

                    if config.verbosity > 1 {
                        println!("Worker {worker_id}: completed");
                    }

                    Ok(())
                })
            })
            .collect();

        // Drop contig receiver so workers know when no more contigs coming
        drop(contig_rx);

        // ================================================================
        // WAIT FOR WORKERS TO COMPLETE
        // ================================================================
        for (idx, handle) in worker_handles.into_iter().enumerate() {
            handle
                .join()
                .map_err(|_| anyhow::anyhow!("Worker thread {} panicked", idx))??;
        }

        let total_segments_processed = total_segments_counter.load(Ordering::SeqCst);

        if self.config.verbosity > 0 {
            println!("All workers complete: {total_segments_processed} segments processed");
        }

        // ================================================================
        // WAIT FOR READER
        // ================================================================
        let (total_contigs, total_bases) = reader_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Reader thread panicked"))??;

        // ================================================================
        // FLUSH REMAINING GROUPS - Send to writer thread
        // ================================================================
        // Take ownership of the groups HashMap
        // Take ownership of DashMap (no need for into_inner() - DashMap isn't wrapped in Mutex)
        let groups_map = Arc::try_unwrap(groups)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap Arc for groups"))?;

        let mut flush_packs_count = 0;
        for (_key, (_gid, mut group_writer)) in groups_map.into_iter() {
            if let Some(uncompressed_pack) = group_writer.prepare_pack(&self.config)? {
                // Compress final pack
                let compressed_data = compress_segment_configured(
                    &uncompressed_pack.uncompressed_data,
                    self.config.compression_level,
                )?;

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

                let compressed_pack = CompressedPack {
                    group_id: uncompressed_pack.group_id,
                    stream_id: uncompressed_pack.stream_id,
                    compressed_data: final_data,
                    uncompressed_size,
                    segments: uncompressed_pack.segments,
                };

                // Send flush pack to writer thread
                pack_tx.send(compressed_pack)
                    .context("Failed to send flush pack to writer")?;
                flush_packs_count += 1;
            }
        }

        if self.config.verbosity > 0 {
            println!("Flush complete: sent {flush_packs_count} final packs to writer");
        }

        // Drop pack_tx so writer thread knows there are no more packs
        drop(pack_tx);

        // ================================================================
        // WAIT FOR WRITER THREAD TO COMPLETE
        // ================================================================
        let packs_written = writer_handle
            .join()
            .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

        if self.config.verbosity > 0 {
            println!("Writer thread complete: {packs_written} total packs written");
        }

        // ================================================================
        // UPDATE STATE
        // ================================================================
        // Restore archive and collection from Arc<Mutex>
        self.archive = Arc::try_unwrap(archive)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap archive Arc"))?
            .into_inner()
            .map_err(|e| anyhow::anyhow!("Failed to get archive mutex: {}", e))?;

        self.collection = Arc::try_unwrap(collection)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap collection Arc"))?
            .into_inner()
            .map_err(|e| anyhow::anyhow!("Failed to get collection mutex: {}", e))?;

        self.next_group_id = next_group_id.load(Ordering::SeqCst);
        self.total_segments = total_segments_processed;
        self.total_bases_processed = total_bases;

        if self.config.verbosity > 0 {
            println!(
                "C++ AGC-style compression complete!\nContigs: {total_contigs}\nSegments: {total_segments_processed}\nBases: {total_bases}\nPacks: {packs_written}"
            );
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
                    let seg1 = final_segment[..seg2_start_pos + kmer_len].to_vec();
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

    /// Prepare segment info for buffering (k-mer normalization, RC, terminator tracking)
    /// Returns the normalized segment ready to add to a group buffer
    #[allow(clippy::too_many_arguments)]
    fn prepare_segment_info(
        config: &StreamingCompressorConfig,
        sample_name: &str,
        contig_name: &str,
        seg_part_no: usize,
        segment: Contig,
        kmer_front: u64,
        kmer_back: u64,
    ) -> Result<PreparedSegment> {
        // K-mer orientation logic (matching add_segment_with_kmers)
        let (final_front, final_back) = if kmer_front != MISSING_KMER && kmer_back == MISSING_KMER {
            let k = config.kmer_length as usize;
            if segment.len() >= k {
                let mut front = Kmer::new(config.kmer_length, KmerMode::Canonical);
                for i in 0..k {
                    if segment[i] > 3 {
                        front.reset();
                        break;
                    }
                    front.insert(segment[i] as u64);
                }

                if front.is_full() && front.data() == kmer_front {
                    if front.is_dir_oriented() {
                        (kmer_front, MISSING_KMER)
                    } else {
                        (MISSING_KMER, kmer_front)
                    }
                } else {
                    (kmer_front, kmer_back)
                }
            } else {
                (kmer_front, kmer_back)
            }
        } else if kmer_front == MISSING_KMER && kmer_back != MISSING_KMER {
            let k = config.kmer_length as usize;
            if segment.len() >= k {
                let mut back = Kmer::new(config.kmer_length, KmerMode::Canonical);
                let start = segment.len() - k;
                for i in 0..k {
                    if segment[start + i] > 3 {
                        back.reset();
                        break;
                    }
                    back.insert(segment[start + i] as u64);
                }

                if back.is_full() && back.data() == kmer_back {
                    if !back.is_dir_oriented() {
                        (kmer_back, MISSING_KMER)
                    } else {
                        (MISSING_KMER, kmer_back)
                    }
                } else {
                    (kmer_front, kmer_back)
                }
            } else {
                (kmer_front, kmer_back)
            }
        } else {
            (kmer_front, kmer_back)
        };

        // Normalize key
        let (key, needs_rc) = SegmentGroupKey::new_normalized(final_front, final_back);
        let final_is_rev_comp = needs_rc; // Original is_rev_comp removed since segments come from splitters (not RC'd yet)

        // Reverse complement if needed
        let final_segment = if needs_rc {
            segment
                .iter()
                .rev()
                .map(|&b| match b {
                    0 => 3,
                    1 => 2,
                    2 => 1,
                    3 => 0,
                    _ => b,
                })
                .collect()
        } else {
            segment
        };

        // NOTE: Terminator updates moved to group creation in parallel worker
        // (matches C++ AGC: lock only when v_segments[group_id] == nullptr)

        // Return prepared segment
        Ok(PreparedSegment {
            key,
            segment: SegmentInfo {
                sample_name: sample_name.to_string(),
                contig_name: contig_name.to_string(),
                seg_part_no,
                data: final_segment,
                is_rev_comp: final_is_rev_comp,
            },
        })
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
                // Try compressing the reference with configured level
                let compressed =
                    compress_segment_configured(&ref_segment.data, self.config.compression_level)?;

                if compressed.len() + 1 < ref_segment.data.len() {
                    let mut compressed_with_marker = compressed;
                    compressed_with_marker.push(0); // Marker byte: 0 = standard ZSTD
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
