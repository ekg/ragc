// Streaming AGC Compressor
// Memory-efficient implementation that flushes groups incrementally

use crate::{
    genome_io::{GenomeIO, parse_sample_from_header},
    kmer::{Kmer, KmerMode},
    lz_diff::LZDiff,
    segment::{split_at_splitters_with_size, MISSING_KMER},
    segment_compression::{compress_segment, compress_segment_configured},
    splitters::determine_splitters,
};
use anyhow::{Context, Result};
use ragc_common::{
    stream_delta_name, stream_ref_name, Archive, CollectionV3, Contig, AGC_FILE_MAJOR, AGC_FILE_MINOR,
    CONTIG_SEPARATOR,
};
use dashmap::DashMap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::sync::{Arc, Mutex, atomic::{AtomicUsize, Ordering}};

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
            // IMPORTANT: To match C++ AGC exactly, we need to keep all segments in memory
            // until finalize(). This is required for the middle splitter optimization
            // which needs to compute coding costs against existing segment data.
            group_flush_threshold: 0, // 0 = disabled (was 100)
            periodic_flush_interval: 0, // 0 = disabled (was 500)
            num_threads,
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
/// IMPORTANT: Keys are normalized so kmer_front <= kmer_back (matching C++ minmax logic)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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
                (SegmentGroupKey { kmer_front, kmer_back }, false)
            } else {
                // Swap and set needs_rc=true (matching C++ lines 1319-1323)
                (SegmentGroupKey { kmer_front: kmer_back, kmer_back: kmer_front }, true)
            }
        } else {
            // One or both are MISSING: don't normalize (C++ lines 1297-1372)
            (SegmentGroupKey { kmer_front, kmer_back }, false)
        }
    }
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

    // K-mer terminator tracking (for missing middle splitter optimization)
    // Maps each k-mer to the list of k-mers it pairs with in existing groups
    // This enables splitting segments to reuse existing groups
    group_terminators: HashMap<u64, Vec<u64>>,

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

            // DEBUG: Output first 20 splitters sorted
            let mut splitter_vec: Vec<u64> = splitters.iter().copied().collect();
            splitter_vec.sort_unstable();
            eprintln!("RUST_SPLITTERS: First 20 splitters:");
            for (i, &spl) in splitter_vec.iter().take(20).enumerate() {
                eprintln!("  {}: {}", i, spl);
            }

            println!();
            println!("=== Pass 2: Segmenting and compressing all samples ===");
        }

        drop(reference_contigs); // Free memory

        // Pass 2: TRUE STREAMING with producer-consumer channels (no buffering!)
        // Producer thread: reads file → sends to channel
        // Worker threads: receive from channel → process immediately
        // This matches C++ AGC's architecture exactly

        if self.config.verbosity > 0 {
            println!("=== Pass 2: Streaming processing with {} worker threads ===", self.config.num_threads);
        }

        use crossbeam::channel::{bounded, Sender, Receiver};
        use std::thread;

        type ContigTask = (String, String, Contig); // (sample_name, contig_name, sequence)

        // Bounded channel matching C++ AGC's CBoundedQueue
        // Queue capacity: allow 100 contigs in flight (balances memory vs throughput)
        let (contig_tx, contig_rx): (Sender<ContigTask>, Receiver<ContigTask>) = bounded(100);

        // Clone path for producer thread
        let fasta_path_clone = fasta_path.to_path_buf();
        let verbosity = self.config.verbosity;

        // Producer thread: reads FASTA and feeds channel (like C++ AGC main thread)
        let producer_handle = thread::spawn(move || -> Result<usize> {
            let mut reader = GenomeIO::<Box<dyn Read>>::open(&fasta_path_clone)
                .context("Failed to open multi-sample FASTA for streaming")?;

            let mut contig_count = 0;
            while let Some((_full_header, sample_name, contig_name, sequence)) = reader.read_contig_with_sample()? {
                if !sequence.is_empty() {
                    // Send to workers (blocks if queue full - backpressure!)
                    contig_tx.send((sample_name, contig_name, sequence))
                        .context("Failed to send contig to workers")?;
                    contig_count += 1;

                    if verbosity > 1 && contig_count % 1000 == 0 {
                        eprintln!("Producer: sent {} contigs to workers", contig_count);
                    }
                }
            }

            // Close channel to signal workers we're done
            drop(contig_tx);

            if verbosity > 0 {
                println!("Producer: finished sending {} contigs", contig_count);
            }

            Ok(contig_count)
        });

        // Make shared state thread-safe for parallel processing
        // Using DashMap for fine-grained locking (per-shard), matching C++ AGC's per-group mutexes
        let segment_groups = {
            let dashmap = DashMap::new();
            for (key, value) in std::mem::take(&mut self.segment_groups) {
                dashmap.insert(key, value);
            }
            Arc::new(dashmap)
        };
        let group_terminators = {
            let dashmap = DashMap::new();
            for (key, value) in std::mem::take(&mut self.group_terminators) {
                dashmap.insert(key, value);
            }
            Arc::new(dashmap)
        };

        // Collect contig registrations for sequential processing after workers finish
        let contig_registrations = Arc::new(DashMap::new());

        let total_segments = Arc::new(AtomicUsize::new(self.total_segments));
        let total_bases = Arc::new(AtomicUsize::new(self.total_bases_processed));
        let config = Arc::new(self.config.clone());
        let splitters_arc = Arc::new(splitters);

        // Worker threads: consume from channel and process (like C++ AGC worker threads)
        let mut worker_handles = Vec::new();
        for worker_id in 0..self.config.num_threads {
            let rx = contig_rx.clone();
            let seg_groups = segment_groups.clone();
            let group_terms = group_terminators.clone();
            let registrations = contig_registrations.clone();
            let total_segs = total_segments.clone();
            let total_b = total_bases.clone();
            let cfg = config.clone();
            let splitters = splitters_arc.clone();

            let worker_handle = thread::spawn(move || -> Result<usize> {
                let mut worker_contigs = 0;

                // Consume from channel until closed
                while let Ok((sample_name, contig_name, sequence)) = rx.recv() {
                    worker_contigs += 1;

                    // Record registration for later (collection registration must be sequential)
                    registrations.insert((sample_name.clone(), contig_name.clone()), ());

                    // Segment at splitters (parallel-safe, no shared state)
                    let segments = split_at_splitters_with_size(
                        &sequence,
                        &splitters,
                        cfg.kmer_length as usize,
                        cfg.segment_size as usize
                    );

                    // Add each segment to groups (DashMap provides fine-grained locking)
                    for (seg_idx, segment) in segments.into_iter().enumerate() {
                        total_b.fetch_add(segment.data.len(), Ordering::Relaxed);

                        Self::add_segment_with_kmers_threadsafe(
                            &cfg,
                            &seg_groups,
                            &group_terms,
                            &total_segs,
                            &sample_name,
                            &contig_name,
                            seg_idx,
                            segment.data,
                            segment.front_kmer,
                            segment.back_kmer,
                            false, // is_rev_comp
                        )?;
                    }

                    if cfg.verbosity > 2 && worker_contigs % 100 == 0 {
                        eprintln!("Worker {}: processed {} contigs", worker_id, worker_contigs);
                    }
                }

                if cfg.verbosity > 1 {
                    println!("Worker {}: processed {} contigs", worker_id, worker_contigs);
                }

                Ok(worker_contigs)
            });

            worker_handles.push(worker_handle);
        }

        // Wait for producer to finish
        let total_contigs = producer_handle.join()
            .map_err(|e| anyhow::anyhow!("Producer thread panicked: {:?}", e))??;

        // Wait for all workers to finish
        let mut total_worker_contigs = 0;
        for (worker_id, handle) in worker_handles.into_iter().enumerate() {
            let worker_contigs = handle.join()
                .map_err(|e| anyhow::anyhow!("Worker {} panicked: {:?}", worker_id, e))??;
            total_worker_contigs += worker_contigs;
        }

        if self.config.verbosity > 0 {
            println!("Streaming complete: {} contigs read, {} processed by workers",
                total_contigs, total_worker_contigs);
        }

        // Register all contigs in collection (sequential, after workers finish)
        let registrations_map = Arc::try_unwrap(contig_registrations)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap contig_registrations Arc"))?;

        if self.config.verbosity > 0 {
            println!("Registering {} contigs in collection...", registrations_map.len());
        }

        for ((sample_name, contig_name), _) in registrations_map.into_iter() {
            self.collection.register_sample_contig(&sample_name, &contig_name)?;
        }

        // Move shared state back from DashMap to HashMap
        {
            let dashmap = Arc::try_unwrap(segment_groups)
                .map_err(|_| anyhow::anyhow!("Failed to unwrap segment_groups Arc"))?;
            self.segment_groups = dashmap.into_iter().collect();
        }
        {
            let dashmap = Arc::try_unwrap(group_terminators)
                .map_err(|_| anyhow::anyhow!("Failed to unwrap group_terminators Arc"))?;
            self.group_terminators = dashmap.into_iter().collect();
        }
        self.total_segments = Arc::try_unwrap(total_segments)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap total_segments Arc"))?
            .into_inner();
        self.total_bases_processed = Arc::try_unwrap(total_bases)
            .map_err(|_| anyhow::anyhow!("Failed to unwrap total_bases Arc"))?
            .into_inner();

        if self.config.verbosity > 0 {
            println!("Processed {} contigs total", total_contigs);
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
            let back_kmer_val = if back.is_full() { back.data() } else { MISSING_KMER };
            let back_is_dir = back.is_dir_oriented();

            // Match C++ logic for one-k-mer cases (lines 1326-1372)
            if front_kmer_val != MISSING_KMER && back_kmer_val == MISSING_KMER {
                // Only front k-mer present (C++ lines 1326-1346)
                // C++ checks is_dir_oriented() directly
                let result = if front_is_dir {
                    (front_kmer_val, MISSING_KMER)  // dir-oriented: (kmer, MISSING)
                } else {
                    (MISSING_KMER, front_kmer_val)  // RC-oriented: (MISSING, kmer)
                };

                // DEBUG: Log a few cases
                if self.total_segments < 20 {
                    eprintln!("RUST_ONE_KMER: FRONT kmer={}, is_dir={}, result=({}, {})",
                        front_kmer_val, front_is_dir, result.0, result.1);
                }

                result
            } else if front_kmer_val == MISSING_KMER && back_kmer_val != MISSING_KMER {
                // Only back k-mer present (C++ lines 1348-1372)
                // CRITICAL: C++ calls kmer.swap_dir_rc() which swaps dir and rc
                // This INVERTS the is_dir_oriented() check!
                // So we need to use !back_is_dir
                let result = if !back_is_dir {
                    (back_kmer_val, MISSING_KMER)   // After swap, dir-oriented: (kmer, MISSING)
                } else {
                    (MISSING_KMER, back_kmer_val)   // After swap, RC-oriented: (MISSING, kmer)
                };

                // DEBUG: Log a few cases
                if self.total_segments < 20 {
                    eprintln!("RUST_ONE_KMER: BACK kmer={}, is_dir={}, !is_dir={}, result=({}, {})",
                        back_kmer_val, back_is_dir, !back_is_dir, result.0, result.1);
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
    fn find_middle_splitter(&self, kmer_front: u64, kmer_back: u64, segment: &Contig) -> Option<(u64, usize)> {
        // Both k-mers must be present (not MISSING_KMER)
        if kmer_front == MISSING_KMER || kmer_back == MISSING_KMER {
            if self.total_segments < 20 {
                eprintln!("RUST_MS_DEBUG: MISSING_KMER check failed");
            }
            return None;
        }

        // Get terminators for both k-mers
        let front_terminators = self.group_terminators.get(&kmer_front);
        let back_terminators = self.group_terminators.get(&kmer_back);

        if self.total_segments < 20 {
            eprintln!("RUST_MS_DEBUG: front_term={:?}, back_term={:?}",
                front_terminators.map(|v| v.len()),
                back_terminators.map(|v| v.len()));
        }

        let front_terminators = front_terminators?;
        let back_terminators = back_terminators?;

        if self.config.verbosity > 2 {
            eprintln!("    Looking for middle splitter: front={} has {} terminators, back={} has {} terminators",
                kmer_front, front_terminators.len(), kmer_back, back_terminators.len());
        }

        // Find intersection of terminators (shared k-mers)
        // C++ uses std::set_intersection on sorted vectors
        let mut shared_kmers: Vec<u64> = front_terminators.iter()
            .filter(|&&k| back_terminators.contains(&k))
            .copied()
            .filter(|&k| k != MISSING_KMER)
            .collect();

        if self.config.verbosity > 2 && !shared_kmers.is_empty() {
            eprintln!("    Found {} shared k-mers between terminators!", shared_kmers.len());
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
            eprintln!("    Found reference groups: ({}, {}) with {} segs, ({}, {}) with {} segs",
                kmer_front, middle_kmer, group1.len(),
                middle_kmer, kmer_back, group2.len());
            eprintln!("    Reference segment sizes: {} and {} bytes", ref_seg1.len(), ref_seg2.len());
        }

        // Compute reverse complement of segment (matching C++ exactly)
        let segment_rc: Vec<u8> = segment.iter().rev().map(|&b| {
            match b {
                0 => 3, // A -> T
                1 => 2, // C -> G
                2 => 1, // G -> C
                3 => 0, // T -> A
                _ => b, // N stays N
            }
        }).collect();

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
            let costs = lz_diff2.get_coding_cost_vector(segment, false);
            costs
        } else {
            // Use RC segment with prefix_costs=true
            let costs = lz_diff2.get_coding_cost_vector(&segment_rc, true);
            costs
        };

        if v_costs2.is_empty() || v_costs1.len() != v_costs2.len() {
            if self.config.verbosity > 2 {
                eprintln!("    Cost vector issue: v_costs2 len={}, v_costs1 len={}", v_costs2.len(), v_costs1.len());
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
            eprintln!("    Best split position: {} (after boundary check) with cost {}", best_pos, best_sum);
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
                        (kmer_front, MISSING_KMER)  // dir-oriented: (kmer, MISSING)
                    } else {
                        (MISSING_KMER, kmer_front)  // RC-oriented: (MISSING, kmer)
                    }
                } else {
                    (kmer_front, kmer_back)  // Fallback if extraction failed
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
                        (kmer_back, MISSING_KMER)   // After swap, dir-oriented: (kmer, MISSING)
                    } else {
                        (MISSING_KMER, kmer_back)   // After swap, RC-oriented: (MISSING, kmer)
                    }
                } else {
                    (kmer_front, kmer_back)  // Fallback if extraction failed
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
            segment.iter().rev().map(|&b| match b {
                0 => 3, // A -> T
                1 => 2, // C -> G
                2 => 1, // G -> C
                3 => 0, // T -> A
                _ => b, // N stays N
            }).collect()
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
            if let Some((middle_kmer, split_pos)) = self.find_middle_splitter(key.kmer_front, key.kmer_back, &final_segment) {
                // DEBUG: Log when middle splitter succeeds
                if self.total_segments < 10 {
                    eprintln!("RUST_MS_SUCCESS: middle_kmer={}, split_pos={}", middle_kmer, split_pos);
                }

                // Matching C++ lines 1400-1444: handle zero-size cases
                let left_size = split_pos;
                let right_size = final_segment.len() - split_pos;

                if left_size == 0 {
                    // C++ lines 1403-1410: left side is empty, just update the key
                    // Don't split - just add with the new key (middle_kmer, kmer_back)
                    if self.config.verbosity > 1 {
                        eprintln!("  Middle splitter: left_size=0, updating key from ({}, {}) to ({}, {})",
                            kmer_front, kmer_back, middle_kmer, kmer_back);
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
                        eprintln!("  Middle splitter: right_size=0, updating key from ({}, {}) to ({}, {})",
                            kmer_front, kmer_back, kmer_front, middle_kmer);
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
                        eprintln!("  Middle splitter optimization: splitting segment at position {}", split_pos);
                        eprintln!("    Original group: ({}, {})", kmer_front, kmer_back);
                        eprintln!("    Split into: ({}, {}) and ({}, {})", kmer_front, middle_kmer, middle_kmer, kmer_back);
                    }

                    // C++ line 1425: Create overlap of kmer_length/2 between segments
                    let kmer_len = self.config.kmer_length as usize;
                    let seg2_start_pos = if split_pos > kmer_len / 2 {
                        split_pos - kmer_len / 2
                    } else {
                        0
                    };

                    // C++ lines 1426-1428:
                    // segment2 = [seg2_start_pos, end)
                    // segment1 = [0, seg2_start_pos + kmer_length)
                    let seg1 = final_segment[..seg2_start_pos + kmer_len].to_vec();
                    let seg2 = final_segment[seg2_start_pos..].to_vec();

                    if self.config.verbosity > 1 {
                        eprintln!("    Overlap: seg1_len={}, seg2_len={}, overlap_size={}",
                            seg1.len(), seg2.len(), seg1.len() + seg2.len() - final_segment.len());
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
            eprintln!("RUST_NEW_GROUP: group_id={}, front_kmer={}, back_kmer={} (normalized)",
                self.segment_groups.len(), key.kmer_front, key.kmer_back);

            // DEBUG: Print first few bases of first segment
            if self.segment_groups.is_empty() {
                let bases_to_print = final_segment.len().min(50);
                let base_str: String = final_segment[..bases_to_print].iter()
                    .map(|&b| match b {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        3 => 'T',
                        _ => 'N',
                    })
                    .collect();
                eprintln!("RUST_FIRST_SEG: sample={}, contig={}, len={}, first_bases={}",
                    sample_name, contig_name, final_segment.len(), base_str);
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
        self.segment_groups.entry(key.clone()).or_default().push(seg_info);
        self.total_segments += 1;

        // Record k-mer pairing for middle splitter optimization (if new group)
        // CRITICAL: Must happen here, not in flush_group, so terminators are available
        // for subsequent segments that need splitting
        // Use NORMALIZED k-mers (matching C++ lines 1018-1025 which use pk.first/pk.second)
        if is_new_group && key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
            // Add back to front's terminator list
            self.group_terminators.entry(key.kmer_front).or_default().push(key.kmer_back);

            // Add front to back's terminator list (if different)
            if key.kmer_front != key.kmer_back {
                self.group_terminators.entry(key.kmer_back).or_default().push(key.kmer_front);
            }

            if self.config.verbosity > 2 {
                println!("    Recorded terminators: {} <-> {}", key.kmer_front, key.kmer_back);
            }
        }

        // Check if this group should be flushed
        if self.config.group_flush_threshold > 0 {
            let group_size = self.segment_groups.get(&key).map(|v| v.len()).unwrap_or(0);
            if group_size >= self.config.group_flush_threshold {
                self.flush_group(&key)?;
            }
        }

        Ok(())
    }

    /// Thread-safe version of add_segment_with_kmers for parallel processing
    /// Matches C++ AGC worker thread behavior with fine-grained locking (DashMap)
    #[allow(clippy::too_many_arguments)]
    fn add_segment_with_kmers_threadsafe(
        config: &StreamingCompressorConfig,
        segment_groups: &Arc<DashMap<SegmentGroupKey, Vec<SegmentInfo>>>,
        group_terminators: &Arc<DashMap<u64, Vec<u64>>>,
        total_segments: &Arc<AtomicUsize>,
        sample_name: &str,
        contig_name: &str,
        seg_part_no: usize,
        segment: Contig,
        kmer_front: u64,
        kmer_back: u64,
        is_rev_comp: bool,
    ) -> Result<()> {
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
        let final_is_rev_comp = is_rev_comp ^ needs_rc;

        // Reverse complement if needed
        let final_segment = if needs_rc {
            segment.iter().rev().map(|&b| match b {
                0 => 3, 1 => 2, 2 => 1, 3 => 0, _ => b,
            }).collect()
        } else {
            segment
        };

        // Check if this is a new group (fine-grained locking with DashMap)
        let is_new_group = !segment_groups.contains_key(&key);

        let seg_info = SegmentInfo {
            sample_name: sample_name.to_string(),
            contig_name: contig_name.to_string(),
            seg_part_no,
            data: final_segment,
            is_rev_comp: final_is_rev_comp,
        };

        // DashMap provides per-shard locking - only locks the shard containing this key!
        // This matches C++ AGC's per-group mutex approach (line 220 in agc_compressor.h)
        segment_groups.entry(key.clone()).or_default().push(seg_info);

        // Update terminators if new group (also fine-grained locking via DashMap)
        if is_new_group && key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
            group_terminators.entry(key.kmer_front).or_default().push(key.kmer_back);
            if key.kmer_front != key.kmer_back {
                group_terminators.entry(key.kmer_back).or_default().push(key.kmer_front);
            }
        }

        // Update counter (lock-free atomic increment)
        total_segments.fetch_add(1, Ordering::Relaxed);

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
                // Try compressing the reference with configured level
                let compressed = compress_segment_configured(&ref_segment.data, self.config.compression_level)?;

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
        let compressed_packs: Vec<(Vec<SegmentInfo>, Vec<u8>, Vec<u8>)> = packs
            .par_iter()
            .map(|(seg_infos, packed_data)| {
                let compressed = compress_segment_configured(packed_data, compression_level)
                    .expect("Compression failed");
                (seg_infos.clone(), packed_data.clone(), compressed)
            })
            .collect();

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
        let compressed_packs: Vec<(Vec<SegmentInfo>, Vec<u8>, Vec<u8>)> = packs
            .par_iter()
            .map(|(seg_infos, packed_data)| {
                let compressed = compress_segment_configured(packed_data, compression_level)
                    .expect("Compression failed");
                (seg_infos.clone(), packed_data.clone(), compressed)
            })
            .collect();

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
                self.archive.add_part(stream_id, &compressed_with_marker, packed_data.len() as u64)?;
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

        // DEBUG: Print first 50 groups sorted by group_id to compare with C++
        if self.config.verbosity > 0 {
            println!("\n=== First 50 unique groups (sorted by group_id) ===");
            let mut groups_with_ids: Vec<(u32, u64, u64)> = self.group_metadata
                .iter()
                .map(|(key, meta)| (meta.group_id, key.kmer_front, key.kmer_back))
                .collect();
            groups_with_ids.sort_by_key(|(id, _, _)| *id);

            for (group_id, front_kmer, back_kmer) in groups_with_ids.iter().take(50) {
                println!("  Group {}: front_kmer={}, back_kmer={}", group_id, front_kmer, back_kmer);
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
