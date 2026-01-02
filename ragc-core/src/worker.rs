// Worker thread implementation for streaming compression pipeline
// Matches C++ AGC's start_compressing_threads (agc_compressor.cpp:1097-1270)

use crate::contig_compression::{compress_contig, CompressionContext};
use crate::priority_queue::{BoundedPriorityQueue, PopResult};
use crate::segment_buffer::BufferedSegments;
use crate::segment_compression::compress_reference_segment;
use crate::task::{ContigProcessingStage, Task};
use crate::zstd_pool::compress_segment_pooled;
use ahash::{AHashMap, AHashSet};
use ragc_common::{Archive, Contig};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Barrier, Mutex};

/// Segment placement info for collection metadata
///
/// Matches C++ AGC's segment placement tracking (agc_compressor.cpp:974-1050)
#[derive(Debug, Clone)]
struct SegmentPlacement {
    sample_name: String,
    contig_name: String,
    place: usize, // Segment index within contig (seg_part_no)
    group_id: u32,
    in_group_id: u32,
    is_rev_comp: bool,
    raw_length: u32, // Uncompressed segment length
}

/// Segment group manager matching C++ AGC's CSegment
///
/// Manages a group of segments with the same (kmer1, kmer2) pair:
/// - Stores reference segment (first in group)
/// - Buffers subsequent delta segments
/// - Writes to archive streams
struct SegmentGroup {
    group_id: u32,
    stream_id: usize,          // Delta stream for packed segments
    ref_stream_id: usize,      // Reference stream for first segment
    reference: Option<Contig>, // First segment (reference for LZ encoding)
    ref_written: bool,         // Whether reference has been written
    in_group_counter: u32,     // Counter for in_group_id assignment
}

impl SegmentGroup {
    fn new(group_id: u32, stream_id: usize, ref_stream_id: usize) -> Self {
        SegmentGroup {
            group_id,
            stream_id,
            ref_stream_id,
            reference: None,
            ref_written: false,
            in_group_counter: 0,
        }
    }

    /// Add a segment to this group
    ///
    /// The first segment becomes the reference, subsequent segments are delta-encoded.
    /// Returns in_group_id (0 for reference, 1+ for delta segments)
    fn add_segment(&mut self, seg_data: &[u8], archive: &mut Archive) -> anyhow::Result<u32> {
        if self.reference.is_none() {
            // First segment - store as reference with adaptive compression
            let seg_vec = seg_data.to_vec();
            self.reference = Some(seg_vec.clone());

            // Compress reference segment (chooses between plain ZSTD or tuple packing)
            let (compressed, marker) = compress_reference_segment(&seg_vec)?;

            // Write compressed reference to archive with marker byte as metadata
            archive.add_part(self.ref_stream_id, &compressed, marker as u64)?;
            self.ref_written = true;

            let in_group_id = self.in_group_counter;
            self.in_group_counter += 1;
            Ok(in_group_id) // in_group_id = 0 for reference
        } else {
            // Subsequent segment - delta encode against reference
            // TODO: Implement LZ encoding - for now just compress with ZSTD
            let seg_vec = seg_data.to_vec();
            let compressed = compress_segment_pooled(&seg_vec, 17)?;

            // Write compressed delta to archive (marker 0 = plain ZSTD)
            archive.add_part(self.stream_id, &compressed, 0)?;

            let in_group_id = self.in_group_counter;
            self.in_group_counter += 1;
            Ok(in_group_id)
        }
    }
}

/// Shared state accessible to all worker threads
///
/// This matches C++ AGC's shared variables captured by lambda (agc_compressor.cpp:1099):
/// - Atomic counters (processed_bases, processed_samples)
/// - Mutexes for shared collections (v_raw_contigs, vv_fallback_minimizers, etc.)
/// - Archive writer, collection descriptor, etc.
pub struct SharedCompressorState {
    /// Total bases processed (for progress tracking)
    pub processed_bases: AtomicUsize,

    /// Number of samples processed (for progress tracking)
    pub processed_samples: AtomicUsize,

    /// Contigs that failed to compress well (need adaptive splitters)
    /// Matches C++ AGC's v_raw_contigs (protected by mtx_raw_contigs)
    pub raw_contigs: Mutex<Vec<(String, String, Vec<u8>)>>,

    /// Verbosity level (0 = quiet, 1 = normal, 2 = verbose)
    pub verbosity: usize,

    /// Buffered segments (KNOWN + NEW)
    /// Matches C++ AGC's CBufferedSegPart buffered_seg_part
    pub buffered_segments: Arc<Mutex<BufferedSegments>>,

    /// Splitter k-mers (exact set)
    /// Matches C++ AGC's hs_splitters
    pub splitters: Arc<Mutex<AHashSet<u64>>>,

    /// Bloom filter for splitters (fast probabilistic check)
    /// Matches C++ AGC's bloom_splitters
    pub bloom_splitters: Arc<Mutex<crate::bloom_filter::BloomFilter>>,

    /// Per-thread vectors for accumulating new splitters (adaptive mode)
    /// Matches C++ AGC's vv_splitters (agc_compressor.h:721)
    /// Each worker thread accumulates splitters, then merges at barrier
    pub vv_splitters: Mutex<Vec<Vec<u64>>>,

    /// Reference genome singleton k-mers (for adaptive mode exclusion)
    /// Matches C++ AGC's v_candidate_kmers (agc_compressor.h:708)
    /// Sorted for efficient set_difference operations
    pub v_candidate_kmers: Vec<u64>,

    /// Reference genome duplicated k-mers (for adaptive mode exclusion)
    /// Matches C++ AGC's v_duplicated_kmers (agc_compressor.h:710)
    /// Sorted for efficient set_difference operations
    pub v_duplicated_kmers: Vec<u64>,

    /// K-mer length for segmentation
    pub kmer_length: usize,

    /// Adaptive compression mode (find new splitters for hard contigs)
    pub adaptive_mode: bool,

    /// Map: (kmer1, kmer2) → group_id
    /// Matches C++ AGC's map_segments
    /// Lower group_id wins when same (k1, k2) pair (earlier samples are reference)
    pub map_segments: Arc<Mutex<AHashMap<(u64, u64), u32>>>,

    /// Map: kmer → Vec<kmer> (sorted terminators)
    /// Matches C++ AGC's map_segments_terminators
    /// Used for split detection: find shared terminators between k1 and k2
    pub map_segments_terminators: Arc<Mutex<AHashMap<u64, Vec<u64>>>>,

    /// Concatenated genomes mode (treat all contigs as one sample)
    pub concatenated_genomes: bool,

    /// Total number of segment groups (assigned group IDs)
    /// Matches C++ AGC's no_segments
    pub no_segments: Arc<Mutex<u32>>,

    /// Number of raw groups (groups without LZ encoding)
    /// Matches C++ AGC's no_raw_groups
    pub no_raw_groups: u32,

    /// Archive writer for segment output
    /// Matches C++ AGC's out_archive
    pub archive: Option<Arc<Mutex<Archive>>>,

    /// Segment groups indexed by group_id
    /// Matches C++ AGC's v_segments (vector<CSegment*>)
    /// Lazily created when first segment arrives for that group
    pub v_segments: Arc<Mutex<Vec<Option<SegmentGroup>>>>,

    /// Collection metadata for samples/contigs/segments
    /// Matches C++ AGC's collection_desc
    pub collection: Option<Arc<Mutex<ragc_common::CollectionV3>>>,

    /// Auxiliary queue for adaptive mode (re-enqueue hard contigs)
    /// Matches C++ AGC's pq_contigs_desc_aux
    pub aux_queue: Arc<BoundedPriorityQueue<Task>>,

    /// Working queue pointer (switches between main and aux)
    /// Matches C++ AGC's pq_contigs_desc_working
    pub working_queue: Mutex<Arc<BoundedPriorityQueue<Task>>>,
    // TODO: Add more shared state as needed:
    // - vv_fallback_minimizers: Mutex<Vec<Vec<Vec<u64>>>>
    // - vv_splitters: Mutex<Vec<Vec<u64>>>
}

impl SharedCompressorState {
    /// Create new shared state
    pub fn new(
        verbosity: usize,
        kmer_length: usize,
        adaptive_mode: bool,
        concatenated_genomes: bool,
        no_raw_groups: u32,
        main_queue: Arc<BoundedPriorityQueue<Task>>,
    ) -> Self {
        // Create bloom filter sized for expected splitters
        // C++ AGC uses ~10K splitters, allocate 8 bits/item = 80K bits = 10KB
        let bloom_filter = crate::bloom_filter::BloomFilter::new(80 * 1024);

        // Create auxiliary queue for adaptive mode (unlimited capacity)
        let aux_queue = Arc::new(BoundedPriorityQueue::new(1, usize::MAX));

        SharedCompressorState {
            processed_bases: AtomicUsize::new(0),
            processed_samples: AtomicUsize::new(0),
            raw_contigs: Mutex::new(Vec::new()),
            verbosity,
            buffered_segments: Arc::new(Mutex::new(BufferedSegments::new(no_raw_groups as usize))),
            splitters: Arc::new(Mutex::new(AHashSet::new())),
            bloom_splitters: Arc::new(Mutex::new(bloom_filter)),
            vv_splitters: Mutex::new(Vec::new()),
            v_candidate_kmers: Vec::new(),
            v_duplicated_kmers: Vec::new(),
            kmer_length,
            adaptive_mode,
            map_segments: Arc::new(Mutex::new(AHashMap::new())),
            map_segments_terminators: Arc::new(Mutex::new(AHashMap::new())),
            concatenated_genomes,
            no_segments: Arc::new(Mutex::new(0)),
            no_raw_groups,
            archive: None, // Set later when archive is created
            v_segments: Arc::new(Mutex::new(Vec::new())),
            collection: None, // Set later when collection is created
            aux_queue,
            working_queue: Mutex::new(main_queue), // Start with main queue
        }
    }
}

/// Worker thread main loop
///
/// This matches C++ AGC's worker lambda (agc_compressor.cpp:1099-1270):
/// ```cpp
/// v_threads.emplace_back([&, i, n_t]() {
///     auto zstd_cctx = ZSTD_createCCtx();
///     auto zstd_dctx = ZSTD_createDCtx();
///     uint32_t thread_id = i;
///
///     while(true) {
///         task_t task;
///         auto q_res = pq_contigs_desc_working->PopLarge(task);
///         // ... process task ...
///     }
///
///     ZSTD_freeCCtx(zstd_cctx);
///     ZSTD_freeDCtx(zstd_dctx);
/// });
/// ```
///
/// # Arguments
/// * `worker_id` - Thread ID (0 to num_workers-1)
/// * `queue` - Priority queue for pulling tasks
/// * `barrier` - Synchronization barrier for registration/new_splitters stages
/// * `shared` - Shared state accessible to all workers
pub fn worker_thread(
    worker_id: usize,
    num_workers: usize,
    _queue: Arc<BoundedPriorityQueue<Task>>, // Legacy parameter, use shared.working_queue instead
    barrier: Arc<Barrier>,
    shared: Arc<SharedCompressorState>,
) {
    // TODO: Create per-thread ZSTD contexts
    // let mut zstd_encoder = ...;
    // let mut zstd_decoder = ...;

    loop {
        // Pop task from working queue (can switch between main and aux)
        // Matches C++ AGC's pq_contigs_desc_working (agc_compressor.cpp:1108)
        let queue = {
            let working_queue_guard = shared.working_queue.lock().unwrap();
            Arc::clone(&*working_queue_guard)
        };
        let (result, task_opt) = queue.pop_large();

        match result {
            // Queue empty but producers still active - retry
            PopResult::Empty => continue,

            // Queue empty and no producers - exit worker loop
            PopResult::Completed => break,

            // Successfully got a task - process it
            PopResult::Normal => {
                let task = task_opt.expect("PopResult::Normal should have task");

                match task.stage {
                    ContigProcessingStage::Registration => {
                        // Handle registration synchronization (agc_compressor.cpp:1114-1185)
                        handle_registration_stage(worker_id, &barrier, &shared);
                        continue; // Return to top of loop
                    }

                    ContigProcessingStage::NewSplitters => {
                        // Handle adaptive splitter finding (agc_compressor.cpp:1187-1237)
                        handle_new_splitters_stage(worker_id, num_workers, &barrier, &shared);
                        continue; // Return to top of loop
                    }

                    ContigProcessingStage::AllContigs => {
                        // Preprocess contig before compression (agc_compressor.cpp:1241)
                        // TODO: preprocess_raw_contig(&task.sequence);
                    }

                    ContigProcessingStage::HardContigs => {
                        // Hard contigs don't need preprocessing (already normalized)
                    }
                }

                // Compress contig (AllContigs and HardContigs both use this path)
                let ctg_size = task.sequence.len();

                if compress_contig_task(&task, worker_id, &barrier, &shared) {
                    // Success: update progress counter (agc_compressor.cpp:1248-1255)
                    let old_pb = shared
                        .processed_bases
                        .fetch_add(ctg_size, Ordering::Relaxed);
                    let new_pb = old_pb + ctg_size;

                    // Print progress every 10 MB
                    if shared.verbosity > 0 && old_pb / 10_000_000 != new_pb / 10_000_000 {
                        eprintln!("Compressed: {} Mb\r", new_pb / 1_000_000);
                    }
                } else {
                    // Failed to compress well: save for reprocessing (agc_compressor.cpp:1258-1261)
                    let mut raw_contigs = shared.raw_contigs.lock().unwrap();
                    raw_contigs.push((
                        task.sample_name.clone(),
                        task.contig_name.clone(),
                        task.sequence.clone(),
                    ));
                }

                // Task memory automatically freed (Rust Drop)
            }
        }
    }

    // ZSTD contexts automatically freed (Rust Drop)
}

/// Handle registration stage synchronization
///
/// Matches C++ AGC's registration logic (agc_compressor.cpp:1114-1185)
/// See BARRIER_USAGE_PATTERN.md for detailed documentation.
fn handle_registration_stage(
    worker_id: usize,
    barrier: &Arc<Barrier>,
    shared: &Arc<SharedCompressorState>,
) {
    // Barrier 1: All workers arrive
    let wait_result = barrier.wait();

    // Leader thread: register segments and process fallbacks
    if wait_result.is_leader() {
        register_segments(shared);
        // TODO: process_fallback_minimizers(shared);
    }

    // Barrier 2: Wait for registration complete
    barrier.wait();

    // All threads: store their segments
    store_segments(worker_id, shared);

    // Barrier 3: Wait for storage complete
    barrier.wait();

    // Thread 0 or 1: cleanup and flush (C++ AGC pattern: BOTH threads do work!)
    // See agc_compressor.cpp:1136-1180
    if worker_id == 0 {
        // TODO: buffered_seg_part.clear(max(1, num_workers - 1));
        // TODO: update processed_samples
        // TODO: store_contig_batch if needed
        // TODO: flush archive buffers
        // TODO: if adaptive_compression: switch back to main queue
    } else if worker_id == 1 {
        // TODO: update processed_samples (same as thread 0!)
        // TODO: store_contig_batch if needed (same as thread 0!)
        // TODO: flush archive buffers (same as thread 0!)
    }

    // Barrier 4: All ready to continue
    barrier.wait();
}

/// Find new splitters for a hard contig
///
/// Matches C++ AGC's find_new_splitters (agc_compressor.cpp:2046-2077)
///
/// Algorithm:
/// 1. Extract all k-mers from contig and find singletons
/// 2. Exclude k-mers present in reference genome (singletons)
/// 3. Exclude k-mers present in reference genome (duplicates)
/// 4. Add remaining k-mers to thread-local splitter vector
///
/// # Arguments
/// * `contig` - The contig sequence that needs new splitters
/// * `thread_id` - Worker thread ID for vv_splitters indexing
/// * `shared` - Shared compressor state with reference k-mers
fn find_new_splitters(
    contig: &ragc_common::Contig,
    thread_id: usize,
    shared: &Arc<SharedCompressorState>,
) {
    use crate::kmer_extract::{enumerate_kmers, remove_non_singletons};

    // Step 1: Extract k-mers from contig and find singletons
    let mut v_contig_kmers = enumerate_kmers(contig, shared.kmer_length);
    v_contig_kmers.sort_unstable();
    remove_non_singletons(&mut v_contig_kmers, 0);

    if shared.verbosity > 1 {
        eprintln!(
            "find_new_splitters: contig has {} singleton k-mers",
            v_contig_kmers.len()
        );
    }

    // Step 2: Exclude k-mers in reference genome (singletons)
    // C++ AGC uses v_candidate_kmers_offset to skip some candidates
    // For now, we use the full v_candidate_kmers (offset=0)
    let mut v_tmp = Vec::with_capacity(v_contig_kmers.len());
    set_difference(&v_contig_kmers, &shared.v_candidate_kmers, &mut v_tmp);

    if shared.verbosity > 1 {
        eprintln!(
            "find_new_splitters: {} k-mers after excluding reference singletons",
            v_tmp.len()
        );
    }

    // Step 3: Exclude k-mers in reference genome (duplicates)
    v_contig_kmers.clear();
    set_difference(&v_tmp, &shared.v_duplicated_kmers, &mut v_contig_kmers);

    if shared.verbosity > 1 {
        eprintln!(
            "find_new_splitters: {} NEW splitters found",
            v_contig_kmers.len()
        );
    }

    // Step 4: Add to thread-local splitter vector
    let mut vv_splitters = shared.vv_splitters.lock().unwrap();
    vv_splitters[thread_id].extend(v_contig_kmers);
}

/// Compute set difference: result = a \ b (elements in a but not in b)
///
/// Both input vectors must be sorted. Matches std::set_difference behavior.
fn set_difference(a: &[u64], b: &[u64], result: &mut Vec<u64>) {
    result.clear();
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            result.push(a[i]);
            i += 1;
        } else if a[i] > b[j] {
            j += 1;
        } else {
            // a[i] == b[j], skip
            i += 1;
            j += 1;
        }
    }

    // Add remaining elements from a
    result.extend_from_slice(&a[i..]);
}

/// Handle new splitters stage synchronization
///
/// Matches C++ AGC's new_splitters logic (agc_compressor.cpp:1187-1237)
/// See BARRIER_USAGE_PATTERN.md for detailed documentation.
fn handle_new_splitters_stage(
    worker_id: usize,
    num_workers: usize,
    barrier: &Arc<Barrier>,
    shared: &Arc<SharedCompressorState>,
) {
    // Barrier 1: All workers arrive
    barrier.wait();

    // bloom_insert logic: Insert new splitters into bloom filter and hash set
    // In C++ AGC, this is a lambda defined inline (lines 1191-1209)
    //
    // Multi-threaded execution: All threads merge in parallel
    // Single-threaded: Only thread 0 does the merge
    if num_workers > 1 || worker_id == 0 {
        let mut splitters = shared.splitters.lock().unwrap();
        let mut bloom = shared.bloom_splitters.lock().unwrap();
        let mut vv_splitters = shared.vv_splitters.lock().unwrap();

        let mut total_new = 0;
        for thread_splitters in vv_splitters.iter_mut() {
            total_new += thread_splitters.len();
            for &kmer in thread_splitters.iter() {
                splitters.insert(kmer);
                bloom.insert(kmer);
            }
            thread_splitters.clear();
        }

        if shared.verbosity > 0 && total_new > 0 {
            eprintln!("Adaptive mode: Added {} new splitters", total_new);
            eprintln!("Total splitters: {}", splitters.len());
        }

        // Check bloom filter filling factor and resize if > 0.3 (matches C++ AGC)
        let filling_factor = bloom.filling_factor();
        if filling_factor > 0.3 {
            if shared.verbosity > 1 {
                eprintln!(
                    "Bloom filter filling factor {:.2}, resizing...",
                    filling_factor
                );
            }

            // Resize to accommodate current + expected growth
            let new_size_bits = (splitters.len() as f64 / 0.25) as usize * 8; // 25% target fill
            bloom.resize(new_size_bits);

            // Re-insert all splitters
            for &kmer in splitters.iter() {
                bloom.insert(kmer);
            }

            if shared.verbosity > 1 {
                eprintln!("Bloom filter resized to {} bits", new_size_bits);
            }
        }
    }

    // Thread 0: Re-enqueue hard contigs and switch queues
    // Matches C++ AGC lines 1216-1236
    if worker_id == 0 {
        let mut raw_contigs = shared.raw_contigs.lock().unwrap();

        if shared.verbosity > 0 && !raw_contigs.is_empty() {
            eprintln!(
                "Adaptive mode: Re-enqueueing {} hard contigs for reprocessing",
                raw_contigs.len()
            );
        }

        // Re-enqueue hard contigs into aux_queue with HardContigs stage
        // Priority = 1 (higher than normal contigs), cost = sequence length
        for (sample_name, contig_name, sequence) in raw_contigs.drain(..) {
            let cost = sequence.len();
            shared.aux_queue.emplace(
                Task::new_contig(
                    sample_name,
                    contig_name,
                    sequence,
                    ContigProcessingStage::HardContigs,
                ),
                1, // priority
                cost,
            );
        }

        // Enqueue registration sync tokens for next round
        // Priority = 0 (lowest), cost = 0
        shared.aux_queue.emplace_many_no_cost(
            Task::new_sync(ContigProcessingStage::Registration),
            0, // priority
            num_workers,
        );

        // Switch working queue to aux queue
        // Matches C++ AGC: pq_contigs_desc_working = pq_contigs_desc_aux
        let mut working_queue = shared.working_queue.lock().unwrap();
        *working_queue = Arc::clone(&shared.aux_queue);

        if shared.verbosity > 1 {
            eprintln!("Switched to auxiliary queue for hard contig reprocessing");
        }
    }

    // Barrier 2: Ready to continue with new splitters
    // All workers synchronized after splitter merge
    barrier.wait();
}

/// Compress a contig task
///
/// Matches C++ AGC's compress_contig (agc_compressor.cpp:2000-2054)
///
/// Returns:
/// - true: Contig compressed successfully
/// - false: Contig didn't compress well (needs adaptive splitters)
fn compress_contig_task(
    task: &Task,
    worker_id: usize,
    _barrier: &Arc<Barrier>,
    shared: &Arc<SharedCompressorState>,
) -> bool {
    // Create compression context from shared state
    let ctx = CompressionContext {
        splitters: Arc::clone(&shared.splitters),
        bloom_splitters: Arc::clone(&shared.bloom_splitters),
        buffered_segments: Arc::clone(&shared.buffered_segments),
        kmer_length: shared.kmer_length,
        adaptive_mode: shared.adaptive_mode,
        map_segments: Arc::clone(&shared.map_segments),
        map_segments_terminators: Arc::clone(&shared.map_segments_terminators),
        concatenated_genomes: shared.concatenated_genomes,
    };

    // Call the actual compression function
    let success = compress_contig(&task.sample_name, &task.contig_name, &task.sequence, &ctx);

    // Adaptive mode: Handle hard contigs (C++ AGC agc_compressor.cpp:2033-2039)
    if shared.adaptive_mode
        && !success
        && task.stage == ContigProcessingStage::AllContigs
        && task.sequence.len() >= 1000
    {
        // This is a hard contig - no splitters found, and it's large enough to matter
        if shared.verbosity > 1 {
            eprintln!(
                "Hard contig detected: {}/{} ({} bp)",
                task.sample_name,
                task.contig_name,
                task.sequence.len()
            );
        }

        // Find new splitters from this contig
        find_new_splitters(&task.sequence, worker_id, shared);

        // Store contig for reprocessing after NewSplitters barrier
        let mut raw_contigs = shared.raw_contigs.lock().unwrap();
        raw_contigs.push((
            task.sample_name.clone(),
            task.contig_name.clone(),
            task.sequence.clone(),
        ));

        // Return false to indicate contig wasn't processed yet
        return false;
    }

    success
}

/// Register segments - assign group IDs to NEW segments
///
/// Matches C++ AGC's register_segments (agc_compressor.cpp:954-971)
///
/// Called by leader thread only (thread 0) during registration barrier.
///
/// # Steps
/// 1. Sort KNOWN segments by (sample, contig, part_no)
/// 2. Process NEW segments - assign group IDs
/// 3. Register archive streams (placeholder for now)
/// 4. Distribute raw group segments (if applicable)
/// 5. Restart read pointer for parallel storage
fn register_segments(shared: &Arc<SharedCompressorState>) {
    let mut buffered = shared.buffered_segments.lock().unwrap();

    // Step 1: Sort KNOWN segments in parallel
    // TODO: Implement parallel sort with rayon
    // For now: sequential sort is already done in BufferedSegments
    buffered.sort_known(1); // Single-threaded for now

    // Step 2: Process NEW segments - assign group IDs
    // Pass global map_segments so new segments can find existing groups
    let no_new = buffered.process_new(&shared.map_segments);

    drop(buffered);

    if no_new > 0 {
        // Step 3: Register archive streams for new groups
        let mut no_segments = shared.no_segments.lock().unwrap();
        let current_no_segments = *no_segments;

        if let Some(archive_mutex) = &shared.archive {
            let mut archive = archive_mutex.lock().unwrap();

            // Register streams for each new group
            // Each group gets two streams: "seg-N" (delta) and "seg_dN" (reference)
            for i in 0..no_new {
                let seg_num = current_no_segments + i;
                // Use AGC v3 format (version 3000)
                archive.register_stream(&ragc_common::stream_delta_name(3000, seg_num));
                archive.register_stream(&ragc_common::stream_ref_name(3000, seg_num));
            }
        }

        *no_segments += no_new;
        drop(no_segments);

        // Step 4: Resize v_segments vector to accommodate new groups
        let no_segments_value = *shared.no_segments.lock().unwrap();
        let mut v_segments = shared.v_segments.lock().unwrap();
        v_segments.resize_with(no_segments_value as usize, || None);
    }

    // Step 5: Distribute raw group segments (if applicable)
    if shared.no_raw_groups > 0 {
        let buffered = shared.buffered_segments.lock().unwrap();
        buffered.distribute_segments(0, 0, shared.no_raw_groups);
    }

    // Step 6: Restart read pointer for parallel storage
    let buffered = shared.buffered_segments.lock().unwrap();
    buffered.restart_read_vec();
}

/// Store segments - write all buffered segments to archive
///
/// Matches C++ AGC's store_segments (agc_compressor.cpp:974-1050)
///
/// Called by ALL worker threads in parallel during registration barrier.
///
/// # Algorithm
/// 1. Atomic get_vec_id() to claim a block of groups
/// 2. For each group in block:
///    - get_part() to pop segments
///    - Create CSegment object on first use (lazy initialization)
///    - Update map_segments and map_segments_terminators (CRITICAL!)
///    - Write segment to archive
/// 3. Batch collection metadata updates (every 32 segments)
fn store_segments(_worker_id: usize, shared: &Arc<SharedCompressorState>) {
    const MAX_BUFF_SIZE: usize = 32;
    let mut buffered_placements: Vec<SegmentPlacement> = Vec::with_capacity(MAX_BUFF_SIZE);

    let buffered = shared.buffered_segments.lock().unwrap();

    loop {
        // Step 1: Atomic block allocation
        let block_group_id = buffered.get_vec_id();

        if block_group_id < 0 {
            break; // No more blocks
        }

        // Step 2: Process groups in block (backwards iteration)
        for group_id in
            (block_group_id - crate::segment_buffer::PART_ID_STEP + 1..=block_group_id).rev()
        {
            if buffered.is_empty_part(group_id) {
                continue;
            }

            // Step 3: Process all segments in this group
            while let Some((
                kmer1,
                kmer2,
                sample_name,
                contig_name,
                seg_data,
                is_rev_comp,
                seg_part_no,
            )) = buffered.get_part(group_id)
            {
                // Step 4: Create SegmentGroup on first use (lazy initialization)
                let mut v_segments = shared.v_segments.lock().unwrap();
                if group_id >= 0
                    && (group_id as usize) < v_segments.len()
                    && v_segments[group_id as usize].is_none()
                {
                    // Get stream IDs for this group
                    if let Some(archive_mutex) = &shared.archive {
                        let archive = archive_mutex.lock().unwrap();
                        let stream_id = archive
                            .get_stream_id(&ragc_common::stream_delta_name(3000, group_id as u32))
                            .unwrap_or(0);
                        let ref_stream_id = archive
                            .get_stream_id(&ragc_common::stream_ref_name(3000, group_id as u32))
                            .unwrap_or(0);
                        drop(archive);

                        v_segments[group_id as usize] =
                            Some(SegmentGroup::new(group_id as u32, stream_id, ref_stream_id));
                    }
                }

                // Step 5: Write segment to archive and get in_group_id
                let mut in_group_id = 0;
                if group_id >= 0 && (group_id as usize) < v_segments.len() {
                    if let Some(segment_group) = &mut v_segments[group_id as usize] {
                        if let Some(archive_mutex) = &shared.archive {
                            let mut archive = archive_mutex.lock().unwrap();
                            if let Ok(id) = segment_group.add_segment(&seg_data, &mut archive) {
                                in_group_id = id;
                            }
                        }
                    }
                }
                drop(v_segments);

                // Step 6: Buffer collection metadata
                buffered_placements.push(SegmentPlacement {
                    sample_name,
                    contig_name,
                    place: seg_part_no as usize,
                    group_id: group_id as u32,
                    in_group_id,
                    is_rev_comp,
                    raw_length: seg_data.len() as u32,
                });

                // Flush batch if full
                if buffered_placements.len() >= MAX_BUFF_SIZE {
                    if let Some(collection_mutex) = &shared.collection {
                        let mut collection = collection_mutex.lock().unwrap();
                        for placement in &buffered_placements {
                            let _ = collection.add_segment_placed(
                                &placement.sample_name,
                                &placement.contig_name,
                                placement.place,
                                placement.group_id,
                                placement.in_group_id,
                                placement.is_rev_comp,
                                placement.raw_length,
                            );
                        }
                    }
                    buffered_placements.clear();
                }

                // Step 6: Update map_segments (CRITICAL REGISTRATION!)
                let key = (kmer1, kmer2);
                let mut map_segments = shared.map_segments.lock().unwrap();
                map_segments
                    .entry(key)
                    .and_modify(|existing| {
                        // Keep lower group_id (earlier samples are reference)
                        if (group_id as u32) < *existing {
                            *existing = group_id as u32;
                        }
                    })
                    .or_insert(group_id as u32);
                drop(map_segments);

                // Step 7: Update map_segments_terminators (CRITICAL!)
                if kmer1 != crate::contig_compression::MISSING_KMER
                    && kmer2 != crate::contig_compression::MISSING_KMER
                {
                    let mut terminators = shared.map_segments_terminators.lock().unwrap();

                    // Add k2 to k1's terminator list
                    let list1 = terminators.entry(kmer1).or_insert_with(Vec::new);
                    if !list1.contains(&kmer2) {
                        list1.push(kmer2);
                        list1.sort_unstable();
                    }

                    // Add k1 to k2's terminator list (if different)
                    if kmer1 != kmer2 {
                        let list2 = terminators.entry(kmer2).or_insert_with(Vec::new);
                        if !list2.contains(&kmer1) {
                            list2.push(kmer1);
                            list2.sort_unstable();
                        }
                    }

                    drop(terminators);
                }
            }
        }
    }

    // Step 9: Final flush of remaining placements
    if !buffered_placements.is_empty() {
        if let Some(collection_mutex) = &shared.collection {
            let mut collection = collection_mutex.lock().unwrap();
            for placement in &buffered_placements {
                let _ = collection.add_segment_placed(
                    &placement.sample_name,
                    &placement.contig_name,
                    placement.place,
                    placement.group_id,
                    placement.in_group_id,
                    placement.is_rev_comp,
                    placement.raw_length,
                );
            }
        }
    }
}

// ============================================================================
// Phase 4: Main Compression Flow
// ============================================================================

/// Create an AGC archive from FASTA files - complete end-to-end pipeline
///
/// This is the top-level API for creating AGC archives using the streaming pipeline.
///
/// # Arguments
/// * `output_path` - Path to write the AGC archive
/// * `sample_files` - Vector of (sample_name, file_path) pairs
/// * `splitters` - Set of splitter k-mers (from determine_splitters)
/// * `kmer_length` - K-mer length for segmentation
/// * `segment_size` - Segment size parameter for collection
/// * `num_threads` - Number of threads to use
/// * `adaptive_mode` - Enable adaptive splitter finding
/// * `concatenated_genomes` - Treat all contigs as one sample
/// * `verbosity` - Verbosity level (0=quiet, 1=normal, 2=verbose)
///
/// # Returns
/// Ok(()) on success, Err on failure
pub fn create_agc_archive(
    output_path: &str,
    sample_files: Vec<(String, String)>,
    splitters: AHashSet<u64>,
    candidate_kmers: AHashSet<u64>,
    duplicated_kmers: AHashSet<u64>,
    kmer_length: usize,
    segment_size: u32,
    num_threads: usize,
    adaptive_mode: bool,
    concatenated_genomes: bool,
    verbosity: usize,
) -> anyhow::Result<()> {
    use ragc_common::CollectionV3;

    if sample_files.is_empty() {
        return Ok(());
    }

    let num_samples = sample_files.len();

    // Create archive for writing
    let mut archive = Archive::new_writer();
    archive.open(output_path)?;

    // Write file_type_info stream (C++ AGC compatibility)
    {
        let mut data = Vec::new();
        let append_str = |data: &mut Vec<u8>, s: &str| {
            data.extend_from_slice(s.as_bytes());
            data.push(0);
        };

        append_str(&mut data, "producer");
        append_str(&mut data, "ragc");

        append_str(&mut data, "producer_version_major");
        append_str(&mut data, &ragc_common::AGC_FILE_MAJOR.to_string());

        append_str(&mut data, "producer_version_minor");
        append_str(&mut data, &ragc_common::AGC_FILE_MINOR.to_string());

        append_str(&mut data, "producer_version_build");
        append_str(&mut data, "0");

        append_str(&mut data, "file_version_major");
        append_str(&mut data, &ragc_common::AGC_FILE_MAJOR.to_string());

        append_str(&mut data, "file_version_minor");
        append_str(&mut data, &ragc_common::AGC_FILE_MINOR.to_string());

        append_str(&mut data, "comment");
        append_str(
            &mut data,
            &format!(
                "RAGC v.{}.{}",
                ragc_common::AGC_FILE_MAJOR,
                ragc_common::AGC_FILE_MINOR
            ),
        );

        let stream_id = archive.register_stream("file_type_info");
        archive.add_part(stream_id, &data, 7)?; // 7 key-value pairs

        if verbosity > 0 {
            eprintln!(
                "Wrote file_type_info: version {}.{}",
                ragc_common::AGC_FILE_MAJOR,
                ragc_common::AGC_FILE_MINOR
            );
        }
    }

    // Write params stream (kmer_length, min_match_len, pack_cardinality, segment_size)
    {
        let params_stream_id = archive.register_stream("params");
        let mut params_data = Vec::new();

        // Standard params format (16 bytes for C++ AGC compatibility)
        params_data.extend_from_slice(&(kmer_length as u32).to_le_bytes());
        params_data.extend_from_slice(&20u32.to_le_bytes()); // min_match_len (default)
        params_data.extend_from_slice(&50u32.to_le_bytes()); // pack_cardinality (default)
        params_data.extend_from_slice(&segment_size.to_le_bytes());

        archive.add_part(params_stream_id, &params_data, 0)?;

        if verbosity > 0 {
            eprintln!(
                "Wrote params: k={}, segment_size={}",
                kmer_length, segment_size
            );
        }
    }

    // Create collection for metadata
    let mut collection = CollectionV3::new();
    collection.set_config(segment_size, kmer_length as u32, None);

    // Register collection streams in archive
    collection.prepare_for_compression(&mut archive)?;

    // Wrap in Arc<Mutex<>> for shared access during compression
    let archive = Arc::new(Mutex::new(archive));
    let collection = Arc::new(Mutex::new(collection));

    // Run compression pipeline
    compress_samples_streaming_with_archive(
        sample_files,
        splitters,
        candidate_kmers,
        duplicated_kmers,
        kmer_length,
        num_threads,
        adaptive_mode,
        concatenated_genomes,
        verbosity,
        Some(archive.clone()),
        Some(collection.clone()),
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    // Serialize collection metadata to archive
    {
        let mut archive_guard = archive.lock().unwrap();
        let mut collection_guard = collection.lock().unwrap();

        if verbosity > 0 {
            eprintln!(
                "Serializing collection metadata for {} samples...",
                num_samples
            );
        }

        // Write sample names
        collection_guard.store_batch_sample_names(&mut archive_guard)?;

        // Write contig names and segment details
        collection_guard.store_contig_batch(&mut archive_guard, 0, num_samples)?;

        if verbosity > 0 {
            eprintln!("Collection metadata serialized successfully");
        }
    }

    // Close archive (writes footer)
    let mut archive = archive.lock().unwrap();
    archive.close()?;

    if verbosity > 0 {
        eprintln!("AGC archive created: {}", output_path);
    }

    Ok(())
}

/// Compress sample files using streaming pipeline with priority queue and worker threads
///
/// Internal function - use create_agc_archive() for complete end-to-end pipeline
///
/// Matches C++ AGC's AddSampleFiles() (agc_compressor.cpp:2121-2270)
///
/// # Arguments
/// * `sample_files` - Vector of (sample_name, file_path) pairs
/// * `splitters` - Set of splitter k-mers
/// * `kmer_length` - K-mer length for segmentation
/// * `num_threads` - Number of threads to use
/// * `adaptive_mode` - Enable adaptive splitter finding
/// * `concatenated_genomes` - Treat all contigs as one sample
/// * `verbosity` - Verbosity level (0=quiet, 1=normal, 2=verbose)
/// * `archive` - Optional archive for writing compressed segments
/// * `collection` - Optional collection for metadata tracking
///
/// # Returns
/// Ok(()) on success, Err on failure
fn compress_samples_streaming_with_archive(
    sample_files: Vec<(String, String)>,
    splitters: AHashSet<u64>,
    candidate_kmers: AHashSet<u64>,
    duplicated_kmers: AHashSet<u64>,
    kmer_length: usize,
    num_threads: usize,
    adaptive_mode: bool,
    concatenated_genomes: bool,
    verbosity: usize,
    archive: Option<Arc<Mutex<Archive>>>,
    collection: Option<Arc<Mutex<ragc_common::CollectionV3>>>,
) -> Result<(), String> {
    use crate::genome_io::GenomeIO;

    if sample_files.is_empty() {
        return Ok(());
    }

    // Step 1: Queue initialization (agc_compressor.cpp:2128-2132)
    let queue_capacity = std::cmp::max(2_u64 << 30, num_threads as u64 * (192_u64 << 20));
    let queue = Arc::new(BoundedPriorityQueue::new(1, queue_capacity as usize));

    // Step 2: Worker thread count (agc_compressor.cpp:2134)
    let no_workers = if num_threads < 8 {
        num_threads
    } else {
        num_threads - 1
    };

    // Convert reference k-mers to sorted vectors (for set_difference in adaptive mode)
    let mut v_candidate_kmers: Vec<u64> = candidate_kmers.into_iter().collect();
    v_candidate_kmers.sort_unstable();

    let mut v_duplicated_kmers: Vec<u64> = duplicated_kmers.into_iter().collect();
    v_duplicated_kmers.sort_unstable();

    // Step 3: Shared state initialization
    // Initialize bloom filter and populate with base splitters
    let mut bloom_filter = crate::bloom_filter::BloomFilter::new(
        splitters.len() * 8, // 8 bits per splitter
    );
    for &kmer in &splitters {
        bloom_filter.insert(kmer);
    }

    // Create auxiliary queue for adaptive mode
    let aux_queue = Arc::new(BoundedPriorityQueue::new(1, usize::MAX));

    let shared = Arc::new(SharedCompressorState {
        processed_bases: AtomicUsize::new(0),
        processed_samples: AtomicUsize::new(0),
        raw_contigs: Mutex::new(Vec::new()),
        verbosity,
        buffered_segments: Arc::new(Mutex::new(BufferedSegments::new(0))),
        splitters: Arc::new(Mutex::new(splitters)),
        bloom_splitters: Arc::new(Mutex::new(bloom_filter)),
        vv_splitters: Mutex::new(vec![Vec::new(); no_workers]),
        v_candidate_kmers,
        v_duplicated_kmers,
        kmer_length,
        adaptive_mode,
        map_segments: Arc::new(Mutex::new(AHashMap::new())),
        map_segments_terminators: Arc::new(Mutex::new(AHashMap::new())),
        concatenated_genomes,
        no_segments: Arc::new(Mutex::new(0)),
        no_raw_groups: 0,
        archive,
        v_segments: Arc::new(Mutex::new(Vec::new())),
        collection,
        aux_queue,
        working_queue: Mutex::new(Arc::clone(&queue)), // Start with main queue
    });

    // Step 4: Thread spawning (agc_compressor.cpp:2136-2141)
    let barrier = Arc::new(Barrier::new(no_workers));
    let mut worker_handles = Vec::with_capacity(no_workers);

    for worker_id in 0..no_workers {
        let q = queue.clone();
        let b = barrier.clone();
        let s = shared.clone();

        let handle = std::thread::spawn(move || {
            worker_thread(worker_id, no_workers, q, b, s);
        });

        worker_handles.push(handle);
    }

    // Step 5: Main processing loop (agc_compressor.cpp:2163-2242)
    let mut sample_priority = usize::MAX;
    let mut _cnt_contigs_in_sample = 0;
    const PACK_CARDINALITY: usize = 50; // TODO: Get from config

    for (sample_name, file_path) in sample_files {
        if verbosity > 0 {
            eprintln!("Processing sample: {} from {}", sample_name, file_path);
        }

        // Open FASTA file
        let mut gio = match GenomeIO::open(&file_path) {
            Ok(g) => g,
            Err(e) => {
                eprintln!("Cannot open file {}: {}", file_path, e);
                continue;
            }
        };

        let mut any_contigs_added = false;

        // Read contigs from file
        while let Ok(Some((contig_name, sequence))) = gio.read_contig_raw() {
            if concatenated_genomes {
                // Concatenated mode: treat all contigs as one sample (empty name)
                // Matches C++ AGC behavior: all genomes concatenated into single sample
                let concat_sample_name = String::from("");

                // Register contig with empty sample name
                if let Some(collection_mutex) = &shared.collection {
                    let mut collection = collection_mutex.lock().unwrap();
                    let _ = collection.register_sample_contig(&concat_sample_name, &contig_name);
                }

                let cost = sequence.len();
                queue.emplace(
                    Task::new_contig(
                        concat_sample_name,
                        contig_name.clone(),
                        sequence,
                        ContigProcessingStage::AllContigs,
                    ),
                    sample_priority,
                    cost,
                );

                any_contigs_added = true;
            } else {
                // Normal mode: one sample per file
                // Register sample/contig with collection
                if let Some(collection_mutex) = &shared.collection {
                    let mut collection = collection_mutex.lock().unwrap();
                    let _ = collection.register_sample_contig(&sample_name, &contig_name);
                }

                let cost = sequence.len();
                queue.emplace(
                    Task::new_contig(
                        sample_name.clone(),
                        contig_name.clone(),
                        sequence,
                        ContigProcessingStage::AllContigs,
                    ),
                    sample_priority,
                    cost,
                );

                any_contigs_added = true;
            }
        }

        // Step 6: Send synchronization tokens after each sample (agc_compressor.cpp:2148-2155)
        if !concatenated_genomes && any_contigs_added {
            let sync_stage = if adaptive_mode {
                ContigProcessingStage::NewSplitters
            } else {
                ContigProcessingStage::Registration
            };

            // Send exactly no_workers sync tokens (0 cost)
            queue.emplace_many_no_cost(Task::new_sync(sync_stage), sample_priority, no_workers);

            sample_priority -= 1;
        }
    }

    // Step 7: Signal completion and wait for workers (agc_compressor.cpp:2254-2256)
    queue.mark_completed();

    for handle in worker_handles {
        handle.join().map_err(|_| "Worker thread panicked")?;
    }

    if verbosity > 0 {
        let total_bases = shared.processed_bases.load(Ordering::Relaxed);
        eprintln!("Compression complete: {} bases processed", total_bases);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shared_state_creation() {
        let queue = Arc::new(BoundedPriorityQueue::new(1, 1000));
        let state = SharedCompressorState::new(1, 21, false, false, 0, queue);
        assert_eq!(state.verbosity, 1);
        assert_eq!(state.kmer_length, 21);
        assert_eq!(state.adaptive_mode, false);
        assert_eq!(state.concatenated_genomes, false);
        assert_eq!(state.no_raw_groups, 0);
        assert_eq!(state.processed_bases.load(Ordering::Relaxed), 0);
        assert_eq!(state.processed_samples.load(Ordering::Relaxed), 0);
        assert_eq!(state.raw_contigs.lock().unwrap().len(), 0);
        assert_eq!(*state.no_segments.lock().unwrap(), 0);
    }

    #[test]
    fn test_worker_thread_completion() {
        use std::thread;

        let queue = Arc::new(BoundedPriorityQueue::new(1, 1000));
        let barrier = Arc::new(Barrier::new(2));
        let shared = Arc::new(SharedCompressorState::new(
            0,
            21,
            false,
            false,
            0,
            queue.clone(),
        ));

        // Mark queue as completed (no tasks)
        queue.mark_completed();

        // Spawn two workers
        let handles: Vec<_> = (0..2)
            .map(|worker_id| {
                let q = queue.clone();
                let b = barrier.clone();
                let s = shared.clone();

                thread::spawn(move || {
                    worker_thread(worker_id, 2, q, b, s);
                })
            })
            .collect();

        // Workers should exit immediately
        for handle in handles {
            handle.join().unwrap();
        }
    }

    #[test]
    fn test_worker_thread_with_task() {
        use std::thread;

        let queue = Arc::new(BoundedPriorityQueue::new(1, 1000));
        let barrier = Arc::new(Barrier::new(1));
        let shared = Arc::new(SharedCompressorState::new(
            0,
            21,
            false,
            false,
            0,
            queue.clone(),
        ));

        // Enqueue a simple task
        queue.emplace(
            Task::new_contig(
                "sample1".to_string(),
                "chr1".to_string(),
                vec![0, 1, 2, 3, 0, 1, 2, 3], // 8 bases
                ContigProcessingStage::AllContigs,
            ),
            100,
            8,
        );
        queue.mark_completed();

        // Spawn single worker
        let q = queue.clone();
        let b = barrier.clone();
        let s = shared.clone();

        let handle = thread::spawn(move || {
            worker_thread(0, 1, q, b, s);
        });

        handle.join().unwrap();

        // Verify task was processed (compress_contig_task returns true)
        assert_eq!(shared.processed_bases.load(Ordering::Relaxed), 8);
    }
}
