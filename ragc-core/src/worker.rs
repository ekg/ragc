// Worker thread implementation for streaming compression pipeline
// Matches C++ AGC's start_compressing_threads (agc_compressor.cpp:1097-1270)

use crate::contig_compression::{compress_contig, CompressionContext};
use crate::priority_queue::{BoundedPriorityQueue, PopResult};
use crate::segment_buffer::BufferedSegments;
use crate::task::{ContigProcessingStage, Task};
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Barrier, Mutex};

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
    pub splitters: Arc<Mutex<HashSet<u64>>>,

    /// Bloom filter for splitters (fast probabilistic check)
    /// Matches C++ AGC's bloom_splitters
    /// TODO: Replace () with actual bloom filter implementation
    pub bloom_splitters: Arc<Mutex<()>>,

    /// K-mer length for segmentation
    pub kmer_length: usize,

    /// Adaptive compression mode (find new splitters for hard contigs)
    pub adaptive_mode: bool,

    /// Map: (kmer1, kmer2) → group_id
    /// Matches C++ AGC's map_segments
    /// Lower group_id wins when same (k1, k2) pair (earlier samples are reference)
    pub map_segments: Arc<Mutex<HashMap<(u64, u64), u32>>>,

    /// Map: kmer → Vec<kmer> (sorted terminators)
    /// Matches C++ AGC's map_segments_terminators
    /// Used for split detection: find shared terminators between k1 and k2
    pub map_segments_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>>,

    /// Concatenated genomes mode (treat all contigs as one sample)
    pub concatenated_genomes: bool,

    /// Total number of segment groups (assigned group IDs)
    /// Matches C++ AGC's no_segments
    pub no_segments: Arc<Mutex<u32>>,

    /// Number of raw groups (groups without LZ encoding)
    /// Matches C++ AGC's no_raw_groups
    pub no_raw_groups: u32,

    // TODO: Add more shared state as needed:
    // - vv_fallback_minimizers: Mutex<Vec<Vec<Vec<u64>>>>
    // - vv_splitters: Mutex<Vec<Vec<u64>>>
    // - out_archive: Mutex<ArchiveWriter>
    // - collection_desc: Mutex<Collection>
    // - v_segments: Mutex<Vec<Option<Arc<Segment>>>>
}

impl SharedCompressorState {
    /// Create new shared state
    pub fn new(
        verbosity: usize,
        kmer_length: usize,
        adaptive_mode: bool,
        concatenated_genomes: bool,
        no_raw_groups: u32,
    ) -> Self {
        SharedCompressorState {
            processed_bases: AtomicUsize::new(0),
            processed_samples: AtomicUsize::new(0),
            raw_contigs: Mutex::new(Vec::new()),
            verbosity,
            buffered_segments: Arc::new(Mutex::new(BufferedSegments::new(no_raw_groups as usize))),
            splitters: Arc::new(Mutex::new(HashSet::new())),
            bloom_splitters: Arc::new(Mutex::new(())),
            kmer_length,
            adaptive_mode,
            map_segments: Arc::new(Mutex::new(HashMap::new())),
            map_segments_terminators: Arc::new(Mutex::new(HashMap::new())),
            concatenated_genomes,
            no_segments: Arc::new(Mutex::new(0)),
            no_raw_groups,
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
    queue: Arc<BoundedPriorityQueue<Task>>,
    barrier: Arc<Barrier>,
    shared: Arc<SharedCompressorState>,
) {
    // TODO: Create per-thread ZSTD contexts
    // let mut zstd_encoder = ...;
    // let mut zstd_decoder = ...;

    loop {
        // Pop task from priority queue (agc_compressor.cpp:1108)
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
    let bloom_insert = || {
        // TODO: Insert splitters from vv_splitters into hs_splitters and bloom_splitters
        // TODO: Clear vv_splitters
        // TODO: If bloom_splitters.filling_factor() > 0.3: resize and rebuild
    };

    // Thread 0: Re-enqueue hard contigs and switch queues
    if worker_id == 0 {
        // Single-threaded: thread 0 does bloom_insert
        if num_workers == 1 {
            bloom_insert();
        }

        // Re-enqueue hard contigs for reprocessing with new splitters
        // TODO: for (sample, contig, seq) in raw_contigs:
        //     aux_queue.emplace(Task::new_contig(sample, contig, seq, HardContigs), priority=1, cost=seq.len())
        // TODO: raw_contigs.clear()

        // Enqueue registration sync tokens for next round
        // TODO: aux_queue.emplace_many_no_cost(Task::new_sync(Registration), priority=0, n_items=num_workers)

        // Switch working queue to aux queue
        // TODO: working_queue = aux_queue
    }
    // Thread 1: Handle bloom_insert for multi-threaded case
    else if worker_id == 1 {
        bloom_insert();
    }

    // Barrier 2: All ready to process hard contigs
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
    _worker_id: usize,
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
    compress_contig(&task.sample_name, &task.contig_name, &task.sequence, &ctx)
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
    let no_new = buffered.process_new();

    drop(buffered);

    if no_new > 0 {
        // Step 3: Register archive streams for new groups
        // TODO: Call out_archive->RegisterStreams() for each new group
        // For now, just update the segment counter

        let mut no_segments = shared.no_segments.lock().unwrap();
        let current_no_segments = *no_segments;

        // Register streams: for i in 0..no_new:
        //   out_archive.register_streams(
        //     format!("seg-{}", current_no_segments + i),
        //     format!("seg-{}", current_no_segments + i)
        //   )

        *no_segments += no_new;
        drop(no_segments);

        // Step 4: Resize v_segments vector
        // TODO: Resize v_segments to accommodate new groups
        // let no_segments_value = *shared.no_segments.lock().unwrap();
        // shared.v_segments.lock().unwrap().resize(no_segments_value, None);
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
    // const MAX_BUFF_SIZE: usize = 32;
    // TODO: Implement buffered placement tracking
    // let mut buffered_placements: Vec<SegmentPlacement> = Vec::with_capacity(MAX_BUFF_SIZE);

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
            while let Some((kmer1, kmer2, _sample_name, _contig_name, _seg_data, _is_rev_comp, _seg_part_no)) = buffered.get_part(group_id) {
                // Step 4: Create CSegment on first use (lazy initialization)
                // TODO: Create v_segments[group_id] if None
                // let mut v_segments = shared.v_segments.lock().unwrap();
                // if v_segments[group_id].is_none() {
                //     v_segments[group_id] = Some(Arc::new(Mutex::new(Segment::new(...))));
                // }

                // Step 5: Update map_segments (CRITICAL REGISTRATION!)
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

                // Step 6: Update map_segments_terminators (CRITICAL!)
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

                // Step 7: Write segment to archive
                // TODO: segment.add() or segment.add_raw()
                // let in_group_id = if group_id < shared.no_raw_groups as i32 {
                //     segment.add_raw(&part.seg_data)
                // } else {
                //     segment.add(&part.seg_data)
                // };

                // Step 8: Buffer collection metadata
                // TODO: buffered_placements.push(SegmentPlacement { ... });
                // if buffered_placements.len() == MAX_BUFF_SIZE {
                //     collection.add_segments_placed(&buffered_placements);
                //     buffered_placements.clear();
                // }
            }
        }
    }

    // Step 9: Final flush
    // if !buffered_placements.is_empty() {
    //     collection.add_segments_placed(&buffered_placements);
    // }
}

// ============================================================================
// Phase 4: Main Compression Flow
// ============================================================================

/// Compress sample files using streaming pipeline with priority queue and worker threads
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
///
/// # Returns
/// Ok(()) on success, Err on failure
pub fn compress_samples_streaming(
    sample_files: Vec<(String, String)>,
    splitters: HashSet<u64>,
    kmer_length: usize,
    num_threads: usize,
    adaptive_mode: bool,
    concatenated_genomes: bool,
    verbosity: usize,
) -> Result<(), String> {
    use crate::genome_io::GenomeIO;

    if sample_files.is_empty() {
        return Ok(());
    }

    // Step 1: Queue initialization (agc_compressor.cpp:2128-2132)
    let queue_capacity = std::cmp::max(2_u64 << 30, num_threads as u64 * (192_u64 << 20));
    let queue = Arc::new(BoundedPriorityQueue::new(1, queue_capacity as usize));

    // TODO: Auxiliary queue for adaptive mode
    // let aux_queue = Arc::new(BoundedPriorityQueue::new(1, usize::MAX));

    // Step 2: Worker thread count (agc_compressor.cpp:2134)
    let no_workers = if num_threads < 8 { num_threads } else { num_threads - 1 };

    // Step 3: Shared state initialization
    let shared = Arc::new(SharedCompressorState {
        processed_bases: AtomicUsize::new(0),
        processed_samples: AtomicUsize::new(0),
        raw_contigs: Mutex::new(Vec::new()),
        verbosity,
        buffered_segments: Arc::new(Mutex::new(BufferedSegments::new(0))),
        splitters: Arc::new(Mutex::new(splitters)),
        bloom_splitters: Arc::new(Mutex::new(())),
        kmer_length,
        adaptive_mode,
        map_segments: Arc::new(Mutex::new(HashMap::new())),
        map_segments_terminators: Arc::new(Mutex::new(HashMap::new())),
        concatenated_genomes,
        no_segments: Arc::new(Mutex::new(0)),
        no_raw_groups: 0,
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
    const PACK_CARDINALITY: usize = 50;  // TODO: Get from config

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
                // TODO: Register with empty sample name
                // TODO: Implement concatenated mode logic
                eprintln!("Concatenated genomes mode not yet implemented");
            } else {
                // Normal mode: one sample per file
                // TODO: Register sample/contig with collection

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
            queue.emplace_many_no_cost(
                Task::new_sync(sync_stage),
                sample_priority,
                no_workers,
            );

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
        let state = SharedCompressorState::new(1, 21, false, false, 0);
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
        let shared = Arc::new(SharedCompressorState::new(0, 21, false, false, 0));

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
        let shared = Arc::new(SharedCompressorState::new(0, 21, false, false, 0));

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
