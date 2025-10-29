// Worker thread implementation for streaming compression pipeline
// Matches C++ AGC's start_compressing_threads (agc_compressor.cpp:1097-1270)

use crate::priority_queue::{BoundedPriorityQueue, PopResult};
use crate::task::{ContigProcessingStage, Task};
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

    // TODO: Add more shared state as needed:
    // - buffered_seg_part: BufferedSegments
    // - vv_fallback_minimizers: Mutex<Vec<Vec<Vec<u64>>>>
    // - vv_splitters: Mutex<Vec<Vec<u64>>>
    // - hs_splitters: Mutex<HashSet<u64>>
    // - bloom_splitters: Mutex<BloomFilter>
    // - out_archive: Mutex<ArchiveWriter>
    // - collection_desc: Mutex<Collection>
}

impl SharedCompressorState {
    /// Create new shared state
    pub fn new(verbosity: usize) -> Self {
        SharedCompressorState {
            processed_bases: AtomicUsize::new(0),
            processed_samples: AtomicUsize::new(0),
            raw_contigs: Mutex::new(Vec::new()),
            verbosity,
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
        // TODO: register_segments(shared);
        // TODO: process_fallback_minimizers(shared);
    }

    // Barrier 2: Wait for registration complete
    barrier.wait();

    // All threads: store their segments
    // TODO: store_segments(worker_id, shared);

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
/// Matches C++ AGC's compress_contig (agc_compressor.cpp:833-1093)
///
/// Returns:
/// - true: Contig compressed successfully
/// - false: Contig didn't compress well (needs adaptive splitters)
fn compress_contig_task(
    task: &Task,
    worker_id: usize,
    barrier: &Arc<Barrier>,
    shared: &Arc<SharedCompressorState>,
) -> bool {
    // TODO: Implement actual compression logic
    // This will involve:
    // 1. Split contig at splitters
    // 2. LZ-encode each segment
    // 3. Check segment sizes (if adaptive mode)
    // 4. ZSTD compress segments
    // 5. Add to buffered_seg_part
    //
    // For now, return true (success) as placeholder
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shared_state_creation() {
        let state = SharedCompressorState::new(1);
        assert_eq!(state.verbosity, 1);
        assert_eq!(state.processed_bases.load(Ordering::Relaxed), 0);
        assert_eq!(state.processed_samples.load(Ordering::Relaxed), 0);
        assert_eq!(state.raw_contigs.lock().unwrap().len(), 0);
    }

    #[test]
    fn test_worker_thread_completion() {
        use std::thread;

        let queue = Arc::new(BoundedPriorityQueue::new(1, 1000));
        let barrier = Arc::new(Barrier::new(2));
        let shared = Arc::new(SharedCompressorState::new(0));

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
        let shared = Arc::new(SharedCompressorState::new(0));

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
