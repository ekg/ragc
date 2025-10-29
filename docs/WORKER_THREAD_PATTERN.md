# Worker Thread Pattern in C++ AGC

**Purpose**: Document C++ AGC's worker thread architecture for Rust implementation

**C++ Reference**:
- Worker loop: `agc_compressor.cpp` lines 1097-1270 (start_compressing_threads)
- compress_contig: `agc_compressor.cpp` lines 833-1093

---

## C++ AGC Worker Thread Structure

### Thread Creation (agc_compressor.cpp:1097-1103)

```cpp
for (uint32_t i = 0; i < n_t; ++i) {
    v_threads.emplace_back([&, i, n_t]() {
        auto zstd_cctx = ZSTD_createCCtx();  // Per-thread compression context
        auto zstd_dctx = ZSTD_createDCtx();  // Per-thread decompression context
        uint32_t thread_id = i;

        // Worker loop...

        ZSTD_freeCCtx(zstd_cctx);
        ZSTD_freeDCtx(zstd_dctx);
    });
}
```

**Key Features**:
- Each worker has unique thread_id (0 to n_t-1)
- Per-thread ZSTD contexts (created once, reused throughout)
- Lambda captures shared state by reference (&)

---

## Worker Loop Pattern

### Main Loop (agc_compressor.cpp:1104-1266)

```cpp
while(true) {
    task_t task;

    // 1. Pop task from priority queue
    auto q_res = pq_contigs_desc_working->PopLarge(task);

    // 2. Handle queue state
    if (q_res == CBoundedPQueue<task_t>::result_t::empty)
        continue;  // Queue empty but producers active, retry
    else if (q_res == CBoundedPQueue<task_t>::result_t::completed)
        break;     // Queue empty and no producers, exit

    // 3. Process task based on stage
    if (get<0>(task) == contig_processing_stage_t::registration) {
        // Registration barrier logic (lines 1115-1185)
        // ... (4 barriers, see BARRIER_USAGE_PATTERN.md)
        continue;
    }

    if (get<0>(task) == contig_processing_stage_t::new_splitters) {
        // Adaptive splitter finding (lines 1187-1237)
        // ... (2 barriers)
        continue;
    }

    if (get<0>(task) == contig_processing_stage_t::all_contigs) {
        preprocess_raw_contig(get<3>(task));
    }

    // 4. Compress contig (all_contigs and hard_contigs)
    size_t ctg_size = get<3>(task).size();

    if (compress_contig(get<0>(task), get<1>(task), get<2>(task),
                       get<3>(task), zstd_cctx, zstd_dctx, thread_id, bar)) {
        // Success: update progress
        auto old_pb = processed_bases.fetch_add(ctg_size);
        auto new_pb = old_pb + ctg_size;

        if (verbosity > 0 && is_app_mode && old_pb / 10'000'000 != new_pb / 10'000'000) {
            cerr << "Compressed: " + to_string(processed_bases / 1'000'000) + " Mb\r";
        }
    } else {
        // Failed to compress well: save for reprocessing with new splitters
        lock_guard<mutex> lck(mtx_raw_contigs);
        v_raw_contigs.emplace_back(get<1>(task), get<2>(task), move(get<3>(task)));
    }

    // 5. Free task memory
    get<3>(task).clear();
    get<3>(task).shrink_to_fit();
}
```

---

## Task Processing Stages

### 1. Registration Stage (lines 1114-1185)

**Purpose**: Synchronize workers to register segments and update terminators

**Pattern**:
```cpp
// Barrier 1: All workers arrive
bar.arrive_and_wait();

// Thread 0: Register segments
if (thread_id == 0)
    register_segments(n_t);

// Thread 0: Process fallback minimizers
if (thread_id == 0)
    for (auto& v_fallback_minimizers : vv_fallback_minimizers) {
        for (auto& x : v_fallback_minimizers)
            add_fallback_mapping(x[0], x[1], x[2], (bool)x[3]);
        v_fallback_minimizers.clear();
    }

// Barrier 2: Wait for registration
bar.arrive_and_wait();

// All threads: Store segments
store_segments(zstd_cctx, zstd_dctx);

// Barrier 3: Wait for storage
bar.arrive_and_wait();

// Thread 0 or 1: Update progress and flush
if (thread_id == 0) {
    buffered_seg_part.clear(max(1u, n_t-1));
    ++processed_samples;  // (simplified)
    out_archive->FlushOutBuffers();

    if (adaptive_compression)
        pq_contigs_desc_working = pq_contigs_desc;  // Switch back to main queue
}
else if (thread_id == 1) {
    ++processed_samples;  // (simplified)
    out_archive->FlushOutBuffers();
}

// Barrier 4: All ready to continue
bar.arrive_and_wait();

continue;  // Return to top of loop
```

### 2. NewSplitters Stage (lines 1187-1237)

**Purpose**: Find adaptive splitters for hard contigs that didn't compress well

**Pattern**:
```cpp
// Barrier 1: All workers arrive
bar.arrive_and_wait();

// Thread 0 or 1: Insert new splitters
auto bloom_insert = [&] {
    for (auto& v : vv_splitters) {
        for (auto& x : v) {
            hs_splitters.insert_fast(x);
            bloom_splitters.insert(x);
        }
        v.clear();
    }

    if (bloom_splitters.filling_factor() > 0.3) {
        bloom_splitters.resize((uint64_t)(hs_splitters.size() / 0.25));
        bloom_splitters.insert(hs_splitters.begin(), hs_splitters.end());
    }
};

if (thread_id == 0) {
    if (n_t == 1)
        bloom_insert();

    // Re-enqueue hard contigs with new splitters
    for (auto& x : v_raw_contigs) {
        auto cost = get<2>(x).size();
        pq_contigs_desc_aux->EmplaceNoLock(
            make_tuple(hard_contigs, get<0>(x), get<1>(x), move(get<2>(x))),
            1, cost);
    }
    v_raw_contigs.clear();

    // Enqueue registration tokens for next round
    pq_contigs_desc_aux->EmplaceManyNoCost(
        make_tuple(registration, "", "", contig_t()), 0, n_t);

    // Switch to aux queue
    pq_contigs_desc_working = pq_contigs_desc_aux;
}
else if (thread_id == 1) {
    bloom_insert();
}

// Barrier 2: All ready to process hard contigs
bar.arrive_and_wait();

continue;  // Return to top of loop
```

### 3. AllContigs Stage (lines 1239-1256)

**Purpose**: Initial contig processing

**Pattern**:
```cpp
if (get<0>(task) == all_contigs) {
    preprocess_raw_contig(get<3>(task));  // Normalize sequence
}

size_t ctg_size = get<3>(task).size();

if (compress_contig(stage, sample_name, contig_name, sequence,
                   zstd_cctx, zstd_dctx, thread_id, bar)) {
    // Success: update processed_bases
    processed_bases.fetch_add(ctg_size);
} else {
    // Failed: save for reprocessing
    lock_guard<mutex> lck(mtx_raw_contigs);
    v_raw_contigs.emplace_back(sample_name, contig_name, move(sequence));
}
```

### 4. HardContigs Stage

**No special handling** - falls through to same compress_contig() path as AllContigs.

**Difference**: Hard contigs are processed with expanded splitter set from adaptive mode.

---

## compress_contig Return Value

```cpp
bool compress_contig(
    contig_processing_stage_t stage,
    const string& sample_name,
    const string& contig_name,
    contig_t& sequence,
    ZSTD_CCtx* zstd_cctx,
    ZSTD_DCtx* zstd_dctx,
    uint32_t thread_id,
    my_barrier& bar
)
```

**Returns**:
- `true`: Contig compressed successfully
  - Segments were small enough (good compression)
  - Can proceed to next task

- `false`: Contig didn't compress well
  - Segments too large (poor compression ratio)
  - Caller saves to v_raw_contigs for adaptive reprocessing
  - Only happens when adaptive_compression enabled

**Key Logic** (agc_compressor.cpp:1082-1093):
```cpp
// Check if segments are acceptable
bool ok = true;
if (adaptive_compression && stage == all_contigs) {
    // ... check segment sizes ...
    if (any segment > threshold)
        ok = false;
}

if (!ok) {
    // Rollback: don't add segments to buffer
    return false;
}

// Add segments to buffered_seg_part
// ...
return true;
```

---

## Rust Implementation Strategy

### Worker Thread Function

```rust
fn worker_thread(
    worker_id: usize,
    queue: Arc<BoundedPriorityQueue<Task>>,
    barrier: Arc<Barrier>,
    shared: Arc<SharedCompressorState>,
) {
    // Create per-thread ZSTD contexts
    let mut zstd_encoder = ZstdEncoder::new();
    let mut zstd_decoder = ZstdDecoder::new();

    loop {
        // Pop task from queue
        let (result, task) = queue.pop_large();

        match result {
            PopResult::Empty => continue,
            PopResult::Completed => break,
            PopResult::Normal => {
                let task = task.unwrap();

                match task.stage {
                    ContigProcessingStage::Registration => {
                        handle_registration(worker_id, &barrier, &shared);
                        continue;
                    }

                    ContigProcessingStage::NewSplitters => {
                        handle_new_splitters(worker_id, &barrier, &shared, &queue);
                        continue;
                    }

                    ContigProcessingStage::AllContigs => {
                        preprocess_contig(&task.sequence);
                    }

                    ContigProcessingStage::HardContigs => {
                        // No preprocessing for hard contigs
                    }
                }

                // Compress contig (AllContigs and HardContigs)
                let ctg_size = task.sequence.len();

                if compress_contig(
                    &task,
                    &mut zstd_encoder,
                    &mut zstd_decoder,
                    worker_id,
                    &barrier,
                    &shared,
                ) {
                    // Success: update progress
                    let old_pb = shared.processed_bases.fetch_add(ctg_size, Ordering::Relaxed);
                    let new_pb = old_pb + ctg_size;

                    if shared.verbosity > 0 && old_pb / 10_000_000 != new_pb / 10_000_000 {
                        eprintln!("Compressed: {} Mb\r", new_pb / 1_000_000);
                    }
                } else {
                    // Failed: save for reprocessing
                    let mut raw_contigs = shared.raw_contigs.lock().unwrap();
                    raw_contigs.push((task.sample_name, task.contig_name, task.sequence));
                }
            }
        }
    }

    // ZSTD contexts automatically dropped (Rust RAII)
}
```

### Key Differences from C++

1. **ZSTD contexts**: Rust uses wrapper types with Drop, no manual free needed
2. **Barrier leader**: Use `barrier.wait().is_leader()` instead of `thread_id == 0`
3. **Mutex locking**: Use `lock().unwrap()` instead of lock_guard
4. **Atomic updates**: Same as C++ (`fetch_add` with `Ordering`)

---

## Testing Strategy

### Unit Test for Worker Loop

```rust
#[test]
fn test_worker_loop_basic() {
    use std::sync::{Arc, Barrier};
    use std::thread;

    let queue = BoundedPriorityQueue::new(1, 1000);
    let barrier = Arc::new(Barrier::new(2));
    let shared = Arc::new(create_test_shared_state());

    // Enqueue test tasks
    queue.emplace(
        Task::new_contig("sample1".into(), "chr1".into(), vec![0, 1, 2, 3],
                        ContigProcessingStage::AllContigs),
        100,
        4
    );
    queue.mark_completed();  // Signal no more tasks

    // Spawn workers
    let handles: Vec<_> = (0..2).map(|worker_id| {
        let q = queue.clone();
        let b = barrier.clone();
        let s = shared.clone();

        thread::spawn(move || {
            worker_thread(worker_id, q, b, s);
        })
    }).collect();

    // Wait for completion
    for handle in handles {
        handle.join().unwrap();
    }

    // Verify results
    assert_eq!(shared.processed_bases.load(Ordering::Relaxed), 4);
}
```

---

## Implementation Status

- [C] **Studied**: C++ AGC worker thread pattern ✓
- [R] **Implementation**: Next step - create worker.rs module
- [✓] **Verified**: Will verify with unit tests

**Note**: This is the core of the streaming architecture. All other phases support this worker loop.
