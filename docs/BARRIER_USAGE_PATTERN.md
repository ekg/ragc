# Barrier Synchronization Pattern in C++ AGC

**Purpose**: Document C++ AGC's barrier usage for Rust implementation

**C++ Reference**:
- Definition: `utils_adv.h` lines 105-174 (CAtomicBarrierWithIncrementing)
- Usage: `agc_compressor.cpp` lines 1114-1180 (registration stage)

---

## C++ AGC Barrier Implementation

### Type Definition (agc_compressor.h:603)
```cpp
using my_barrier = CAtomicBarrierWithIncrementing;
```

### CAtomicBarrierWithIncrementing (utils_adv.h:105-174)

```cpp
class CAtomicBarrierWithIncrementing {
    std::atomic<int32_t> a_count;
    std::atomic<int32_t> a_generation;
    int32_t count_reset_value;

public:
    explicit CAtomicBarrierWithIncrementing(int32_t count);

    void arrive_and_wait() {
        int32_t old_generation = a_generation.load();

        if (!a_count.fetch_sub(1, memory_order_relaxed)) {
            // Last thread to arrive - reset for next use
            a_count = count_reset_value;
            ++a_generation;
            a_generation.notify_all();
            return;
        }

        // Wait for generation to change
        a_generation.wait(old_generation);
    }

    // Additional methods for dynamic thread adjustment
    bool try_increment(int32_t inc = 1);
    int32_t try_increment_max(int32_t inc_req);
    void decrement(int32_t dec = 1);
};
```

**Key Features**:
- Reusable barrier (generation counter increments)
- Wait using atomic wait/notify (C++20 feature)
- Can dynamically add/remove threads (try_increment/decrement)

---

## Usage Pattern in Worker Threads

### Registration Stage (agc_compressor.cpp:1114-1180)

```cpp
if (get<0>(task) == contig_processing_stage_t::registration) {
    // 1. BARRIER: All workers arrive
    bar.arrive_and_wait();

    // 2. CRITICAL SECTION: Only thread 0 executes
    if (thread_id == 0)
        register_segments(n_t);

    // 3. CRITICAL SECTION: Only thread 0 processes fallbacks
    if (thread_id == 0)
        for (auto& v_fallback_minimizers : vv_fallback_minimizers) {
            for (auto& x : v_fallback_minimizers)
                add_fallback_mapping(x[0], x[1], x[2], (bool)x[3]);
            v_fallback_minimizers.clear();
        }

    // 4. BARRIER: Wait for registration complete
    bar.arrive_and_wait();

    // 5. PARALLEL: All threads store their segments
    store_segments(zstd_cctx, zstd_dctx);

    // 6. BARRIER: Wait for all storage complete
    bar.arrive_and_wait();

    // 7. CLEANUP: Thread 0 or 1 does cleanup
    if (thread_id == 0) {
        buffered_seg_part.clear(max(1u, n_t-1));
        // ... progress tracking, flushing
    }
    else if (thread_id == 1) {
        // ... progress tracking, flushing
    }

    // 8. BARRIER: Final sync before continuing
    bar.arrive_and_wait();

    continue;  // Process next task
}
```

### New Splitters Stage (agc_compressor.cpp:1187-1240)

```cpp
if (get<0>(task) == contig_processing_stage_t::new_splitters) {
    // 1. BARRIER: All workers arrive
    bar.arrive_and_wait();

    // 2. CRITICAL SECTION: Thread 0 or 1 inserts new splitters
    auto bloom_insert = [&] {
        for (auto& v : vv_splitters) {
            for (auto& x : v) {
                hs_splitters.insert_fast(x);
                bloom_splitters.insert(x);
            }
            v.clear();
        }

        if (bloom_splitters.filling_factor() > 0.3) {
            bloom_splitters.resize(...);
            bloom_splitters.insert(hs_splitters.begin(), hs_splitters.end());
        }
    };

    if (thread_id == 0) {
        if (n_t == 1)
            bloom_insert();

        // Re-enqueue hard contigs
        for (auto& x : v_raw_contigs) {
            auto cost = get<2>(x).size();
            pq_contigs_desc_aux->EmplaceNoLock(
                make_tuple(hard_contigs, get<0>(x), get<1>(x), move(get<2>(x))),
                1, cost);
        }

        v_raw_contigs.clear();

        // Enqueue sync tokens
        pq_contigs_desc_aux->EmplaceManyNoCost(
            make_tuple(registration, "", "", contig_t()),
            0, n_t);

        // Switch to aux queue
        pq_contigs_desc_working = pq_contigs_desc_aux;
    }
    else if (thread_id == 1) {
        bloom_insert();
    }

    // 3. BARRIER: All ready to process hard contigs
    bar.arrive_and_wait();

    continue;
}
```

---

## Pattern Summary

**Common Sequence**:
1. `bar.arrive_and_wait()` - All threads arrive
2. Thread 0 (or sometimes thread 0 OR 1) executes critical section
3. `bar.arrive_and_wait()` - Wait for critical section complete
4. All threads execute parallel work (optional)
5. `bar.arrive_and_wait()` - Wait for parallel work complete
6. Thread 0/1 does cleanup
7. `bar.arrive_and_wait()` - Final sync
8. Continue to next task

**Key Observations**:
- Barriers are reusable (generation-based)
- Only thread 0 (or 0/1) executes single-threaded critical sections
- All threads participate in parallel sections
- Multiple barriers used to separate phases of work
- Barriers created once with N workers: `my_barrier bar(no_workers);`

---

## Rust Implementation Strategy

### Use std::sync::Barrier

Rust's standard library provides a reusable barrier:

```rust
use std::sync::{Arc, Barrier};

// Create barrier for N workers
let barrier = Arc::new(Barrier::new(num_workers));

// In worker thread
let result = barrier.wait();

// Check if this is the leader thread
if result.is_leader() {
    // Execute critical section (only one thread)
}

// All threads wait here
barrier.wait();
```

**Advantages**:
- `Barrier::wait()` returns `BarrierWaitResult`
- `is_leader()` identifies ONE thread to execute critical sections
- Automatically reusable (like C++ version)
- Simpler API than C++ (no manual generation tracking)

**Mapping**:
- C++ `bar.arrive_and_wait()` → Rust `barrier.wait()`
- C++ `if (thread_id == 0)` → Rust `if result.is_leader()`

### Example Worker Code

```rust
fn worker_thread(
    worker_id: usize,
    queue: Arc<BoundedPriorityQueue<Task>>,
    barrier: Arc<Barrier>,
    shared: Arc<SharedState>,
) {
    loop {
        let (result, task) = queue.pop_large();

        match result {
            PopResult::Completed => break,
            PopResult::Empty => continue,
            PopResult::Normal => {
                let task = task.unwrap();

                match task.stage {
                    ContigProcessingStage::Registration => {
                        // 1. All arrive
                        let wait_result = barrier.wait();

                        // 2. Leader executes critical section
                        if wait_result.is_leader() {
                            register_segments(&shared);
                            process_fallback_minimizers(&shared);
                        }

                        // 3. Wait for registration
                        barrier.wait();

                        // 4. All store segments
                        store_segments(worker_id, &shared);

                        // 5. Wait for storage
                        let wait_result = barrier.wait();

                        // 6. Leader cleanup
                        if wait_result.is_leader() {
                            flush_and_update_progress(&shared);
                        }

                        // 7. Final sync
                        barrier.wait();
                    }
                    // ... other stages
                }
            }
        }
    }
}
```

---

## Testing Strategy

### Unit Test for Barrier

```rust
#[test]
fn test_barrier_pattern() {
    use std::sync::{Arc, Barrier};
    use std::thread;
    use std::sync::atomic::{AtomicUsize, Ordering};

    let barrier = Arc::new(Barrier::new(4));
    let counter = Arc::new(AtomicUsize::new(0));
    let leader_count = Arc::new(AtomicUsize::new(0));

    let mut handles = vec![];

    for _ in 0..4 {
        let b = Arc::clone(&barrier);
        let c = Arc::clone(&counter);
        let l = Arc::clone(&leader_count);

        let handle = thread::spawn(move || {
            // All increment counter
            c.fetch_add(1, Ordering::SeqCst);

            // Wait at barrier
            let result = b.wait();

            // Only leader increments leader_count
            if result.is_leader() {
                l.fetch_add(1, Ordering::SeqCst);
            }

            // Wait at second barrier
            b.wait();
        });

        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    assert_eq!(counter.load(Ordering::SeqCst), 4);
    assert_eq!(leader_count.load(Ordering::SeqCst), 1);
}
```

---

## Implementation Status

- [C] **Studied**: C++ AGC barrier pattern ✓
- [R] **Implementation**: Will use `std::sync::Barrier` in worker threads (Phase 2)
- [✓] **Verified**: Documentation complete, test pattern defined

**Note**: No separate barrier module needed - Rust stdlib provides everything required.
Actual usage will be implemented in Phase 2 (Worker Thread Architecture).
