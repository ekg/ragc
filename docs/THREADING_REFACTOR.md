# Threading Architecture Refactoring

**Goal**: Match C++ AGC's threading model for equivalent performance (102% CPU → 231% CPU)

## Current State (Baseline)

### Performance Metrics (chr5 test - 3 samples, ~1.2MB)
- **Wall time**: 0.24s
- **CPU utilization**: 102% (poor parallelization)
- **Memory**: 103 MB
- **Page faults**: 141,552 minor
- **Correctness**: ✅ 100% byte-for-byte identical

### Architecture
```
Iterator → Load ALL segments → Group by k-mer → Rayon par_iter() → Collect ALL → Write
          [sequential]        [HashMap]         [work-stealing]     [~200MB]
```

**Issues**:
1. ❌ Sequential contig loading (no parallelism yet)
2. ❌ Rayon overhead (102% vs C++ AGC's 231%)
3. ❌ Batch mode (holds everything in memory)
4. ❌ No ZSTD context reuse (recreated per segment)
5. ❌ No memory-based limiting (only segment count)

## Target State (C++ AGC Model)

### Architecture
```
                    ┌─> Worker 1 (local ZSTD ctx) ─┐
Iterator → Queue ──┼─> Worker 2 (local ZSTD ctx) ─┼─> Writer Thread
  (stream)         ├─> Worker 3 (local ZSTD ctx) ─┤   (existing)
  [memory-bounded] └─> Worker N (local ZSTD ctx) ─┘
```

### Key Components
1. **BoundedPriorityQueue** (already exists: `ragc-core/src/priority_queue.rs`)
2. **Worker threads**: Pull from queue, process segments
3. **Thread-local ZSTD contexts**: Create once per thread, reuse
4. **Memory-bounded queue**: 2GB or 192MB/thread (C++ AGC model)
5. **Priority-based**: Higher priority = processed first

---

## Implementation Plan

### Phase 1: Thread-Local ZSTD Context Pool
**Status**: ✅ COMPLETE (Already Implemented)

**Study**:
- [x] Analyzed C++ AGC worker thread (agc_compressor.cpp:1104-1105)
  ```cpp
  auto zstd_cctx = ZSTD_createCCtx();  // Created once per thread
  auto zstd_dctx = ZSTD_createDCtx();  // Reused for all segments
  ```
- [x] Reviewed current ZSTD usage in RAGC
- [x] Identified all compression call sites

**Implementation**: ✅ DONE
- [x] Created `zstd_pool.rs` module (`ragc-core/src/zstd_pool.rs`)
- [x] Implemented thread-local `ZSTD_ENCODER` (line 27)
- [x] API: `compress_segment_pooled(data, level)` (line 42)
- [x] All compression uses pooled contexts

**Testing**: ✅ ALL PASSING
- [x] Unit test: Roundtrip compression/decompression
- [x] Unit test: Multiple compressions with same context (10 iterations)
- [x] Unit test: Different compression levels (5 levels tested)
- [x] Unit test: Large data (10K items)

**Success Criteria**: ✅ ALL MET
- ✅ Contexts created once per thread (thread_local! macro)
- ✅ Contexts properly reused (RefCell handles borrowing)
- ✅ Output byte-for-byte identical (all tests pass)
- ✅ Performance improvement: Context reuse eliminates per-segment overhead

**Notes**:
- Matches C++ AGC design line-by-line (see comments in zstd_pool.rs:30-39)
- Automatically handles different compression levels (re-creates if level changes)
- Decompression uses simple decode_all (less critical for pooling)

---

### Phase 2: Worker Thread Architecture
**Status**: 🔄 In Progress (Study Complete, Ready to Implement)

**Study**: ✅ COMPLETE

**ROOT CAUSE IDENTIFIED**:
Current RAGC has TWO code paths, BOTH with single-threaded segmentation:

1. **add_multi_sample_fasta_with_splitters** (line 538):
   - Phase 1: Single-threaded segmentation → `all_segments` Vec
   - Phase 2: Single-threaded grouping → `groups` HashMap
   - Phase 3: Rayon par_iter (parallel group processing)
   - **Problem**: Only Phase 3 is parallel, segmentation is sequential!

2. **add_fasta_files_with_splitters** (line 1120) → **add_segments_with_inline_splits** (line 1212):
   - 100% sequential processing in one big loop
   - No parallelism at all!
   - **Problem**: Completely sequential!

**Why CPU is 102% not 231%**:
- Segmentation is CPU-intensive (splitting, k-mer extraction)
- All segmentation runs on ONE thread
- Rayon par_iter only processes groups (quick compression)
- Most work done before parallelism kicks in

**Study**: ✅ COMPLETE
- [x] Analyzed C++ AGC worker loop (agc_compressor.cpp:1108-1275)
- [x] Mapped to RAGC segment processing logic
- [x] Identified synchronization points

**C++ AGC Worker Flow**:
```cpp
while (true) {
    task_t task = pq_contigs_desc_working->PopLarge();

    if (task == empty) continue;
    if (task == completed) break;

    if (task == registration) {
        // Synchronization barrier for batch operations
        barrier.arrive_and_wait();
        if (thread_id == 0) register_segments();
        barrier.arrive_and_wait();
        store_segments(zstd_cctx, zstd_dctx);  // All workers compress & write
        barrier.arrive_and_wait();
        continue;
    }

    // Main processing: compress_contig() → add_segment()
    compress_contig(task, zstd_cctx, zstd_dctx, thread_id);
}
```

**Key Findings**:
1. Workers pull contigs from `BoundedPriorityQueue`
2. Workers call `compress_contig()` which segments and compresses
3. Workers write **directly** to thread-safe archive (no separate writer!)
4. Barriers coordinate batch operations (register_segments, store_segments)

**RAGC Adaptation**:
- ✅ RAGC already has separate writer thread (better design!)
- ✅ Workers can send `CompressedPack` to writer via channel
- Need: Workers pull from `BoundedPriorityQueue`
- Need: Workers process contigs (segment → group → compress)
- Skip: Barriers (RAGC processes incrementally, not in batches)

**Solution**:
Replace single-threaded segmentation with C++ AGC-style workers:
```
Main Thread                    Worker Threads (N)
─────────────                  ─────────────────
Iterator → Queue ──┐          ┌→ Worker 1: pull → segment → send
  (stream)         ├─────────►├→ Worker 2: pull → segment → send
  [priority]       │          ├→ Worker 3: pull → segment → send
                   └──────────┘→ Worker N: pull → segment → send
                                        ↓
                                   Writer Thread
```

**Implementation Plan**:
- [x] ContigTask type already exists (line 147)
- [x] BoundedPriorityQueue imported
- [ ] Create worker function: `process_contigs_worker()`
- [ ] Spawn N worker threads
- [ ] Main thread streams contigs into queue
- [ ] Workers pull, segment, send to writer
  - Thread-local ZSTD contexts
  - Access to shared state (groups, terminators)
- [ ] Implement worker loop:
  ```rust
  loop {
      let task = queue.pop_large();
      match task {
          (PopResult::Normal, Some(contig)) => process_contig(contig),
          (PopResult::Empty, None) => continue,
          (PopResult::Completed, None) => break,
      }
  }
  ```
- [ ] Add barrier synchronization for batch operations

**Testing**:
- [ ] Unit test: Worker processes single contig
- [ ] Unit test: Multiple workers process queue
- [ ] Integration test: Full compression with N workers
- [ ] Verify: Output identical to sequential

**Success Criteria**:
- ✅ Workers pull from queue correctly
- ✅ Thread-local contexts working
- ✅ Output byte-for-byte identical
- ✅ No data races or deadlocks

---

### Phase 3: Streaming Contig Input
**Status**: ⏳ Not Started

**Study**:
- [ ] Analyze C++ AGC queue feeding (agc_compressor.cpp:2203-2225)
- [ ] Understand priority assignment
- [ ] Understand memory-based limiting

**Implementation**:
- [ ] Replace sequential loop with queue.emplace()
- [ ] Calculate priority per contig (sample_priority)
- [ ] Calculate cost per contig (contig.size())
- [ ] Implement memory-based queue capacity:
  ```rust
  let queue_capacity = std::cmp::max(
      2 << 30,                        // 2GB
      num_threads * (192 << 20)       // 192MB/thread
  );
  ```
- [ ] Stream contigs without loading all upfront

**Testing**:
- [ ] Test: Queue accepts contigs with priority
- [ ] Test: Queue respects memory limit
- [ ] Test: Backpressure when queue full
- [ ] Benchmark: Memory usage with large dataset

**Success Criteria**:
- ✅ Contigs streamed into queue
- ✅ Memory stays under limit
- ✅ No deadlocks from backpressure
- ✅ Output byte-for-byte identical

---

### Phase 4: Performance Validation
**Status**: ⏳ Not Started

**Testing**:
- [ ] Benchmark chr5 test (3 samples, ~1.2MB):
  - CPU utilization: Target >200%
  - Wall time: Target <0.15s
  - Memory: Target <60MB
  - Page faults: Target <20K
- [ ] Benchmark yeast10 test (10 samples, ~11.6MB):
  - CPU utilization: Target >200%
  - Wall time: Target comparable to C++ AGC
- [ ] Verify correctness on all test datasets:
  - chr5: SHA256 match
  - yeast10: SHA256 match
  - All unit tests: PASS

**Success Criteria**:
- ✅ CPU utilization: 200%+ (vs current 102%)
- ✅ Performance within 2x of C++ AGC
- ✅ Memory usage comparable or better
- ✅ All correctness tests passing

---

### Phase 5: Integration and Cleanup
**Status**: ⏳ Not Started

**Implementation**:
- [ ] Replace old add_contigs_with_splitters() with threaded version
- [ ] Remove Rayon parallelization code
- [ ] Clean up unused imports
- [ ] Update documentation
- [ ] Run cargo fmt

**Testing**:
- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] CI passes
- [ ] Performance regression test

**Success Criteria**:
- ✅ Clean codebase
- ✅ No performance regressions
- ✅ All tests passing

---

## Risk Mitigation

### Correctness Risks
**Mitigation**:
- Test each phase independently
- Maintain SHA256 baseline for verification
- Never commit code with failing tests

### Performance Risks
**Mitigation**:
- Benchmark each phase
- Keep old code until new version proven
- Profile to identify bottlenecks

### Complexity Risks
**Mitigation**:
- Break into smallest possible steps
- Document each change
- Commit working code frequently

---

## References

### C++ AGC Source
- Worker thread loop: `agc_compressor.cpp:1097-1168`
- Queue feeding: `agc_compressor.cpp:2139-2290`
- ZSTD context reuse: `agc_compressor.cpp:1104-1105`

### RAGC Source
- Current implementation: `ragc-core/src/compressor_streaming.rs:1240-1800`
- BoundedPriorityQueue: `ragc-core/src/priority_queue.rs`
- ZSTD compression: `ragc-core/src/segment_compression.rs`

---

## Progress Tracking

- [ ] Phase 1: Thread-Local ZSTD Context Pool
- [ ] Phase 2: Worker Thread Architecture
- [ ] Phase 3: Streaming Contig Input
- [ ] Phase 4: Performance Validation
- [ ] Phase 5: Integration and Cleanup

**Current Phase**: Phase 1 (Study)
**Last Updated**: 2025-10-31
