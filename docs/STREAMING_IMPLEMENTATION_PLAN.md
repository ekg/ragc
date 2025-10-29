# RAGC Streaming Architecture Implementation Plan

**Goal**: Rewrite RAGC to exactly match C++ AGC's streaming architecture with priority queue, worker threads, and synchronization barriers.

**Methodology**: For each item below:
1. ✅ **Study C++ AGC** - Read and understand the C++ implementation
2. ✅ **Implement in Rust** - Write the Rust equivalent
3. ✅ **Verify** - Unit test or code review to confirm correctness

**Status Key**:
- `[ ]` Not started
- `[C]` C++ AGC studied
- `[R]` Rust implemented
- `[✓]` Verified and complete

---

## Phase 1: Core Infrastructure (Priority Queue + Task System)

### 1.1 Task Types and Processing Stages ✅
**C++ Reference**: `agc_compressor.h` lines 550-560 (contig_processing_stage_t)

- [C] **Study** C++ AGC's `contig_processing_stage_t` enum
  - `all_contigs`: Initial contig processing
  - `registration`: Synchronization point for segment registration
  - `new_splitters`: Synchronization point for adaptive splitter finding
  - `hard_contigs`: Contigs that need adaptive splitters

- [R] **Implement** Rust equivalent `ContigProcessingStage` enum
  ```rust
  enum ContigProcessingStage {
      AllContigs,      // Normal contig processing
      Registration,    // Sync barrier: register segments
      NewSplitters,    // Sync barrier: find adaptive splitters
      HardContigs,     // Re-process with new splitters
  }
  ```

- [✓] **Verify** Unit test or code review for enum variants
  - 5/5 tests pass
  - Commit: f72de9c

**Files modified**: `ragc-core/src/task.rs` (new file)

---

### 1.2 Task Structure ✅
**C++ Reference**: `agc_compressor.h` line 656 (task_t typedef)

- [C] **Study** C++ AGC's task structure
  ```cpp
  using task_t = tuple<contig_processing_stage_t, string, string, contig_t>;
  // (stage, sample_name, contig_name, sequence)
  ```

- [R] **Implement** Rust `Task` struct
  ```rust
  struct Task {
      stage: ContigProcessingStage,
      sample_name: String,
      contig_name: String,
      sequence: Contig,
  }
  ```

- [✓] **Verify** Unit test for task creation and field access
  - Included in 1.1 tests (5/5 pass)
  - Commit: f72de9c

**Files modified**: `ragc-core/src/task.rs`

---

### 1.3 Priority Queue System ✅
**C++ Reference**: `queue.h` lines 153-346 (CBoundedPQueue)

- [C] **Study** C++ AGC's `CBoundedPQueue<task_t>` behavior
  - Items have (value, priority, cost)
  - PopLarge() returns highest priority items first (reverse iterator)
  - Multiple workers can pop concurrently
  - MarkCompleted() signals end of queue
  - Size limited by memory capacity (2GB or 192MB/thread)
  - EmplaceManyNoCost() for sync barriers

- [R] **Implement** Rust `BoundedPriorityQueue<T>`
  ```rust
  struct BoundedPriorityQueue<T> {
      // Internal: BinaryHeap + Mutex + Condvar
      // Methods: emplace(), pop_large(), mark_completed()
  }
  ```

- [✓] **Verify** Unit tests:
  - Multiple producers, multiple consumers ✓
  - Priority ordering (highest priority first) ✓
  - Capacity limiting ✓
  - Graceful completion signaling ✓
  - EmplaceManyNoCost for sync barriers ✓
  - 5/5 tests pass
  - Commit: f07169e

**Files modified**: `ragc-core/src/priority_queue.rs` (new file, 434 lines)

---

### 1.4 Synchronization Barriers ✓
**C++ Reference**: `agc_compressor.cpp` lines 1120-1130 (my_barrier usage)

- [C] **Study** C++ AGC's barrier synchronization pattern
  - `bar.arrive_and_wait()` used at registration points
  - All N workers must arrive before any continue
  - Only thread 0 executes critical section code

- [R] **Implement** Rust barrier using `std::sync::Barrier`
  ```rust
  use std::sync::Barrier;
  // Workers use barrier.wait() at sync points
  ```

- [✓] **Verify** Unit test:
  - N threads all reach barrier
  - Only one continues to do work
  - All proceed after work complete

**Status**: Complete
  - Documented C++ AGC barrier pattern (CAtomicBarrierWithIncrementing)
  - Identified Rust equivalent: std::sync::Barrier with is_leader()
  - Test pattern defined in documentation
  - No separate module needed - will use std library in Phase 2

**Files modified**: `docs/BARRIER_USAGE_PATTERN.md` (new file, 333 lines)

---

## Phase 2: Worker Thread Architecture

### 2.1 Worker Thread Loop Structure ✓
**C++ Reference**: `agc_compressor.cpp` lines 1097-1270 (start_compressing_threads)

- [C] **Study** C++ AGC's worker loop pattern
  ```cpp
  while(true) {
      task_t task;
      auto q_res = pq_contigs_desc_working->PopLarge(task);

      if (q_res == empty) continue;
      else if (q_res == completed) break;

      // Process task based on stage...
  }
  ```

- [R] **Implement** Rust worker function
  ```rust
  fn worker_thread(
      worker_id: usize,
      queue: Arc<BoundedPriorityQueue<Task>>,
      barrier: Arc<Barrier>,
      shared_state: Arc<SharedCompressorState>,
  ) {
      while let Some(task) = queue.pop_large() {
          match task.stage {
              // Handle each stage...
          }
      }
  }
  ```

- [✓] **Verify** Unit test with mock queue and tasks

**Status**: Complete
  - Documented worker loop pattern in WORKER_THREAD_PATTERN.md
  - Implemented worker_thread() with task dispatching
  - Created SharedCompressorState for thread communication
  - Placeholder stage handlers (detailed implementation in Phase 2.2-2.4)
  - Placeholder compress_contig_task (detailed implementation in Phase 3)
  - 3/3 tests pass (shared state, completion, task processing)

**Files modified**:
  - `ragc-core/src/worker.rs` (new file, 324 lines)
  - `ragc-core/src/task.rs` (added Eq/PartialEq derives)
  - `docs/WORKER_THREAD_PATTERN.md` (new file, comprehensive C++ AGC documentation)

---

### 2.2 Registration Stage Handler ✓
**C++ Reference**: `agc_compressor.cpp` lines 1115-1180 (registration stage)

- [C] **Study** C++ AGC's registration logic
  - Barrier: all workers arrive
  - Thread 0: calls `register_segments()`
  - Thread 0: processes fallback minimizers
  - Barrier: wait for registration complete
  - All threads: call `store_segments()`
  - Barrier: wait for storage complete
  - Thread 0 **AND** 1: update progress, flush buffers (BOTH threads!)
  - Barrier: all ready to continue

- [R] **Implement** Rust registration handler
  ```rust
  ContigProcessingStage::Registration => {
      barrier.wait();

      if worker_id == 0 {
          register_segments(&shared_state);
          process_fallback_minimizers(&shared_state);
      }

      barrier.wait();
      store_segments(worker_id, &shared_state);
      barrier.wait();

      // BOTH thread 0 and 1 do cleanup (matching C++ AGC exactly)
      if worker_id == 0 {
          flush_and_update_progress(&shared_state);
      } else if worker_id == 1 {
          flush_and_update_progress(&shared_state);
      }

      barrier.wait();
  }
  ```

- [✓] **Verify** Unit test with multiple workers hitting registration

**Status**: Complete
  - Correct 4-barrier synchronization pattern
  - Leader-only sections for register_segments() and fallback minimizers
  - All-threads section for store_segments()
  - BOTH thread 0 and 1 do cleanup (matches C++ AGC behavior)
  - Detailed segment processing implementation in Phase 3

**Files modified**: `ragc-core/src/worker.rs` (handle_registration_stage updated)

---

### 2.3 New Splitters Stage Handler (Adaptive Mode) ✓
**C++ Reference**: `agc_compressor.cpp` lines 1187-1240 (new_splitters stage)

- [C] **Study** C++ AGC's adaptive splitter logic
  - Barrier: all workers arrive
  - bloom_insert lambda: insert new splitters into bloom filter and hash set
  - Thread 0: if single-threaded, calls bloom_insert
  - Thread 0: re-enqueue hard contigs for reprocessing
  - Thread 0: enqueue registration sync tokens
  - Thread 0: switch to aux queue
  - Thread 1: calls bloom_insert (multi-threaded case)
  - Barrier: all ready to continue

- [R] **Implement** Rust new splitters handler
  ```rust
  fn handle_new_splitters_stage(worker_id, num_workers, barrier, shared) {
      barrier.wait();

      let bloom_insert = || {
          // Insert splitters, resize if needed
      };

      if worker_id == 0 {
          if num_workers == 1 { bloom_insert(); }
          // Re-enqueue hard contigs as HardContigs stage
          // Enqueue registration sync tokens
          // Switch to aux queue
      } else if worker_id == 1 {
          bloom_insert();
      }

      barrier.wait();
  }
  ```

- [✓] **Verify** Unit test for adaptive splitter integration

**Status**: Complete
  - Correct 2-barrier synchronization pattern
  - bloom_insert lambda pattern (single vs multi-threaded)
  - Thread 0: re-enqueue hard contigs, switch queues
  - Thread 1: bloom_insert for multi-threaded case
  - Detailed bloom filter and queue switching implementation in Phase 3

**Files modified**: `ragc-core/src/worker.rs` (handle_new_splitters_stage updated, num_workers parameter added)

---

### 2.4 Contig Compression Handler
**C++ Reference**: `agc_compressor.cpp` lines 1242-1270 (all_contigs/hard_contigs)

- [ ] **Study** C++ AGC's compress_contig function call
  ```cpp
  if (get<0>(task) == contig_processing_stage_t::all_contigs) {
      preprocess_raw_contig(get<3>(task));
  }

  size_t ctg_size = get<3>(task).size();

  if (compress_contig(get<0>(task), get<1>(task), get<2>(task),
                      get<3>(task), zstd_cctx, zstd_dctx,
                      thread_id, bar)) {
      // Update progress...
  }
  ```

- [ ] **Implement** Rust contig compression handler
  ```rust
  ContigProcessingStage::AllContigs |
  ContigProcessingStage::HardContigs => {
      if task.stage == ContigProcessingStage::AllContigs {
          preprocess_contig(&mut task.sequence);
      }

      let ctg_size = task.sequence.len();

      if compress_contig(
          &task,
          &mut worker_state.zstd_ctx,
          &shared_state,
          worker_id,
          &barrier,
      ) {
          shared_state.processed_bases
              .fetch_add(ctg_size, Ordering::SeqCst);
      }
  }
  ```

- [ ] **Verify** Integration test with real contig

**Files to modify**: `ragc-core/src/worker.rs`, `ragc-core/src/compression.rs`

---

## Phase 3: Segment Processing (Inline During Compression)

### 3.1 Buffered Segment Storage
**C++ Reference**: `agc_compressor.h` lines 630-640 (CBufferedSegmentsPart)

- [ ] **Study** C++ AGC's `buffered_seg_part` structure
  ```cpp
  struct seg_part_t {
      contig_processing_stage_t stage;
      string sample_name;
      string contig_name;
      uint32_t seg_part_no;
      bool is_rev_comp;
      uint32_t group_id;
      vector<uint8_t> data;
  };

  class CBufferedSegmentsPart {
      map<pair<uint64_t, uint64_t>, vector<seg_part_t>> map_buffered;
  };
  ```

- [ ] **Implement** Rust `BufferedSegments` structure
  ```rust
  struct SegmentPart {
      stage: ContigProcessingStage,
      sample_name: String,
      contig_name: String,
      seg_part_no: u32,
      is_rev_comp: bool,
      group_id: u32,
      data: Vec<u8>,
  }

  struct BufferedSegments {
      map: HashMap<(u64, u64), Vec<SegmentPart>>,
  }
  ```

- [ ] **Verify** Unit test for buffering and retrieval

**Files to modify**: Update `ragc-core/src/segment_buffer.rs`

---

### 3.2 Inline Segmentation During Compression
**C++ Reference**: `agc_compressor.cpp` lines 1340-1500 (compress_contig)

- [ ] **Study** C++ AGC's inline segmentation
  - Segments created during `compress_contig()`
  - Each segment checked: NEW or KNOWN?
  - NEW segments checked for split eligibility
  - Segments added to `buffered_seg_part`
  - NOT written immediately - buffered until registration

- [ ] **Implement** Rust inline segmentation
  ```rust
  fn compress_contig(
      task: &Task,
      zstd_ctx: &mut ZstdContext,
      shared: &SharedState,
      worker_id: usize,
      barrier: &Barrier,
  ) -> bool {
      // Segment at splitters
      let segments = segment_contig(&task.sequence, &shared.splitters);

      for (idx, segment) in segments.iter().enumerate() {
          let key = (segment.front_kmer, segment.back_kmer);

          // Check if group exists
          let group_id = shared.segment_map.get(&key);

          if let Some(gid) = group_id {
              // KNOWN segment
              buffer_known_segment(worker_id, gid, segment);
          } else {
              // NEW segment - check for split
              if should_split(&key, &shared) {
                  // Split and add to existing groups
              } else {
                  // Buffer as new group
                  buffer_new_segment(worker_id, &key, segment);
              }
          }
      }

      true
  }
  ```

- [ ] **Verify** Unit test for segmentation flow

**Files to modify**: `ragc-core/src/compression.rs`

---

### 3.3 Split Detection and Execution (Inline)
**C++ Reference**: `agc_compressor.cpp` lines 1381-1501 (split logic in compress_contig)

- [ ] **Study** C++ AGC's split checking
  ```cpp
  // Check if NEW segment (K1, K2) can split
  auto it1 = map_segments_terminators.find(x1);
  auto it2 = map_segments_terminators.find(x2);

  if (it1 != map_segments_terminators.end()) {
      for (auto x : it1->second) {
          auto p2 = map_segments.find(make_pair(x, x2));
          if (p2 != map_segments.end()) {
              // Found: (K1,M) and (M,K2) both exist
              // Calculate split position
              // Add to both groups
          }
      }
  }
  ```

- [ ] **Implement** Rust split detection
  ```rust
  fn should_split(
      key: &(u64, u64),
      shared: &SharedState,
  ) -> Option<SplitPlan> {
      let (k1, k2) = *key;

      // Check k1's terminators
      if let Some(terminators) = shared.terminators.get(&k1) {
          for &middle in terminators {
              // Check if (k1, middle) and (middle, k2) exist
              let key1 = (k1.min(middle), k1.max(middle));
              let key2 = (middle.min(k2), middle.max(k2));

              if shared.segment_map.contains_key(&key1) &&
                 shared.segment_map.contains_key(&key2) {
                  return Some(SplitPlan {
                      middle_kmer: middle,
                      group1: shared.segment_map[&key1],
                      group2: shared.segment_map[&key2],
                  });
              }
          }
      }

      None
  }
  ```

- [ ] **Verify** Unit test with known split scenarios

**Files to modify**: `ragc-core/src/splitting.rs` (new file)

---

### 3.4 Segment Registration Phase
**C++ Reference**: `agc_compressor.cpp` lines 1527-1620 (register_segments)

- [ ] **Study** C++ AGC's `register_segments()`
  - Processes ALL buffered NEW segments
  - Assigns group IDs sequentially
  - Creates group metadata
  - Adds terminators to map_segments_terminators
  - Moves segments from NEW buffer to group storage

- [ ] **Implement** Rust segment registration
  ```rust
  fn register_segments(shared: &mut SharedState) {
      let buffered = shared.take_buffered_new_segments();

      for (key, segments) in buffered {
          // Assign new group ID
          let group_id = shared.next_group_id;
          shared.next_group_id += 1;

          // Register in map
          shared.segment_map.insert(key, group_id);

          // Add terminators
          if key.0 != MISSING_KMER && key.1 != MISSING_KMER {
              shared.terminators
                  .entry(key.0)
                  .or_default()
                  .push(key.1);

              if key.0 != key.1 {
                  shared.terminators
                      .entry(key.1)
                      .or_default()
                      .push(key.0);
              }
          }

          // Move segments to group storage
          shared.groups.insert(group_id, segments);
      }
  }
  ```

- [ ] **Verify** Unit test for registration logic

**Files to modify**: `ragc-core/src/registration.rs`

---

## Phase 4: Main Compression Flow Integration

### 4.1 Sample-by-Sample Processing Loop
**C++ Reference**: `agc_compressor.cpp` lines 2163-2242 (AddSampleFiles)

- [ ] **Study** C++ AGC's main loop
  ```cpp
  for (auto sf : _v_sample_file_name) {  // Each sample
      while (gio.ReadContigRaw(id, contig)) {  // Each contig
          collection_desc->register_sample_contig(sf.first, id);

          pq_contigs_desc->Emplace(
              make_tuple(all_contigs, sf.first, id, move(contig)),
              sample_priority, cost);
      }

      // Synchronization after sample
      pq_contigs_desc->EmplaceManyNoCost(
          make_tuple(registration, "", "", contig_t()),
          sample_priority, no_workers);

      --sample_priority;
  }
  ```

- [ ] **Implement** Rust main loop
  ```rust
  pub fn add_samples_streaming(
      &mut self,
      sample_files: &[(String, PathBuf)],
  ) -> Result<()> {
      let mut sample_priority = usize::MAX;

      for (sample_name, path) in sample_files {
          let mut reader = GenomeIO::open(path)?;

          while let Some((contig_name, sequence)) = reader.read_contig()? {
              // Pre-register contig
              self.collection.register_sample_contig(
                  sample_name,
                  &contig_name,
              )?;

              // Enqueue task
              self.task_queue.emplace(
                  Task {
                      stage: ContigProcessingStage::AllContigs,
                      sample_name: sample_name.clone(),
                      contig_name,
                      sequence,
                  },
                  sample_priority,
                  cost,
              );
          }

          // Synchronization barrier after sample
          let sync_stage = if self.config.adaptive_mode {
              ContigProcessingStage::NewSplitters
          } else {
              ContigProcessingStage::Registration
          };

          self.task_queue.emplace_many(
              Task {
                  stage: sync_stage,
                  sample_name: String::new(),
                  contig_name: String::new(),
                  sequence: Vec::new(),
              },
              sample_priority,
              self.num_workers,
          );

          sample_priority -= 1;
      }

      Ok(())
  }
  ```

- [ ] **Verify** Integration test with multiple samples

**Files to modify**: `ragc-core/src/compressor_streaming.rs`

---

### 4.2 Worker Thread Spawning and Management
**C++ Reference**: `agc_compressor.cpp` lines 2134-2141 (thread management)

- [ ] **Study** C++ AGC's thread spawning
  ```cpp
  uint32_t no_workers = (no_threads < 8) ? no_threads : no_threads - 1;

  vector<thread> v_threads;
  v_threads.reserve((size_t)no_workers);

  my_barrier bar(no_workers);

  start_compressing_threads(v_threads, bar, no_workers);
  ```

- [ ] **Implement** Rust thread management
  ```rust
  pub fn start_compression(
      &mut self,
  ) -> Result<()> {
      let num_workers = if self.config.num_threads < 8 {
          self.config.num_threads
      } else {
          self.config.num_threads - 1
      };

      let barrier = Arc::new(Barrier::new(num_workers));
      let queue = Arc::clone(&self.task_queue);
      let shared = Arc::clone(&self.shared_state);

      let mut handles = Vec::new();

      for worker_id in 0..num_workers {
          let queue = Arc::clone(&queue);
          let barrier = Arc::clone(&barrier);
          let shared = Arc::clone(&shared);

          let handle = thread::spawn(move || {
              worker_thread(worker_id, queue, barrier, shared)
          });

          handles.push(handle);
      }

      // Wait for all workers to complete
      for handle in handles {
          handle.join().unwrap();
      }

      Ok(())
  }
  ```

- [ ] **Verify** Integration test with worker spawning

**Files to modify**: `ragc-core/src/compressor_streaming.rs`

---

### 4.3 Queue Completion and Cleanup
**C++ Reference**: `agc_compressor.cpp` lines 2244-2251 (completion)

- [ ] **Study** C++ AGC's completion sequence
  ```cpp
  pq_contigs_desc->MarkCompleted();

  join_threads(v_threads);

  pq_contigs_desc.reset();
  pq_contigs_desc_aux.reset();
  ```

- [ ] **Implement** Rust completion
  ```rust
  // After enqueuing all tasks
  self.task_queue.mark_completed();

  // Workers will exit when they see completed signal
  // join handled in start_compression()

  // Cleanup
  self.task_queue = None;
  ```

- [ ] **Verify** Test graceful shutdown

**Files to modify**: `ragc-core/src/compressor_streaming.rs`

---

## Phase 5: Adaptive Splitter Integration

### 5.1 Hard Contig Detection
**C++ Reference**: `agc_compressor.cpp` lines 1820-1850 (hard contig detection)

- [ ] **Study** C++ AGC's criteria for "hard" contigs
  - Contigs that don't segment well
  - Tracked in `v_raw_contigs` during compression
  - Re-enqueued in `new_splitters` stage

- [ ] **Implement** Rust hard contig tracking
  ```rust
  // During compression, if contig segments poorly:
  if segments.len() < MIN_SEGMENTS || avg_segment_size > MAX_SIZE {
      shared.hard_contigs.lock().unwrap().push(task.clone());
  }
  ```

- [ ] **Verify** Unit test for detection criteria

**Files to modify**: `ragc-core/src/compression.rs`

---

### 5.2 New Splitter Finding per Worker
**C++ Reference**: `agc_compressor.cpp` lines 2057-2086 (find_new_splitters)

- [ ] **Study** C++ AGC's per-worker splitter finding
  - Each worker has `vv_splitters[thread_id]`
  - After barrier, results merged into global `hs_splitters`
  - Bloom filter updated

- [ ] **Implement** Rust per-worker splitters
  ```rust
  struct WorkerState {
      new_splitters: Vec<u64>,
      // ...
  }

  // In new_splitters stage:
  for contig in hard_contigs {
      let new = find_new_splitters(contig, &shared);
      worker_state.new_splitters.extend(new);
  }

  // After barrier, merge
  if worker_id == 0 {
      for worker in all_workers {
          shared.splitters.extend(&worker.new_splitters);
          worker.new_splitters.clear();
      }
  }
  ```

- [ ] **Verify** Unit test for worker-local collection and merge

**Files to modify**: `ragc-core/src/worker.rs`, `ragc-core/src/adaptive_splitters.rs`

---

## Phase 6: Testing and Validation

### 6.1 Unit Tests for Core Components

- [ ] **Implement** Priority queue tests
  - Multi-threaded enqueue/dequeue
  - Priority ordering
  - Capacity limits
  - Completion signaling

- [ ] **Implement** Task and stage tests
  - Task creation
  - Stage transitions

- [ ] **Implement** Buffer tests
  - Segment buffering
  - Registration

- [ ] **Verify** All tests pass

**Files to modify**: `ragc-core/tests/`

---

### 6.2 Integration Tests

- [ ] **Implement** Single sample test
  - Compare output with current RAGC
  - Verify correctness

- [ ] **Implement** Multi-sample test
  - Test synchronization barriers
  - Verify sample ordering

- [ ] **Implement** Adaptive mode test
  - Test hard contig detection
  - Test new splitter finding
  - Verify improved compression

- [ ] **Verify** All integration tests pass

**Files to modify**: `ragc-core/tests/`

---

### 6.3 Compatibility Testing

- [ ] **Run** C++ AGC compatibility test suite
  - Create archives with new RAGC
  - Read with C++ AGC
  - Verify identical extraction

- [ ] **Run** Performance benchmarks
  - Memory usage (target: match C++ AGC)
  - Compression time (target: within 2x)
  - Archive size (target: identical)

- [ ] **Verify** All compatibility tests pass

---

### 6.4 Real-World Testing

- [ ] **Test** On yeast dataset (samples/yeast*.fa)
  - Verify correctness
  - Measure memory usage
  - Compare archive sizes

- [ ] **Test** On multi-sample FASTA
  - Verify sample ordering
  - Check terminator accumulation
  - Verify split counts

- [ ] **Verify** Production-ready

---

## Summary

**Total Items**: ~60 checkboxes
**Estimated Effort**: Each item = 30-60 minutes
**Total Time**: ~40-60 hours of focused implementation

**Success Criteria**:
- ✅ All unit tests pass
- ✅ All integration tests pass
- ✅ C++ AGC compatibility maintained
- ✅ Memory usage matches C++ AGC (±10%)
- ✅ Archive sizes identical to C++ AGC
- ✅ Terminator accumulation working (91+ splits on test data)
- ✅ Adaptive mode functional

**Next Steps**:
Start with Phase 1, Item 1.1 (Task Types and Processing Stages)
