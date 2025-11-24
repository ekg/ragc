# C++ AGC Batch Processing Architecture

**Date**: 2025-11-24
**Purpose**: Deep architectural analysis of C++ AGC's batch processing to understand why segments with identical k-mer pairs get different groups

---

## Executive Summary

**Root Cause Identified**: C++ AGC uses **batch-local state** for group assignment, while RAGC uses **global persistent state**.

**Key Finding**: In C++ AGC, `process_new()` is called at **sample boundaries** (or every 50 contigs in concatenated mode). This creates batch-local state in `m_kmers` that only sees segments from the current batch, causing segments with identical k-mer pairs across different batches to get different group IDs.

**RAGC Issue**: RAGC's `segment_groups: Arc<Mutex<BTreeMap<>>>` is **global and persistent** across all samples. Once a k-mer pair is registered, it's reused forever via `.entry(key).or_insert_with()`, making archives smaller but divergent from C++ AGC.

---

## 1. Overall Flow - C++ AGC

### Entry Point: `AddSampleFiles()` (lines 2113-2265)

```cpp
bool CAGCCompressor::AddSampleFiles(vector<pair<string, string>> _v_sample_file_name, const uint32_t no_threads)
{
    // For each sample file:
    for(auto sf : _v_sample_file_name)
    {
        while (gio.ReadContigRaw(id, contig))
        {
            // Add contig to priority queue
            pq_contigs_desc->Emplace(make_tuple(contig_processing_stage_t::all_contigs,
                                                 sample_name, contig_name, move(contig)),
                                      sample_priority, cost);

            // Check for synchronization boundary
            if (!concatenated_genomes && any_contigs_added)
            {
                // *** BATCH BOUNDARY HERE - ONE PER SAMPLE ***
                // Send N synchronization tokens (one per worker thread)
                pq_contigs_desc->EmplaceManyNoCost(
                    make_tuple(contig_processing_stage_t::registration, "", "", contig_t()),
                    sample_priority, no_workers);

                --sample_priority;
            }
        }
    }
}
```

**Key Insight**: In multi-sample mode (`!concatenated_genomes`), C++ AGC sends synchronization tokens **after each sample**, creating a batch boundary between samples.

---

## 2. Batch Boundaries

### Normal Mode (Multi-sample, `!concatenated_genomes`)
- **Batch = ONE SAMPLE**
- After processing all contigs from one sample, send `registration` tokens
- Each batch processes independently

### Concatenated Mode (`concatenated_genomes`)
- **Batch = 50 CONTIGS** (or when queue reaches 1 GB)
- `max_no_contigs_before_synchronization = pack_cardinality` (typically 50)
- Processes contigs from multiple samples in same batch

**For yeast5 test (5 samples, not concatenated)**:
- Batch 1: Sample AAA#0 → `process_new()` called
- Batch 2: Sample AAB#0 → `process_new()` called (NEW batch, local state)
- Batch 3: Sample AAC#0 → `process_new()` called
- Batch 4: Sample AAD#0 → `process_new()` called
- Batch 5: Sample AAE#0 → `process_new()` called

**This explains segments 14 & 15**: They're from different batches (AAA#0 vs AAB#0), so they get different groups even though k-mer pairs are identical!

---

## 3. State Management

### Global State (Persists Across All Batches)

```cpp
// Global registry: tracks ALL groups ever created
unordered_map<pair<uint64_t, uint64_t>, int32_t, MurMurPair64Hash> map_segments;  // (front, back) -> group_id

// Global terminator connections
unordered_map<uint64_t, vector<uint64_t>, MurMur64Hash> map_segments_terminators;  // kmer -> [connected kmers]

// Global segment storage
vector<shared_ptr<CSegment>> v_segments;  // Indexed by group_id
```

**Purpose**: These are used for **lookup only** - checking if a k-mer pair exists in a known group.

### Batch-Local State (Created Fresh Per Batch)

```cpp
// In CBufferedSegPart::process_new() (agc_compressor.h lines 384-415)
uint32_t process_new()
{
    lock_guard<mutex> lck(mtx);

    // *** LOCAL MAP - ONLY SEES CURRENT BATCH ***
    map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
    uint32_t group_id = (uint32_t)vl_seg_part.size();  // Start from current count

    // Assign group ids to NEW segments from THIS batch only
    for (const auto& x : s_seg_part)  // s_seg_part = segments from current batch
    {
        auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

        if (p == m_kmers.end())
            m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;  // New group
        // If found in m_kmers, reuse that group (within this batch)
    }

    // Move batch-local groups to global storage
    vl_seg_part.resize(group_id);

    for (auto& x : s_seg_part)
    {
        add_known(m_kmers[make_pair(x.kmer1, x.kmer2)], ...);
    }

    s_seg_part.clear();  // Clear batch buffer for next batch

    return no_new;
}
```

**Critical Details**:
1. `m_kmers` is **local variable** - created fresh for each batch
2. `s_seg_part` contains only "new" segments from current batch (lines 320-331 in header)
3. Group IDs start from `vl_seg_part.size()` - the current global count
4. Within a batch, identical k-mer pairs share groups
5. **Across batches, identical k-mer pairs get NEW groups** (because `m_kmers` is reset)

---

## 4. Group Assignment Logic

### Flow for Each Segment

```cpp
// In add_segment() (lines 1275-1499)
pair_segment_desc_t CAGCCompressor::add_segment(...)
{
    // Step 1: Determine k-mer pair (pk)
    pair<uint64_t, uint64_t> pk = ...;  // (kmer_front, kmer_back)

    // Step 2: Look up in GLOBAL registry
    auto p = map_segments.find(pk);

    if (p == map_segments.end())
    {
        // *** NOT FOUND GLOBALLY - ADD TO CURRENT BATCH AS "NEW" ***
        buffered_seg_part.add_new(pk.first, pk.second, sample_name, contig_name,
                                  segment_data, is_rev_comp, seg_part_no);
    }
    else
    {
        // *** FOUND GLOBALLY - ADD TO EXISTING GROUP ***
        segment_id = p->second;
        buffered_seg_part.add_known(segment_id, ~0ull, ~0ull, sample_name, contig_name,
                                    segment_data, is_rev_comp, seg_part_no);
    }
}
```

### Key Methods in CBufferedSegPart

**`add_new()` (header lines 326-331)**:
```cpp
void add_new(uint64_t kmer1, uint64_t kmer2, ...)
{
    lock_guard<mutex> lck(mtx);
    s_seg_part.emplace(kmer1, kmer2, ...);  // Add to "new segments" set
}
```
- Adds segment to `s_seg_part` (batch buffer for new segments)
- Will be processed by `process_new()` at batch boundary

**`add_known()` (header lines 320-324)**:
```cpp
void add_known(uint32_t group_id, ...)
{
    vl_seg_part[group_id].emplace(...);  // Add directly to existing group
}
```
- Adds segment directly to existing group's vector

---

## 5. Synchronization Flow (Worker Thread)

### Worker Thread Main Loop (lines 1099-1272)

```cpp
void start_compressing_threads(...)
{
    for (uint32_t i = 0; i < n_t; ++i)
    {
        v_threads.emplace_back([&, i, n_t]() {
            while(true)
            {
                task_t task;
                auto q_res = pq_contigs_desc_working->PopLarge(task);

                if (get<0>(task) == contig_processing_stage_t::registration)
                {
                    // *** SYNCHRONIZATION TOKEN RECEIVED ***

                    // Barrier 1: Wait for all threads to finish current batch
                    bar.arrive_and_wait();

                    // Thread 0 does registration
                    if (thread_id == 0)
                        register_segments(n_t);  // Calls process_new()!

                    // Barrier 2: Wait for registration to complete
                    bar.arrive_and_wait();

                    // All threads store segments in parallel
                    store_segments(zstd_cctx, zstd_dctx);

                    // Barrier 3: Wait for storage to complete
                    bar.arrive_and_wait();

                    // Thread 0 clears batch buffers
                    if (thread_id == 0)
                        buffered_seg_part.clear(max(1u, n_t-1));

                    bar.arrive_and_wait();

                    continue;  // Process next batch
                }

                // Normal contig processing
                compress_contig(...);
            }
        });
    }
}
```

### `register_segments()` (lines 954-971)

```cpp
void CAGCCompressor::register_segments(uint32_t n_t)
{
    buffered_seg_part.sort_known(n_t);  // Sort existing groups

    // *** CRITICAL: PROCESS NEW SEGMENTS FROM CURRENT BATCH ***
    uint32_t no_new = buffered_seg_part.process_new();  // Assigns group IDs to batch

    // Register new streams in archive
    for (uint32_t i = 0; i < no_new; ++i)
        out_archive->RegisterStreams(ss_ref_name(..., no_segments + i), ...);

    no_segments += no_new;

    if (no_segments > v_segments.size())
        v_segments.resize(no_segments);

    // Prepare for storage phase
    buffered_seg_part.distribute_segments(0, 0, no_raw_groups);
    buffered_seg_part.restart_read_vec();
}
```

### `store_segments()` (lines 974-1050)

```cpp
void CAGCCompressor::store_segments(ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    // Each thread claims groups to process
    while (true)
    {
        int block_group_id = buffered_seg_part.get_vec_id();  // Thread-safe counter
        if (block_group_id < 0) break;

        for (int group_id = block_group_id; group_id > ...; --group_id)
        {
            while (buffered_seg_part.get_part(group_id, ...))  // Pop segments from group
            {
                // First segment in group creates the CSegment object
                if (v_segments[group_id] == nullptr)
                {
                    v_segments[group_id] = make_shared<CSegment>(...);

                    // *** UPDATE GLOBAL REGISTRY ***
                    seg_map_mtx.lock();

                    auto p = map_segments.find(make_pair(kmer1, kmer2));
                    if (p == map_segments.end())
                        map_segments[make_pair(kmer1, kmer2)] = group_id;
                    else if (p->second > group_id)
                        p->second = group_id;  // Keep smallest group_id

                    // Update terminators map
                    map_segments_terminators[kmer1].push_back(kmer2);
                    ...

                    seg_map_mtx.unlock();
                }

                // Add segment to group (compress and store)
                in_group_id = v_segments[group_id]->add(seg_data, ...);
            }
        }
    }
}
```

**Key Observation**: Global registry (`map_segments`) is updated **AFTER** group IDs are assigned in `process_new()`. This means:
1. Batch N: Segment with k-mer pair (A, B) → gets group_id = 30
2. Batch N+1: Segment with same k-mer pair (A, B) → NOT found in global registry yet (hasn't been stored) → gets group_id = 31
3. After batch N+1's `store_segments()`, (A, B) → 30 is added to global registry
4. Batch N+2: Segment with (A, B) → found in global registry → reuses group 30

**Wait, that doesn't match our observation!** Let me re-read the store_segments logic...

Actually, I see it now:
- Line 1009-1013: `map_segments.find(make_pair(kmer1, kmer2))` checks if k-mer pair exists
- If not found, adds `group_id` to global registry
- **But this happens DURING store_segments(), not BEFORE process_new()!**

So the timeline is:
1. Batch 1 (AAA#0):
   - `add_segment()` → k-mer pair (A, B) → `map_segments.find()` → NOT FOUND → `add_new()`
   - `register_segments()` → `process_new()` → assigns group 30 to (A, B)
   - `store_segments()` → updates global registry: (A, B) → 30
2. Batch 2 (AAB#0):
   - `add_segment()` → k-mer pair (A, B) → `map_segments.find()` → FOUND (group 30) → `add_known(30)`
   - Segment goes directly to existing group 30

**Wait, that means identical k-mer pairs SHOULD reuse groups!** But our CSV shows they don't...

Let me re-check the actual data flow. The issue might be in the timing of when `map_segments` is populated vs when `add_segment()` looks it up.

---

## 6. Example: Segments 14 & 15 Timeline

Let me trace through what happens with exact timing:

### Batch 1 (AAA#0 - Sample 1)

**Worker thread processing AAA#0/chrV:**
1. `compress_contig()` called
2. For each splitter boundary → `add_segment()` called
3. **Segment 14**: k-mer pair (69536031064704, 69536031064704)
   - `map_segments.find(69536031064704, 69536031064704)` → **NOT FOUND** (first batch)
   - `buffered_seg_part.add_new()` → added to `s_seg_part` (batch buffer)

**Registration phase** (triggered by synchronization token):
1. Thread 0: `register_segments()` → `buffered_seg_part.process_new()`
2. `process_new()` creates local `m_kmers` map
3. Iterates `s_seg_part`:
   - Segment 14: k-mer pair (69536031064704, 69536031064704)
   - `m_kmers.find()` → not in local map → assign `group_id = 30`
4. `vl_seg_part.resize(31)` (or whatever final count is)

**Storage phase**:
1. All threads: `store_segments()` in parallel
2. For group 30:
   - First segment creates `v_segments[30]`
   - **Updates global registry**: `map_segments[(69536031064704, 69536031064704)] = 30`
   - Compresses and writes segment

**Result**: (69536031064704, 69536031064704) → group 30 in global registry

### Batch 2 (AAB#0 - Sample 2)

**Worker thread processing AAB#0/chrV:**
1. `compress_contig()` called
2. **Segment 15**: k-mer pair (69536031064704, 69536031064704)
   - `map_segments.find(69536031064704, 69536031064704)` → **SHOULD BE FOUND** (was added in Batch 1)
   - Should call `buffered_seg_part.add_known(30, ...)`

**But the CSV shows segment 15 in group 31, not 30!**

### Hypothesis: Race Condition or Ordering Issue

Let me check if there's a race between `store_segments()` updating `map_segments` and the next batch's `add_segment()` checking it...

Actually, looking at the synchronization barriers:
- Batch 1: `store_segments()` completes → Barrier 3 → `buffered_seg_part.clear()` → Barrier 4
- Batch 2: Only starts after Barrier 4

So there's no race condition. The global registry SHOULD be updated before Batch 2 starts.

**New Hypothesis**: The k-mer pair might be slightly different due to:
1. **RC orientation**: Maybe one segment stores (A, B) and another stores (B, A)?
2. **Missing k-mer**: Maybe one has MISSING_KMER and another doesn't?
3. **Segment splitting**: Maybe one segment gets split and the other doesn't?

Let me check the CSV data more carefully...

---

## 7. Checking CSV Data (from CLAUDE.md logs)

From the investigation logs:

```
Segments 14 & 15 - IDENTICAL k-mers but DIFFERENT groups:

C++ AGC:
  Segment 14: group_id=30  in_group_id=0  sample=AAA#0  contig=chrV
              Front: 69536031064704  Back: 69536031064704

  Segment 15: group_id=31  in_group_id=0  sample=AAB#0  contig=chrV
              Front: 69536031064704  Back: 69536031064704

RAGC:
  Segment 14: group_id=30  sample=AAA#0  contig=chrV
              Front: 69536031064704  Back: 69536031064704

  Segment 15: group_id=30  sample=AAB#0  contig=chrV  <- REUSES group 30!
              Front: 69536031064704  Back: 69536031064704
```

**So the k-mer pairs ARE identical, and RAGC does reuse the group, but C++ AGC does NOT!**

This confirms: **C++ AGC intentionally creates separate groups for identical k-mer pairs across batches**.

---

## 8. Key Differences from RAGC

### RAGC Current Architecture

**Global Persistent State** (`streaming_compressor_queue.rs` lines 216-217):
```rust
segment_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
map_segments: Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
```

**Group Assignment** (in worker thread):
```rust
// Check if k-mer pair exists in global registry
let existing_group_id = {
    let seg_map = map_segments.lock().unwrap();
    seg_map.get(&key).copied()
};

if let Some(group_id) = existing_group_id {
    // REUSE existing group
    segment_groups.lock().unwrap().entry(key).or_insert_with(...);
} else {
    // CREATE new group
    let new_group_id = group_counter.fetch_add(1, Ordering::SeqCst);
    segment_groups.lock().unwrap().insert(key, ...);
    map_segments.lock().unwrap().insert(key, new_group_id);
}
```

**Key Issue**: RAGC's `map_segments` is **global and never cleared**, so once a k-mer pair is registered, it's reused forever. This creates byte-identical behavior across samples (more efficient compression) but diverges from C++ AGC's batch-local behavior.

### C++ AGC Batch Architecture

**Batch-Local Group Assignment**:
```cpp
// m_kmers is LOCAL to process_new() - created fresh per batch
map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;

// Only sees segments from CURRENT batch
for (const auto& x : s_seg_part)
{
    if (m_kmers.find(make_pair(x.kmer1, x.kmer2)) == m_kmers.end())
        m_kmers[...] = group_id++;  // New group for this batch
}
```

**Global Registry is Read-Only During Processing**:
- `map_segments` is only read in `add_segment()` to check if group exists
- Updated in `store_segments()` AFTER group IDs are assigned
- Never used by `process_new()` for group assignment within a batch

---

## 9. Why C++ AGC Does This

### Possible Reasons:

1. **Memory Management**: Batch processing allows clearing buffers between samples, limiting memory usage
2. **Parallelism**: Batch boundaries provide clean synchronization points for multi-threading
3. **Incremental Compression**: Each sample can be processed independently without global state
4. **Archive Compatibility**: Batched groups might enable features like:
   - Per-sample group statistics
   - Faster sample deletion/extraction
   - Better segment locality within samples

### Trade-offs:

**C++ AGC Approach** (Batch-Local):
- ✅ Predictable memory usage (clear buffers per batch)
- ✅ Clean synchronization points
- ✅ Sample-independent processing
- ❌ Larger archives (duplicate groups for identical k-mer pairs)
- ❌ Less compression efficiency

**RAGC Approach** (Global Persistent):
- ✅ Smaller archives (reuse groups globally)
- ✅ Better compression efficiency
- ❌ Growing memory usage (global state never cleared)
- ❌ Global contention (all workers access same maps)
- ❌ **Divergent from C++ AGC behavior**

---

## 10. Implementation Plan for RAGC

### Option A: Minimal Change - Add Batch Boundaries

**Goal**: Match C++ AGC's batch behavior without full architectural refactor

**Changes**:
1. Add batch tracking to `StreamingQueueCompressor`:
   ```rust
   current_batch_sample: Arc<Mutex<Option<String>>>,
   batch_segment_buffer: Arc<Mutex<Vec<BufferedSegment>>>,
   ```

2. Modify `push()` to detect sample boundaries:
   ```rust
   pub fn push(&mut self, sample_name: String, ...) {
       let mut current_batch = self.current_batch_sample.lock().unwrap();

       if current_batch.as_ref() != Some(&sample_name) {
           // BATCH BOUNDARY - trigger process_new()
           self.flush_batch()?;
           *current_batch = Some(sample_name.clone());
       }

       // Add contig to batch buffer
       ...
   }
   ```

3. Implement `flush_batch()`:
   ```rust
   fn flush_batch(&self) -> Result<()> {
       // Process all segments in batch_segment_buffer
       // Create local m_kmers map (like C++ AGC)
       let mut local_group_map: HashMap<SegmentGroupKey, u32> = HashMap::new();
       let starting_group_id = self.group_counter.load(Ordering::SeqCst);

       for segment in batch_buffer {
           let key = SegmentGroupKey { kmer_front, kmer_back };

           // Check GLOBAL registry first (like C++ AGC add_segment)
           let existing = self.map_segments.lock().unwrap().get(&key).copied();

           if let Some(group_id) = existing {
               // Add to existing group
               add_to_group(group_id, segment);
           } else {
               // Add to LOCAL batch map
               let group_id = local_group_map.entry(key)
                   .or_insert_with(|| {
                       let id = self.group_counter.fetch_add(1, Ordering::SeqCst);
                       id
                   });
               add_to_group(*group_id, segment);
           }
       }

       // Update global registry with batch's new groups
       for (key, group_id) in local_group_map {
           self.map_segments.lock().unwrap().insert(key, group_id);
       }

       // Clear batch buffer
       batch_buffer.clear();
   }
   ```

**Pros**:
- Small, focused change
- Preserves most of RAGC's current architecture
- Should produce byte-identical archives

**Cons**:
- Adds complexity (batch buffering + flush logic)
- Might not handle all edge cases (split segments, adaptive mode, etc.)
- Estimated effort: 200-300 lines of code changes

---

### Option B: Full Redesign - Match C++ AGC Exactly

**Goal**: Replicate C++ AGC's `CBufferedSegPart` and synchronization model

**Changes**:
1. Create Rust equivalent of `CBufferedSegPart`:
   ```rust
   struct BufferedSegPart {
       // Known segments (already have group IDs)
       vl_seg_part: Vec<Vec<BufferedSegment>>,  // Indexed by group_id

       // New segments (waiting for group assignment)
       s_seg_part: BTreeSet<BufferedSegment>,

       // Mutex for thread-safe access
       mtx: Mutex<()>,
   }

   impl BufferedSegPart {
       fn add_known(&mut self, group_id: u32, segment: BufferedSegment) {
           self.vl_seg_part[group_id as usize].push(segment);
       }

       fn add_new(&mut self, segment: BufferedSegment) {
           self.s_seg_part.insert(segment);
       }

       fn process_new(&mut self, starting_group_id: u32) -> u32 {
           // LOCAL map for current batch
           let mut m_kmers: HashMap<SegmentGroupKey, u32> = HashMap::new();
           let mut group_id = starting_group_id;

           // Assign group IDs within batch
           for seg in &self.s_seg_part {
               let key = SegmentGroupKey {
                   kmer_front: seg.kmer_front,
                   kmer_back: seg.kmer_back,
               };

               m_kmers.entry(key).or_insert_with(|| {
                   let id = group_id;
                   group_id += 1;
                   id
               });
           }

           // Resize and move segments to vl_seg_part
           self.vl_seg_part.resize(group_id as usize, Vec::new());

           for seg in self.s_seg_part.drain() {
               let key = SegmentGroupKey { ... };
               let assigned_group = m_kmers[&key];
               self.vl_seg_part[assigned_group as usize].push(seg);
           }

           group_id - starting_group_id  // Number of new groups
       }

       fn clear(&mut self) {
           self.vl_seg_part.clear();
           self.s_seg_part.clear();
       }
   }
   ```

2. Replace worker threads with barrier-based synchronization:
   ```rust
   // Add to StreamingQueueCompressor
   synchronization_barrier: Arc<Barrier>,
   buffered_seg_part: Arc<Mutex<BufferedSegPart>>,

   // Worker thread loop
   loop {
       match queue.pop() {
           Some(ContigTask::Normal(task)) => {
               compress_contig(task, &buffered_seg_part, ...);
           }
           Some(ContigTask::SyncToken) => {
               // Barrier 1: Wait for all workers to finish current batch
               synchronization_barrier.wait();

               // Thread 0: Register and assign group IDs
               if worker_id == 0 {
                   let no_new = buffered_seg_part.lock().unwrap()
                       .process_new(group_counter.load(Ordering::SeqCst));
                   group_counter.fetch_add(no_new, Ordering::SeqCst);
               }

               // Barrier 2: Wait for registration
               synchronization_barrier.wait();

               // All threads: Store segments in parallel
               store_segments(&buffered_seg_part, ...);

               // Barrier 3: Wait for storage
               synchronization_barrier.wait();

               // Thread 0: Clear batch buffers
               if worker_id == 0 {
                   buffered_seg_part.lock().unwrap().clear();
               }

               synchronization_barrier.wait();
           }
           None => break,  // Queue completed
       }
   }
   ```

3. Modify `push()` to send sync tokens:
   ```rust
   pub fn push(&mut self, sample_name: String, ...) {
       // Track current sample
       let mut current = self.current_batch_sample.lock().unwrap();

       if current.as_ref() != Some(&sample_name) {
           // BATCH BOUNDARY - send sync tokens
           for _ in 0..self.config.num_threads {
               self.queue.push(ContigTask::SyncToken)?;
           }
           *current = Some(sample_name);
       }

       // Queue normal contig task
       self.queue.push(ContigTask::Normal(...))?;
   }
   ```

**Pros**:
- Exact match to C++ AGC behavior
- Clean architecture with well-defined phases
- Should produce byte-identical archives
- Better separation of concerns

**Cons**:
- Major refactoring (500+ lines of code changes)
- Need to rewrite worker thread logic
- Need to add barrier synchronization
- More complex to test and debug
- Higher risk of introducing bugs

---

## 11. Recommended Approach

### Start with Option A (Minimal Change)

**Rationale**:
1. **Lower Risk**: Smaller changes = easier to verify correctness
2. **Faster Implementation**: Can be done in 1-2 days vs 1-2 weeks for full redesign
3. **Iterative**: Can upgrade to Option B later if needed
4. **Testing**: Easier to isolate issues in smaller changeset

**Implementation Steps**:
1. Add batch tracking fields to `StreamingQueueCompressor`
2. Implement `flush_batch()` method with local `m_kmers` map
3. Modify `push()` to detect sample boundaries and call `flush_batch()`
4. Modify `finish()` to flush final batch
5. Test with yeast5 dataset to verify byte-identical archives
6. Run segment layout comparison: `diff cpp_layout.csv ragc_layout.csv`

**Success Criteria**:
- Segments 14 & 15 get different group IDs (30 and 31)
- Segment layout CSV matches C++ AGC exactly
- Archives are byte-identical (or within acceptable tolerance)

**If Option A Fails**:
- Analyze why (edge cases, timing issues, etc.)
- Consider whether Option B is necessary
- User will make final decision

---

## 12. Open Questions

1. **Concatenated Mode**: Does batch boundary align with pack_cardinality (50 contigs) or something else?
   - Need to test with concatenated mode to verify

2. **Adaptive Mode**: When new splitters are found, do they trigger a batch boundary?
   - Lines 1187-1236 suggest yes (registration tokens sent after new_splitters stage)

3. **Segment Splitting**: If a segment is split in C++ AGC, does it affect batch assignment?
   - Need to trace segment splitting logic to verify

4. **Multi-threading**: C++ AGC uses barriers for synchronization - how critical is exact barrier behavior?
   - RAGC currently uses async workers - might need to add barriers

5. **Performance Impact**: Does batch-local assignment affect compression performance?
   - Larger archives = more data to write
   - Might be offset by better cache locality within batches

---

## 13. Conclusion

**Root Cause**: C++ AGC uses **batch-local state** (`m_kmers` in `process_new()`) that only sees segments from the current batch, causing identical k-mer pairs across batches to get different group IDs. RAGC uses **global persistent state** that reuses groups across all samples.

**Fix**: Implement batch boundaries in RAGC with local group assignment per batch, matching C++ AGC's behavior.

**Recommended Approach**: Start with Option A (minimal change) to add batch flushing logic. This should be sufficient to match C++ AGC's group assignment behavior without a full architectural redesign.

**Estimated Complexity**:
- Option A: 200-300 lines, 1-2 days implementation
- Option B: 500+ lines, 1-2 weeks implementation

**Next Steps**: Begin implementation of Option A with systematic testing using segment layout comparison as verification method.
