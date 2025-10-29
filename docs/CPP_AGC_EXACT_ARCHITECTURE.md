# C++ AGC Exact Architecture: Complete Flow Analysis

**Date**: 2025-10-29
**Purpose**: Document C++ AGC's exact processing flow to understand why it achieves 91 split attempts vs RAGC's 37
**Goal**: Rewrite RAGC to match C++ AGC's behavior exactly

---

## Executive Summary

**Root Cause of Split Difference**: C++ AGC uses a **batched processing model with synchronization points** where terminators from batch N become available for split checking in batch N+1 and later. RAGC's inline approach processes everything in one pass, missing opportunities.

**Key Finding**: C++ AGC's architecture creates a "terminator accumulation effect" where:
- Sample AAA's segments create terminators during its `registration` phase
- Sample AAB can then use AAA's terminators for split checking
- Sample AAC can use terminators from both AAA and AAB
- This explains why all 91 split attempts come from AAB and AAC (none from AAA)

---

## C++ AGC Complete Processing Flow

### Phase 1: Initialization (Lines 350-380)

```cpp
// If loading from existing archive:
for (each stored segment group) {
    map_segments[make_pair(x1, x2)] = group_id;

    // ADD TERMINATORS FROM EXISTING ARCHIVE
    if (x1 != ~0ull && x2 != ~0ull) {
        map_segments_terminators[x1].push_back(x2);
        if (x1 != x2)
            map_segments_terminators[x2].push_back(x1);
    }
}
```

**Key**: If appending to an archive, terminators are loaded upfront. For new archives, `map_segments_terminators` starts empty.

### Phase 2: Main Processing Loop (Lines 2163-2242)

```cpp
for (auto sf : _v_sample_file_name) {  // For each sample file
    while (gio.ReadContigRaw(id, contig)) {  // For each contig in sample
        // Register contig in collection
        collection_desc->register_sample_contig(sf.first, id);

        // Enqueue contig to priority queue with sample_priority
        pq_contigs_desc->Emplace(
            make_tuple(contig_processing_stage_t::all_contigs,
                      sf.first, id, move(contig)),
            sample_priority, cost);
    }

    // After ALL contigs from this sample are enqueued:
    // Send synchronization tokens (one per worker thread)
    pq_contigs_desc->EmplaceManyNoCost(
        make_tuple(adaptive_compression ?
                  contig_processing_stage_t::new_splitters :
                  contig_processing_stage_t::registration,
                  "", "", contig_t()),
        sample_priority, no_workers);

    --sample_priority;  // Lower priority for next sample
}
```

**Critical Behavior**:
- Contigs enqueued with `contig_processing_stage_t::all_contigs`
- After each **sample** (not contig), synchronization tokens sent
- Priority decreases: sample 1 gets highest priority, sample 2 next, etc.
- Workers process higher-priority items first

### Phase 3: Worker Thread Processing (Lines 1097-1270)

Each worker thread runs this loop:

```cpp
v_threads.emplace_back([&, i, n_t]() {
    auto zstd_cctx = ZSTD_createCCtx();
    auto zstd_dctx = ZSTD_createDCtx();

    while(true) {
        task_t task;
        auto q_res = pq_contigs_desc_working->PopLarge(task);

        if (q_res == completed) break;
        if (q_res == empty) continue;

        // Handle synchronization token
        if (get<0>(task) == contig_processing_stage_t::registration) {
            bar.arrive_and_wait();  // Synchronize all workers

            if (thread_id == 0)
                register_segments(n_t);  // Register new groups

            bar.arrive_and_wait();

            store_segments(zstd_cctx, zstd_dctx);  // CRITICAL: Adds terminators!

            bar.arrive_and_wait();
            // ... cleanup and flushing
            bar.arrive_and_wait();

            continue;  // Go back to processing contigs
        }

        // Handle new splitters phase
        if (get<0>(task) == contig_processing_stage_t::new_splitters) {
            // Add discovered splitters to bloom filter
            // Re-queue hard contigs that couldn't be compressed
            continue;
        }

        // Process regular contig
        if (get<0>(task) == contig_processing_stage_t::all_contigs) {
            preprocess_raw_contig(get<3>(task));
        }

        compress_contig(get<0>(task), get<1>(task), get<2>(task),
                       get<3>(task), zstd_cctx, zstd_dctx, thread_id, bar);
    }
});
```

**Critical Flow**:
1. Workers pop contigs from priority queue
2. Process contigs from sample 1 (highest priority) first
3. When synchronization token hit → all workers synchronize
4. Thread 0 calls `register_segments()` and `store_segments()`
5. `store_segments()` **adds terminators to `map_segments_terminators`**
6. Workers continue, now processing sample 2 with sample 1's terminators available

### Phase 4: Contig Compression (Lines 2020-2053)

```cpp
bool CAGCCompressor::compress_contig(...) {
    // Split contig at splitter k-mers
    for (each k-mer in contig) {
        if (bloom_splitters.check(d) && hs_splitters.check(d)) {
            // Found a splitter → call add_segment for this piece
            auto seg_id = add_segment(sample_name, id, seg_part_no,
                move(get_part(contig, split_pos, pos + 1 - split_pos)),
                split_kmer, kmer, zstd_cctx, zstd_dctx, thread_id, bar);

            ++seg_part_no;
            split_kmer = kmer;
            split_pos = pos + 1;
        }
    }

    // Add final segment
    add_segment(sample_name, id, seg_part_no, ...);
}
```

### Phase 5: Segment Addition with Split Checking (Lines 1275-1504)

**This is where split checking happens:**

```cpp
pair_segment_desc_t CAGCCompressor::add_segment(
    const string& sample_name, const string& contig_name, uint32_t seg_part_no,
    contig_t &&segment, CKmer kmer_front, CKmer kmer_back, ...) {

    // Normalize k-mer pair to canonical form
    pk = make_pair(kmer_front.canonical(), kmer_back.canonical());

    // Check if this segment group already exists
    auto p = map_segments.find(pk);

    // SPLIT CHECKING: If segment doesn't exist AND both endpoint k-mers
    // are terminators, try to split
    if (!concatenated_genomes &&
        p == map_segments.end() &&
        pk.first != ~0ull && pk.second != ~0ull &&
        map_segments_terminators.count(pk.first) &&     // ← CRITICAL CHECK
        map_segments_terminators.count(pk.second)) {    // ← CRITICAL CHECK

        // Attempt to find middle splitter k-mer that connects both endpoints
        auto split_match = find_cand_segment_with_missing_middle_splitter(
            kmer1, kmer2, segment, segment_rc, zstd_dctx, bar);

        if (split_match.first != ~0ull) {
            // SPLIT SUCCESSFUL
            // Split segment into two parts:
            // 1. (kmer_front → middle_kmer)
            // 2. (middle_kmer → kmer_back)

            segment_id = map_segments.at(make_pair(kmer_front, middle_kmer));
            segment_id2 = map_segments.at(make_pair(middle_kmer, kmer_back));

            // Both groups MUST exist (that's why split succeeded)
        }
    }

    // If segment group doesn't exist, buffer as NEW
    if (p == map_segments.end()) {
        buffered_seg_part.add_new(pk.first, pk.second, sample_name,
                                  contig_name, segment, ...);
    }
    // If segment group exists (or split succeeded), buffer as KNOWN
    else {
        buffered_seg_part.add_known(segment_id, sample_name,
                                    contig_name, segment, ...);

        if (segment_id2 >= 0)  // If split into 2 segments
            buffered_seg_part.add_known(segment_id2, ...);
    }
}
```

**Critical Logic**:
- Split checking requires BOTH `pk.first` and `pk.second` to be in `map_segments_terminators`
- This check happens DURING segment processing (inline)
- Terminators are only present if they were added by previous batches

### Phase 6: Store Segments (Lines 995-1049)

**This is where terminators are ADDED:**

```cpp
void CAGCCompressor::store_segments(...) {
    for (each group_id with buffered segments) {
        // Get segments from buffer
        while (buffered_seg_part.get_part(group_id, kmer1, kmer2, ...)) {
            if (v_segments[group_id] == nullptr) {
                // Create new segment group
                v_segments[group_id] = make_shared<CSegment>(...);

                seg_map_mtx.lock();

                // Register group in map_segments
                auto p = map_segments.find(make_pair(kmer1, kmer2));
                if (p == map_segments.end())
                    map_segments[make_pair(kmer1, kmer2)] = group_id;

                // ADD TERMINATORS ← THIS IS THE KEY!
                if (kmer1 != ~0ull && kmer2 != ~0ull) {
                    map_segments_terminators[kmer1].push_back(kmer2);
                    sort(map_segments_terminators[kmer1].begin(),
                         map_segments_terminators[kmer1].end());

                    if (kmer1 != kmer2) {
                        map_segments_terminators[kmer2].push_back(kmer1);
                        sort(map_segments_terminators[kmer2].begin(),
                             map_segments_terminators[kmer2].end());
                    }
                }

                seg_map_mtx.unlock();
            }

            // Add segment data to group
            in_group_id = v_segments[group_id]->add(seg_data, ...);
        }
    }
}
```

**Critical Behavior**:
- Called during `registration` synchronization phase
- Processes buffered segments from the PREVIOUS batch
- Creates new groups and **adds their k-mers to `map_segments_terminators`**
- These terminators are then available for split checking in NEXT batch

---

## Timeline Example: 3 Samples (AAA, AAB, AAC)

### Batch 1: Sample AAA Processing

1. **Enqueue Phase**: All AAA contigs enqueued to priority queue
2. **Process Phase**: Workers pop AAA contigs, call `compress_contig()` → `add_segment()`
3. **Split Check**: `map_segments_terminators` is EMPTY → no splits possible
4. **Buffer**: All AAA segments buffered as NEW
5. **Sync Token**: Registration token hit, all workers synchronize
6. **Store Phase**: `store_segments()` creates groups from AAA segments
7. **Add Terminators**: AAA's k-mer pairs added to `map_segments_terminators`

**Result**: AAA segments processed, 0 splits, terminators now available

### Batch 2: Sample AAB Processing

1. **Enqueue Phase**: All AAB contigs enqueued
2. **Process Phase**: Workers pop AAB contigs, call `add_segment()`
3. **Split Check**: `map_segments_terminators` contains AAA's terminators!
   - If AAB segment has endpoints (k1, k2) where:
     - k1 exists in AAA terminators
     - k2 exists in AAA terminators
     - Middle k-mer k_mid connects them
     - Groups (k1→k_mid) and (k_mid→k2) both exist
   - Then split succeeds!
4. **Buffer**: AAB segments (some split, some not) buffered
5. **Sync Token**: Registration token hit
6. **Store Phase**: `store_segments()` creates groups from AAB segments
7. **Add Terminators**: Both AAA and AAB terminators now available

**Result**: AAB benefits from AAA's terminators, many splits!

### Batch 3: Sample AAC Processing

1. **Process Phase**: AAC contigs processed
2. **Split Check**: `map_segments_terminators` contains terminators from BOTH AAA and AAB
3. **Even more split opportunities** than AAB had!

**Result**: AAC benefits from both AAA and AAB terminators

---

## Why C++ AGC Gets 91 Attempts, RAGC Gets 37

### C++ AGC's Advantage

```
Batch 1 (AAA):
  - Process contigs
  - Terminators available: []
  - Split attempts: 0
  - Store segments → add terminators

Batch 2 (AAB):
  - Process contigs
  - Terminators available: [AAA's k-mers]
  - Split attempts: ~40-50
  - Store segments → add more terminators

Batch 3 (AAC):
  - Process contigs
  - Terminators available: [AAA's + AAB's k-mers]
  - Split attempts: ~40-50

Total: ~80-100 attempts
Actual: 91 attempts
```

### RAGC's Inline Approach (Current)

```
Process all 3 samples sequentially:
  For each contig:
    - Check splits (terminators available from PREVIOUS contigs only)
    - Add segment inline
    - Add terminators inline

Problem:
  - AAA's early contigs have NO terminators available
  - AAB's contigs only see terminators from AAA contigs processed BEFORE them
  - AAC's contigs see some AAA+AAB terminators, but not all

Total: 37 attempts (missing 54 opportunities)
```

**The Difference**: C++ AGC ensures ALL of sample N's terminators are available before processing sample N+1. RAGC processes contigs one-by-one, so terminators trickle in gradually.

---

## Detailed Function Call Stack

### When Processing a Contig

```
compress_genomes_mt()                       // Line 2163: main loop
  ↓
[Priority Queue] ← all_contigs tasks
  ↓
Worker Thread Loop                          // Line 1104
  ↓
pq_contigs_desc_working->PopLarge(task)     // Line 1108
  ↓
preprocess_raw_contig()                     // Line 1241 (optional)
  ↓
compress_contig()                           // Line 1246
  ↓
add_segment()                               // Line 2023, 2050
  ↓
[SPLIT CHECK at line 1366-1464]
  ↓
buffered_seg_part.add_new() or add_known()  // Line 1493 or 1500
```

### When Synchronization Token Hit

```
Worker Thread Loop
  ↓
PopLarge() returns registration token       // Line 1114
  ↓
bar.arrive_and_wait()                       // All workers sync
  ↓
register_segments()                         // Thread 0 only, line 1119
  ↓
bar.arrive_and_wait()
  ↓
store_segments()                            // ALL workers, line 1132
  ↓
  For each buffered segment:
    - Create CSegment group
    - Add to map_segments
    - ADD TO map_segments_terminators       // Lines 1017-1024
  ↓
bar.arrive_and_wait()
  ↓
Continue processing next batch
```

---

## Key Data Structures

### `map_segments_terminators`

```cpp
map<uint64_t, vector<uint64_t>> map_segments_terminators;
```

**Purpose**: For each k-mer, stores list of k-mers it's connected to in existing groups

**Example**:
```
If groups exist:
  (k1 → k2)
  (k1 → k3)
  (k3 → k4)

Then:
  map_segments_terminators[k1] = [k2, k3]
  map_segments_terminators[k2] = [k1]
  map_segments_terminators[k3] = [k1, k4]
  map_segments_terminators[k4] = [k3]
```

**Split Checking**: Segment with endpoints (ka, kb) can split if:
1. `map_segments_terminators.count(ka) > 0` (ka is a terminator)
2. `map_segments_terminators.count(kb) > 0` (kb is a terminator)
3. There exists k_mid in both `map_segments_terminators[ka]` and `map_segments_terminators[kb]`
4. Groups (ka→k_mid) and (k_mid→kb) both exist

### `buffered_seg_part`

```cpp
CBufferedSegmentsPart buffered_seg_part;
```

**Purpose**: Buffers segments during processing phase, flushed during `store_segments()`

**Behavior**:
- `add_new()`: Buffer segment with unknown group (will create new group)
- `add_known()`: Buffer segment with known group_id
- `get_part()`: Retrieve buffered segments during store phase

---

## Critical Differences vs RAGC

| Aspect | C++ AGC | RAGC (Current Inline) |
|--------|---------|----------------------|
| **Processing** | Batched by sample | Sequential by contig |
| **Synchronization** | Barriers between samples | None |
| **Terminators Added** | After each sample batch | Inline with each segment |
| **Terminator Availability** | All of sample N before sample N+1 | Trickles in contig-by-contig |
| **Split Opportunities** | High (91 attempts) | Low (37 attempts) |
| **Buffering** | Segments buffered until sync point | No buffering (inline write) |

---

## Next Steps: Rewrite Plan

See `RAGC_REWRITE_PLAN.md` for detailed implementation strategy to match C++ AGC's architecture.
