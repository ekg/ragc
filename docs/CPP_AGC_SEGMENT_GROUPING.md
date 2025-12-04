# C++ AGC Segment Grouping: Complete Technical Specification

**Date**: 2025-12-04
**Purpose**: Exact documentation of C++ AGC's segment grouping mechanism for comparison with RAGC

---

## Executive Summary

C++ AGC uses a **streaming priority queue architecture** where:
1. Segments are created during contig processing and buffered in `CBufferedSegPart`
2. NEW segments (unknown k-mer pairs) are grouped in batches at synchronization points
3. KNOWN segments are immediately assigned to existing groups
4. The grouping key is `(kmer1, kmer2)` where both can be `~0ull` (MISSING)
5. Groups are created with monotonically increasing IDs starting from `vl_seg_part.size()`
6. Processing order is determined by priority queue with priority based on sample insertion order

---

## 1. WHEN ARE SEGMENTS CREATED?

**Location**: `agc_compressor.cpp`, lines 2114-2168, function `compress_contig()`

### Trigger Points

Segments are created at **splitter boundaries**:

```cpp
// Line 2124-2153
for (auto x : contig)
{
    if (x >> 2)         // x > 3 (non-ACGT)
        kmer.Reset();
    else
    {
        kmer.insert_canonical(x);

        if (kmer.is_full())
        {
            uint64_t d = kmer.data_canonical();

            // Check if this k-mer is a splitter
            if (bloom_splitters.check(d) && hs_splitters.check(d))
            {
                // CREATE SEGMENT: from split_pos to pos+1
                auto seg_id = add_segment(sample_name, id, seg_part_no,
                    move(get_part(contig, split_pos, pos + 1 - split_pos)),
                    split_kmer, kmer, zstd_cctx, zstd_dctx, thread_id, bar);

                ++seg_part_no;
                if (seg_id.contains_second)
                    ++seg_part_no;

                split_pos = pos + 1 - kmer_length;
                split_kmer = kmer;
                kmer.Reset();
            }
        }
    }
    ++pos;
}
```

### Final Segment

After the loop, if any data remains:

```cpp
// Line 2163-2165
if (split_pos < contig.size())
    add_segment(sample_name, id, seg_part_no,
        move(get_part(contig, split_pos, contig.size() - split_pos)),
        split_kmer, CKmer(kmer_length, kmer_mode_t::canonical), ...);
```

**Key Point**: The final segment has `kmer_back = CKmer()` (empty/MISSING).

---

## 2. K-MER EXTRACTION

**Location**: `agc_compressor.cpp`, lines 1298-1567, function `add_segment()`

### Rules

K-mers are **extracted during segment creation** in `compress_contig()`:
- **Front k-mer**: `split_kmer` from previous splitter (or empty for first segment)
- **Back k-mer**: Current splitter k-mer (or empty for last segment)

### States

K-mers can be in 3 states:
1. **FULL**: `is_full() == true`, has valid canonical k-mer data
2. **EMPTY/MISSING**: `is_full() == false`, represents `~0ull` (all bits set)
3. **Data**: `kmer.data()` returns the 64-bit canonical k-mer value

### Canonical Orientation

**CRITICAL**: K-mers are stored in **canonical form** during extraction:

```cpp
// Line 2130
kmer.insert_canonical(x);

// Line 2134
uint64_t d = kmer.data_canonical();
```

This means **k-mer orientation is determined during extraction**, not during grouping.

---

## 3. GROUPING KEY CONSTRUCTION

**Location**: `agc_compressor.cpp`, lines 1298-1567, function `add_segment()`

### Case Analysis

The grouping key `pk = (kmer1, kmer2)` is constructed based on which k-mers are present:

#### Case 1: Both k-mers MISSING
```cpp
// Lines 1309-1323
if (!kmer_front.is_full() && !kmer_back.is_full())
{
    // Try fallback minimizers (if enabled)
    if (fallback_filter)
        tie(pk, store_rc) = find_cand_segment_using_fallback_minimizers(segment, 1);
    else
        pk = pk_empty;  // (~0ull, ~0ull)
}
```

**Key**: `(~0ull, ~0ull)`

#### Case 2: Both k-mers FULL
```cpp
// Lines 1325-1336
else if (kmer_front.is_full() && kmer_back.is_full())
{
    if (kmer_front.data() < kmer_back.data())
        pk = make_pair(kmer_front.data(), kmer_back.data());
    else
    {
        pk = make_pair(kmer_back.data(), kmer_front.data());
        reverse_complement_copy(segment, segment_rc);
        store_rc = true;
    }
}
```

**Key**: `(min(kmer_front, kmer_back), max(kmer_front, kmer_back))`
**Orientation**: Segment is stored in RC if `kmer_front > kmer_back`

#### Case 3a: Only front k-mer FULL
```cpp
// Lines 1338-1365
else if (kmer_front.is_full())
{
    CKmer kmer = kmer_front;
    reverse_complement_copy(segment, segment_rc);

    tie(pk, store_rc) = find_cand_segment_with_one_splitter(
        kmer, segment, segment_rc, zstd_dctx, bar);

    // If no match found, try fallback minimizers
    if (pk.first == ~0ull || pk.second == ~0ull)
    {
        tie(pk_alt, store_rc_alt) =
            find_cand_segment_using_fallback_minimizers(segment, 5);
        if (pk_alt != pk_empty)
        {
            pk = pk_alt;
            store_rc = store_rc_alt;
        }
    }
}
```

**Key**: Result from `find_cand_segment_with_one_splitter()` which tries to match against existing segments

#### Case 3b: Only back k-mer FULL
```cpp
// Lines 1367-1401
else if (kmer_back.is_full())
{
    CKmer kmer = kmer_back;
    kmer.swap_dir_rc();  // Swap orientation
    reverse_complement_copy(segment, segment_rc);

    tie(pk, store_dir) = find_cand_segment_with_one_splitter(
        kmer, segment_rc, segment, zstd_dctx, bar);
    store_rc = !store_dir;

    // Fallback minimizers if no match
    if (pk.first == ~0ull || pk.second == ~0ull)
    {
        tie(pk_alt, store_dir_alt) =
            find_cand_segment_using_fallback_minimizers(segment_rc, 5);
        if (pk_alt != pk_empty)
        {
            pk = pk_alt;
            store_rc = !store_dir_alt;
        }
    }
}
```

**Key**: Result from `find_cand_segment_with_one_splitter()` on RC segment

### Segment Splitting Logic

**Lines 1419-1528**: If initial lookup fails AND both splitters exist AND both splitters are known terminators, C++ AGC attempts to **split the segment** into two parts by finding a middle splitter.

---

## 4. GROUP LOOKUP LOGIC

**Location**: `agc_compressor.cpp`, line 1416

### Lookup

```cpp
// Line 1416
auto p = map_segments.find(pk);
```

### Data Structure

```cpp
// agc_compressor.h, line 628
unordered_map<pair<uint64_t, uint64_t>, int32_t, MurMurPair64Hash> map_segments;
```

**Type**: `unordered_map<(kmer1, kmer2) -> group_id>`

### Synchronization

```cpp
// agc_compressor.h, line 608
shared_mutex seg_map_mtx;
```

Reads use shared lock, writes use exclusive lock.

---

## 5. GROUP CREATION LOGIC

**Location**: `agc_compressor.h`, lines 384-415, function `CBufferedSegPart::process_new()`

### When Groups Are Created

Groups are created during **batch registration** at synchronization points:

```cpp
// Lines 388-397
map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
uint32_t group_id = (uint32_t)vl_seg_part.size();

// Assign group ids to new segments
for (const auto& x : s_seg_part)
{
    auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

    if (p == m_kmers.end())
        m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
}
```

### Counter/ID Assignment

**Starting ID**: `group_id = vl_seg_part.size()` (current number of groups)
**Increment**: `group_id++` for each NEW unique k-mer pair
**Result**: Monotonically increasing group IDs

### Adding to Groups

```cpp
// Lines 406-410
for (auto& x : s_seg_part)
{
    add_known(m_kmers[make_pair(x.kmer1, x.kmer2)], x.kmer1, x.kmer2,
        x.sample_name, x.contig_name, const_cast<contig_t&&>(x.seg_data),
        x.is_rev_comp, x.seg_part_no);
}
```

**Key Point**: All segments with the same `(kmer1, kmer2)` get the same `group_id`.

---

## 6. PROCESSING ORDER

**Location**: `agc_compressor.cpp`, lines 2244-2388, function `AddSampleFiles()`

### Priority Queue

**Type**: `CBoundedPQueue<task_t>` (bounded priority queue)

```cpp
// Lines 2244-2246
pq_contigs_desc = make_shared<CBoundedPQueue<task_t>>(1, queue_capacity);
pq_contigs_desc_aux = make_shared<CBoundedPQueue<task_t>>(1, ~0ull);
pq_contigs_desc_working = pq_contigs_desc;
```

### Priority Assignment

**Per-sample mode** (lines 2322-2353):
```cpp
size_t sample_priority = ~0ull;  // Start at maximum

// For each sample file:
for(auto sf : _v_sample_file_name)
{
    while (gio.ReadContigRaw(id, contig))
    {
        // Emplace contig with current priority
        pq_contigs_desc->Emplace(
            make_tuple(contig_processing_stage_t::all_contigs, sf.first, id, move(contig)),
            sample_priority, cost);
    }

    // Send synchronization token
    pq_contigs_desc->EmplaceManyNoCost(
        make_tuple(contig_processing_stage_t::registration, "", "", contig_t()),
        sample_priority, no_workers);

    --sample_priority;  // Decrease for next sample
}
```

**Priority order**: HIGHER priority value = processed FIRST (reverse numerical order)

### Queue Behavior

**Insertion** (queue.h, lines 238-244):
```cpp
void Emplace(T&& data, const size_t priority, const size_t cost)
{
    unique_lock<mutex> lck(mtx);
    cv_queue_full.wait(lck, [this] {return current_cost < max_cost; });

    bool was_empty = n_elements == 0;
    q.emplace(make_pair(priority, cost), move(data));  // multimap key
    ++n_elements;
    current_cost += cost;
}
```

**Extraction** (queue.h, lines 284-313):
```cpp
CBoundedPQueue::result_t PopLarge(T& data)
{
    unique_lock<mutex> lck(mtx);
    cv_queue_empty.wait(lck, [this] {return !this->q.empty() || !this->n_producers; });

    if (n_elements == 0)
        return n_producers ? result_t::empty : result_t::completed;

    data.swap(q.rbegin()->second);  // Get LARGEST (last) element
    size_t cost = q.rbegin()->first.second;

    q.erase(--q.end());
    --n_elements;
    current_cost -= cost;
}
```

**Key behavior**: `PopLarge()` gets the **LAST** element in the multimap, which is the **LARGEST** `(priority, cost)` pair.

### Processing Order

1. Contigs with **HIGHER priority** values are processed first
2. Within same priority, **LARGER cost** (contig size) processed first
3. Samples are processed in **insertion order** (first sample has highest priority)

---

## 7. BATCH BOUNDARIES

**Location**: `agc_compressor.cpp`, lines 2264-2367

### Synchronization Tokens

**Per-sample batching**:
```cpp
// Line 2349-2351
pq_contigs_desc->EmplaceManyNoCost(
    make_tuple(contig_processing_stage_t::registration, "", "", contig_t()),
    sample_priority, no_workers);
```

This sends `no_workers` (number of threads) synchronization tokens at the **end of each sample**.

### Registration Phase

**Location**: `agc_compressor.cpp`, lines 1137-1158

```cpp
if (get<0>(task) == contig_processing_stage_t::registration)
{
    // synchronization token received
    bar.arrive_and_wait();  // Barrier 1

    if (thread_id == 0)
        register_segments(n_t);  // Process new segments

    if (thread_id == 0)
        for (auto& v_fallback_minimizers : vv_fallback_minimizers)
        {
            for (auto& x : v_fallback_minimizers)
                add_fallback_mapping(x[0], x[1], x[2], (bool)x[3]);
            v_fallback_minimizers.clear();
        }

    bar.arrive_and_wait();  // Barrier 2

    store_segments(zstd_cctx, zstd_dctx);  // Write segments

    bar.arrive_and_wait();  // Barrier 3

    if (thread_id == 0)
    {
        buffered_seg_part.clear(max(1u, n_t-1));  // Clear buffers
        ++processed_samples;
    }
}
```

### Order of Operations

1. **All threads arrive at barrier 1**
2. **Thread 0 only**: Call `register_segments()` which:
   - Calls `buffered_seg_part.process_new()` to assign group IDs
   - Registers new streams in archive
   - Resizes `v_segments` vector
3. **All threads arrive at barrier 2**
4. **All threads**: Call `store_segments()` to write segment data
5. **All threads arrive at barrier 3**
6. **Thread 0 only**: Clear buffers

**Critical**: Batch boundaries are **per-sample** (or per-`pack_cardinality` contigs in concatenated mode).

---

## 8. CRITICAL DETAILS

### 8.1 MISSING K-mer Handling

**Value**: `~0ull` (all bits set, max uint64_t)

**Special group**:
```cpp
// agc_compressor.cpp, line 358
map_segments[make_pair(~0ull, ~0ull)] = 0;
```

Group ID 0 is **reserved** for segments with both k-mers missing.

### 8.2 First vs Subsequent Samples

**No special handling** - all samples follow the same logic.

**However**: The first sample typically has **more NEW segments** because `map_segments` starts empty.

### 8.3 Synchronization

**Segment buffer access**:
```cpp
// agc_compressor.h, lines 320-330
void add_known(uint32_t group_id, ...)
{
    vl_seg_part[group_id].emplace(...);  // Has internal mutex
}

void add_new(...)
{
    lock_guard<mutex> lck(mtx);
    s_seg_part.emplace(...);
}
```

**Map updates** (during `store_segments()`):
```cpp
// agc_compressor.cpp, lines 1025-1045
seg_map_mtx.lock();

auto p = map_segments.find(make_pair(kmer1, kmer2));
if (p == map_segments.end())
    map_segments[make_pair(kmer1, kmer2)] = group_id;
else if (p->second > group_id)
    p->second = group_id;  // Keep LOWER group ID

if (kmer1 != ~0ull && kmer2 != ~0ull)
{
    map_segments_terminators[kmer1].push_back(kmer2);
    sort(map_segments_terminators[kmer1].begin(),
         map_segments_terminators[kmer1].end());

    if (kmer1 != kmer2)
    {
        map_segments_terminators[kmer2].push_back(kmer1);
        sort(map_segments_terminators[kmer2].begin(),
             map_segments_terminators[kmer2].end());
    }
}

seg_map_mtx.unlock();
```

### 8.4 Raw Groups

**Special handling** for first `no_raw_groups` groups:

```cpp
// agc_compressor.cpp, lines 1048-1053
if (group_id < (int)no_raw_groups)
{
    in_group_id = v_segments[group_id]->add_raw(seg_data, zstd_cctx, zstd_dctx);
}
else
    in_group_id = v_segments[group_id]->add(seg_data, zstd_cctx, zstd_dctx);
```

**Raw groups** use simpler compression (no LZ encoding).

---

## 9. EXACT ALGORITHM FLOW

### Step-by-step for a Single Sample

1. **Main thread reads FASTA file** → emplace contigs into priority queue with priority `~0ull - sample_idx`

2. **Worker threads** pop contigs from queue (largest priority first)

3. **For each contig** (`compress_contig()`):
   - Iterate through bases
   - Find splitter k-mers
   - Create segments between splitters
   - Call `add_segment()` for each segment

4. **In `add_segment()`**:
   - Determine grouping key `pk = (kmer1, kmer2)` based on which k-mers are present
   - Check `map_segments.find(pk)`
   - If **FOUND**: `buffered_seg_part.add_known(group_id, ...)`
   - If **NOT FOUND**: `buffered_seg_part.add_new(kmer1, kmer2, ...)`

5. **After all contigs in sample**, main thread emplace synchronization tokens

6. **Worker threads** receive synchronization token → arrive at barrier

7. **Thread 0** calls `register_segments()`:
   - `buffered_seg_part.sort_known(n_t)` - sort known segments
   - `buffered_seg_part.process_new()`:
     - Iterate through `s_seg_part` (NEW segments)
     - For each unique `(kmer1, kmer2)`, assign next available group ID
     - Move segments from `s_seg_part` to `vl_seg_part[group_id]`
   - Register new archive streams
   - `buffered_seg_part.distribute_segments(0, 0, no_raw_groups)` - distribute group 0 segments
   - `buffered_seg_part.restart_read_vec()` - prepare for reading

8. **All threads** arrive at barrier 2

9. **All threads** call `store_segments()`:
   - Get next group ID from `buffered_seg_part.get_vec_id()` (reverse order)
   - For each segment in group:
     - Create `CSegment` object if needed
     - Update `map_segments` with actual group ID
     - Add segment data to `CSegment`
     - Record placement in collection

10. **All threads** arrive at barrier 3

11. **Thread 0** clears buffers

12. **Repeat** for next sample

---

## 10. KEY DIFFERENCES FROM RAGC

Based on this analysis, potential divergence points:

1. **Grouping order**: C++ AGC processes segments in priority queue order, RAGC may use different order
2. **Batch boundaries**: C++ AGC groups per-sample, RAGC may group differently
3. **Group ID assignment**: C++ AGC uses `vl_seg_part.size()` as starting ID, RAGC may differ
4. **Missing k-mer handling**: Both should use `~0ull` but implementation may differ
5. **Segment splitting**: C++ AGC has complex splitting logic when initial lookup fails
6. **Map updates timing**: C++ AGC updates `map_segments` during `store_segments()`, not during `process_new()`

---

## CONCLUSION

C++ AGC's segment grouping is a **two-phase buffered system**:
1. **Buffering phase**: Segments classified as KNOWN (existing group) or NEW (needs group)
2. **Registration phase**: NEW segments get group IDs, all segments written to archive

The grouping key is strictly `(kmer1, kmer2)` where both are canonical k-mer values or `~0ull` for missing.

Group IDs are assigned **sequentially** starting from the current number of groups, ensuring **deterministic ordering** when processing segments in the same order.

The **priority queue** ensures samples are processed in **insertion order**, and **per-sample batching** ensures all segments from a sample are grouped before the next sample starts.
