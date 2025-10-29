# Segment Registration and Storage Pattern (C++ AGC)

This document analyzes C++ AGC's segment registration and storage logic to guide Rust implementation.

**Purpose**: Understand how NEW segments are assigned group IDs and how ALL segments (NEW + KNOWN) are stored to the archive with proper metadata updates.

**Key Files**:
- `agc_compressor.cpp` lines 954-971 (register_segments)
- `agc_compressor.cpp` lines 974-1050 (store_segments)
- `agc_compressor.h` lines 384-415 (CBufferedSegPart::process_new)

---

## Overview

The registration and storage phase happens at synchronization barriers (ContigProcessingStage::Registration). It involves two main operations:

1. **register_segments()** - Called once by thread 0
   - Sort KNOWN segments
   - Assign group IDs to NEW segments
   - Register archive streams
   - Prepare for parallel storage

2. **store_segments()** - Called by ALL worker threads in parallel
   - Atomically read buffered segments
   - Create CSegment objects for new groups
   - Update map_segments and map_segments_terminators
   - Write segments to archive
   - Track collection metadata

**Critical Insight**: Despite the name "register_segments", the actual semantic registration (updating `map_segments` and `map_segments_terminators`) happens in **store_segments**, NOT register_segments.

---

## 1. register_segments(n_t)

**Source**: `agc_compressor.cpp` lines 954-971

```cpp
void CAGCCompressor::register_segments(uint32_t n_t)
{
    buffered_seg_part.sort_known(n_t);

    uint32_t no_new = buffered_seg_part.process_new();

    for (uint32_t i = 0; i < no_new; ++i)
        out_archive->RegisterStreams(ss_ref_name(archive_version, no_segments + i),
                                      ss_delta_name(archive_version, no_segments + i));

    no_segments += no_new;

    if (no_segments > v_segments.size())
        v_segments.resize(no_segments);

    buffered_seg_part.distribute_segments(0, 0, no_raw_groups);

    buffered_seg_part.restart_read_vec();
}
```

### Step-by-Step Analysis

#### Step 1: Sort KNOWN Segments
```cpp
buffered_seg_part.sort_known(n_t);
```
- Parallel sort (with n_t threads) of all KNOWN segments in `vl_seg_part`
- Ordering: (sample_name, contig_name, seg_part_no)
- **Why**: Ensures segments are written in canonical order for archive compatibility

#### Step 2: Process NEW Segments
```cpp
uint32_t no_new = buffered_seg_part.process_new();
```
- Assign group IDs to all NEW segments in `s_seg_part`
- Move NEW segments to KNOWN groups in `vl_seg_part`
- Returns count of newly created groups
- **Details**: See Section 2 below

#### Step 3: Register Archive Streams
```cpp
for (uint32_t i = 0; i < no_new; ++i)
    out_archive->RegisterStreams(ss_ref_name(archive_version, no_segments + i),
                                  ss_delta_name(archive_version, no_segments + i));
```
- For each new group, register two streams in the archive:
  - Reference stream: "seg-NN" or "seg_dNN" (depending on archive_version)
  - Delta stream: "seg-NN" or "seg_dNN" (delta encoding)
- Group IDs are sequential: `no_segments`, `no_segments + 1`, ...

#### Step 4: Update Segment Counter
```cpp
no_segments += no_new;
```
- Total groups in archive now includes the new groups

#### Step 5: Resize v_segments
```cpp
if (no_segments > v_segments.size())
    v_segments.resize(no_segments);
```
- Ensure `v_segments` vector has capacity for all groups
- Each element is `shared_ptr<CSegment>` (initially nullptr)

#### Step 6: Distribute Raw Group Segments
```cpp
buffered_seg_part.distribute_segments(0, 0, no_raw_groups);
```
- Special handling for group 0 (raw segments)
- Round-robin distribution across first `no_raw_groups` groups
- **Note**: If `no_raw_groups == 0`, this is a no-op

#### Step 7: Restart Read Pointer
```cpp
buffered_seg_part.restart_read_vec();
```
- Reset atomic read counter `a_v_part_id` to enable parallel reading
- Workers will use `get_vec_id()` and `get_part()` to consume segments

---

## 2. CBufferedSegPart::process_new()

**Source**: `agc_compressor.h` lines 384-415

```cpp
uint32_t process_new()
{
    lock_guard<mutex> lck(mtx);

    map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
    uint32_t group_id = (uint32_t)vl_seg_part.size();

    // Assign group ids to new segments
    for (const auto& x : s_seg_part)
    {
        auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

        if (p == m_kmers.end())
            m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
    }

    uint32_t no_new = group_id - (uint32_t)vl_seg_part.size();

    if (vl_seg_part.capacity() < group_id)
        vl_seg_part.reserve((uint64_t)(group_id * 1.2));
    vl_seg_part.resize(group_id);

    for (auto& x : s_seg_part)
    {
        add_known(m_kmers[make_pair(x.kmer1, x.kmer2)], x.kmer1, x.kmer2,
            x.sample_name, x.contig_name, const_cast<contig_t&&>(x.seg_data),
            x.is_rev_comp, x.seg_part_no);
    }

    s_seg_part.clear();

    return no_new;
}
```

### Step-by-Step Analysis

#### Step 1: Create Local Group ID Map
```cpp
map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
uint32_t group_id = (uint32_t)vl_seg_part.size();
```
- `m_kmers`: Maps (kmer1, kmer2) → group_id
- Start group IDs from current size of `vl_seg_part` (all existing groups)

#### Step 2: Assign Group IDs to Unique (kmer1, kmer2) Pairs
```cpp
for (const auto& x : s_seg_part)
{
    auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

    if (p == m_kmers.end())
        m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
}
```
- Iterate through all NEW segments
- For each unique (kmer1, kmer2) pair, assign a new sequential group_id
- Multiple segments with same (kmer1, kmer2) will share the same group_id

**Example**:
```
NEW segments:
  - Segment A: (k1=10, k2=20) → group_id=100
  - Segment B: (k1=10, k2=20) → group_id=100 (same as A)
  - Segment C: (k1=15, k2=25) → group_id=101
  - Segment D: (k1=15, k2=25) → group_id=101 (same as C)

Result: 2 new groups (100, 101)
```

#### Step 3: Calculate Number of New Groups
```cpp
uint32_t no_new = group_id - (uint32_t)vl_seg_part.size();
```
- Count of unique (kmer1, kmer2) pairs found

#### Step 4: Resize vl_seg_part
```cpp
if (vl_seg_part.capacity() < group_id)
    vl_seg_part.reserve((uint64_t)(group_id * 1.2));
vl_seg_part.resize(group_id);
```
- Reserve 20% extra capacity to avoid frequent reallocations
- Resize to accommodate all new groups

#### Step 5: Move NEW Segments to KNOWN Groups
```cpp
for (auto& x : s_seg_part)
{
    add_known(m_kmers[make_pair(x.kmer1, x.kmer2)], x.kmer1, x.kmer2,
        x.sample_name, x.contig_name, const_cast<contig_t&&>(x.seg_data),
        x.is_rev_comp, x.seg_part_no);
}
```
- For each NEW segment, look up its assigned group_id in `m_kmers`
- Call `add_known()` to add to the appropriate `vl_seg_part[group_id]`

#### Step 6: Clear NEW Segments
```cpp
s_seg_part.clear();
```
- All NEW segments have been moved to KNOWN groups
- Clear the NEW segment set

#### Step 7: Return Count
```cpp
return no_new;
```

---

## 3. store_segments(zstd_cctx, zstd_dctx)

**Source**: `agc_compressor.cpp` lines 974-1050

```cpp
void CAGCCompressor::store_segments(ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    string sample_name;
    string contig_name;
    contig_t seg_data;
    bool is_rev_comp;
    uint32_t seg_part_no;
    uint64_t kmer1;
    uint64_t kmer2;

    const size_t max_buff_size = 32;

    vector<segments_to_place_t> buffered_coll_insertions;

    int no_parts = buffered_seg_part.get_no_parts();

    while (true)
    {
        int block_group_id = buffered_seg_part.get_vec_id();
        int in_group_id;

        if (block_group_id < 0)
            break;

        for (int group_id = block_group_id; group_id > block_group_id - CBufferedSegPart::part_id_step; --group_id)
        {
            if (!buffered_seg_part.is_empty_part(group_id))
                while (buffered_seg_part.get_part(group_id, kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no))
                {
                    if (v_segments[group_id] == nullptr)
                    {
                        v_segments[group_id] = make_shared<CSegment>(ss_base(archive_version, group_id), nullptr, out_archive, pack_cardinality, min_match_len, concatenated_genomes, archive_version);

                        seg_map_mtx.lock();

                        auto p = map_segments.find(make_pair(kmer1, kmer2));
                        if (p == map_segments.end())
                            map_segments[make_pair(kmer1, kmer2)] = group_id;
                        else if (p->second > group_id)
                            p->second = group_id;

                        if (kmer1 != ~0ull && kmer2 != ~0ull)
                        {
                            map_segments_terminators[kmer1].push_back(kmer2);
                            sort(map_segments_terminators[kmer1].begin(), map_segments_terminators[kmer1].end());

                            if (kmer1 != kmer2)
                            {
                                map_segments_terminators[kmer2].push_back(kmer1);
                                sort(map_segments_terminators[kmer2].begin(), map_segments_terminators[kmer2].end());
                            }
                        }

                        seg_map_mtx.unlock();
                    }

                    if (group_id < (int)no_raw_groups)
                    {
                        in_group_id = v_segments[group_id]->add_raw(seg_data, zstd_cctx, zstd_dctx);
                    }
                    else
                        in_group_id = v_segments[group_id]->add(seg_data, zstd_cctx, zstd_dctx);

                    if (buffered_coll_insertions.size() == max_buff_size)
                    {
                        collection_desc->add_segments_placed(buffered_coll_insertions);
                        buffered_coll_insertions.clear();
                    }

                    buffered_coll_insertions.emplace_back(sample_name, contig_name, seg_part_no, group_id, in_group_id, is_rev_comp, (uint32_t)seg_data.size());
                }
        }
    }

    collection_desc->add_segments_placed(buffered_coll_insertions);
}
```

### Step-by-Step Analysis

#### Step 1: Declare Local Variables
```cpp
string sample_name;
string contig_name;
contig_t seg_data;
bool is_rev_comp;
uint32_t seg_part_no;
uint64_t kmer1;
uint64_t kmer2;
```
- Reused for each segment to avoid allocations

#### Step 2: Atomic Block Iteration
```cpp
while (true)
{
    int block_group_id = buffered_seg_part.get_vec_id();
    int in_group_id;

    if (block_group_id < 0)
        break;
```
- `get_vec_id()` atomically returns next block of `part_id_step` groups (typically 64)
- Returns -1 when no more blocks available
- **Parallel-safe**: Multiple workers can call this concurrently

#### Step 3: Process Group Block
```cpp
for (int group_id = block_group_id; group_id > block_group_id - CBufferedSegPart::part_id_step; --group_id)
{
    if (!buffered_seg_part.is_empty_part(group_id))
```
- Iterate backwards through the block (group_id, group_id-1, ..., group_id-63)
- Skip empty groups

#### Step 4: Process Each Segment in Group
```cpp
while (buffered_seg_part.get_part(group_id, kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no))
{
```
- `get_part()` pops one segment from `vl_seg_part[group_id]`
- Uses virt_begin optimization (no actual removal, just increment pointer)
- Returns false when group is exhausted

#### Step 5: Create CSegment on First Use
```cpp
if (v_segments[group_id] == nullptr)
{
    v_segments[group_id] = make_shared<CSegment>(ss_base(archive_version, group_id), nullptr, out_archive, pack_cardinality, min_match_len, concatenated_genomes, archive_version);
```
- Lazy initialization: only create CSegment when first segment for this group is processed
- CSegment handles all LZ encoding and ZSTD compression for this group

#### Step 6: Update map_segments (CRITICAL)
```cpp
seg_map_mtx.lock();

auto p = map_segments.find(make_pair(kmer1, kmer2));
if (p == map_segments.end())
    map_segments[make_pair(kmer1, kmer2)] = group_id;
else if (p->second > group_id)
    p->second = group_id;
```
- **This is the actual "registration" of segment keys**
- Maps (kmer1, kmer2) → group_id
- If key already exists, keep the LOWER group_id (earlier samples have priority)
- **Mutex-protected**: Multiple workers may try to register same key

**Example**:
```
Sample 1, Contig A, Segment 1: (k1=10, k2=20) → group_id=100
Sample 2, Contig B, Segment 1: (k1=10, k2=20) → group_id=105

map_segments[(10, 20)] = 100 (keep lower group_id)

Later compression of Sample 2 will match against group 100 (from Sample 1)
```

#### Step 7: Update map_segments_terminators (CRITICAL)
```cpp
if (kmer1 != ~0ull && kmer2 != ~0ull)
{
    map_segments_terminators[kmer1].push_back(kmer2);
    sort(map_segments_terminators[kmer1].begin(), map_segments_terminators[kmer1].end());

    if (kmer1 != kmer2)
    {
        map_segments_terminators[kmer2].push_back(kmer1);
        sort(map_segments_terminators[kmer2].begin(), map_segments_terminators[kmer2].end());
    }
}
```
- **Terminator tracking for split detection**
- For each segment with both terminal k-mers:
  - Add k2 to the list of terminators for k1
  - Add k1 to the list of terminators for k2 (unless k1 == k2)
  - Keep lists sorted
- **Why**: Used in add_segment() to find middle splitters for split detection

**Example**:
```
Segment with (k1=10, k2=20):
  map_segments_terminators[10] = [..., 20]
  map_segments_terminators[20] = [..., 10]

Later, when processing NEW segment (k1=10, k2=20):
  - Check: is k1 in map_segments_terminators? YES
  - Check: is k2 in map_segments_terminators? YES
  - Try to find middle splitter k_mid where:
    - k_mid is in map_segments_terminators[10]
    - k_mid is in map_segments_terminators[20]
    - (10, k_mid) exists in map_segments
    - (k_mid, 20) exists in map_segments
```

**Note**: The sort() after every push_back is inefficient for large lists. C++ AGC does this for simplicity. Rust implementation could use BTreeSet or insert in sorted order.

#### Step 8: Add Segment to Group
```cpp
if (group_id < (int)no_raw_groups)
{
    in_group_id = v_segments[group_id]->add_raw(seg_data, zstd_cctx, zstd_dctx);
}
else
    in_group_id = v_segments[group_id]->add(seg_data, zstd_cctx, zstd_dctx);
```
- **add_raw()**: For raw groups (no LZ encoding, just ZSTD)
- **add()**: For normal groups (LZ differential encoding + ZSTD)
- Returns `in_group_id`: index within the group

#### Step 9: Buffer Collection Metadata
```cpp
if (buffered_coll_insertions.size() == max_buff_size)
{
    collection_desc->add_segments_placed(buffered_coll_insertions);
    buffered_coll_insertions.clear();
}

buffered_coll_insertions.emplace_back(sample_name, contig_name, seg_part_no, group_id, in_group_id, is_rev_comp, (uint32_t)seg_data.size());
```
- Batch collection metadata updates (flush every 32 segments)
- **segments_to_place_t**: (sample, contig, part_no, group_id, in_group_id, is_rev_comp, size)
- **Why batch**: collection_desc operations may be synchronized, so batching reduces contention

#### Step 10: Final Flush
```cpp
collection_desc->add_segments_placed(buffered_coll_insertions);
```
- After loop, flush any remaining buffered insertions

---

## 4. Worker Thread Synchronization Flow

**Source**: `agc_compressor.cpp` lines 1114-1182 (registration stage handler)

```cpp
if (get<0>(task) == contig_processing_stage_t::registration)
{
    // Barrier 1: All workers arrive
    bar.arrive_and_wait();

    // Leader-only: register segments
    if (thread_id == 0)
        register_segments(n_t);

    // Leader-only: process fallback minimizers
    if(thread_id == 0)
        for (auto& v_fallback_minimizers : vv_fallback_minimizers)
        {
            for (auto& x : v_fallback_minimizers)
                add_fallback_mapping(x[0], x[1], x[2], (bool)x[3]);
            v_fallback_minimizers.clear();
        }

    // Barrier 2: Wait for registration complete
    bar.arrive_and_wait();

    // All workers: store segments
    store_segments(zstd_cctx, zstd_dctx);

    // Barrier 3: Wait for storage complete
    bar.arrive_and_wait();

    // Thread 0: cleanup and progress tracking
    if (thread_id == 0)
    {
        buffered_seg_part.clear(max(1u, n_t-1));

        if (n_t == 1)
        {
            // Update processed_samples counter
            // Flush archive buffers
        }

        if(adaptive_compression)
            pq_contigs_desc_working = pq_contigs_desc;
    }
    // Thread 1: cleanup and progress tracking (ALSO!)
    else if (thread_id == 1)
    {
        // Same logic as thread 0 (redundant, but that's C++ AGC)
    }

    // Barrier 4: All ready to continue
    bar.arrive_and_wait();

    continue;
}
```

### Barrier Sequence

1. **Barrier 1**: All workers arrive
2. **Leader work**: register_segments() + fallback minimizers
3. **Barrier 2**: Wait for registration
4. **All workers**: store_segments() (parallel)
5. **Barrier 3**: Wait for storage
6. **Thread 0 & 1**: cleanup and progress (both do same work!)
7. **Barrier 4**: Ready to continue

**Critical**: BOTH thread 0 and thread 1 do cleanup work (update processed_samples, flush buffers). This appears redundant but is the C++ AGC behavior. Rust should match this exactly.

---

## 5. Rust Implementation Strategy

### Data Structures

```rust
/// Segment placement metadata for collection
#[derive(Debug, Clone)]
pub struct SegmentPlacement {
    pub sample_name: String,
    pub contig_name: String,
    pub seg_part_no: u32,
    pub group_id: u32,
    pub in_group_id: u32,
    pub is_rev_comp: bool,
    pub seg_size: u32,
}

/// Shared state for registration and storage
pub struct SharedCompressorState {
    /// Buffered segments (KNOWN + NEW)
    pub buffered_segments: BufferedSegments,

    /// Segment objects indexed by group_id
    pub segments: Vec<Option<Arc<Mutex<Segment>>>>,

    /// Map: (kmer1, kmer2) → group_id
    pub map_segments: Mutex<HashMap<(u64, u64), u32>>,

    /// Map: kmer → Vec<kmer> (sorted terminators)
    pub map_segments_terminators: Mutex<HashMap<u64, Vec<u64>>>,

    /// Total number of groups
    pub no_segments: AtomicU32,

    /// Archive writer
    pub archive: Arc<Mutex<Archive>>,

    /// Collection metadata
    pub collection: Arc<Mutex<Collection>>,

    // ... other fields
}
```

### register_segments()

```rust
pub fn register_segments(
    shared: &SharedCompressorState,
    num_threads: usize,
) {
    // Step 1: Sort KNOWN segments in parallel
    shared.buffered_segments.sort_known(num_threads);

    // Step 2: Process NEW segments
    let no_new = shared.buffered_segments.process_new();

    // Step 3: Register archive streams
    let current_no_segments = shared.no_segments.load(Ordering::Relaxed);

    for i in 0..no_new {
        let group_id = current_no_segments + i;
        shared.archive.lock().unwrap().register_streams(
            segment_ref_name(group_id),
            segment_delta_name(group_id),
        );
    }

    // Step 4: Update segment counter
    shared.no_segments.fetch_add(no_new, Ordering::Relaxed);

    // Step 5: Resize segments vector
    let new_total = current_no_segments + no_new;
    let mut segments = shared.segments.lock().unwrap();
    if new_total as usize > segments.len() {
        segments.resize(new_total as usize, None);
    }
    drop(segments);

    // Step 6: Distribute raw group segments
    if shared.no_raw_groups > 0 {
        shared.buffered_segments.distribute_segments(0, 0, shared.no_raw_groups);
    }

    // Step 7: Restart read pointer
    shared.buffered_segments.restart_read_vec();
}
```

### store_segments()

```rust
pub fn store_segments(
    worker_id: usize,
    shared: &SharedCompressorState,
    zstd_pool: &ZstdBufferPool,
) {
    const MAX_BUFF_SIZE: usize = 32;
    let mut buffered_placements = Vec::with_capacity(MAX_BUFF_SIZE);

    loop {
        // Atomic block allocation
        let block_group_id = shared.buffered_segments.get_vec_id();

        if block_group_id < 0 {
            break;
        }

        // Process groups in block (backwards)
        for group_id in (block_group_id - PART_ID_STEP + 1..=block_group_id).rev() {
            if shared.buffered_segments.is_empty_part(group_id) {
                continue;
            }

            // Process all segments in this group
            while let Some(part) = shared.buffered_segments.get_part(group_id) {
                // Create Segment on first use
                {
                    let mut segments = shared.segments.lock().unwrap();
                    if segments[group_id as usize].is_none() {
                        let segment = Segment::new(
                            group_id,
                            Arc::clone(&shared.archive),
                            shared.pack_cardinality,
                            shared.min_match_len,
                            shared.concatenated_genomes,
                        );
                        segments[group_id as usize] = Some(Arc::new(Mutex::new(segment)));

                        // Update map_segments
                        let key = (part.kmer1, part.kmer2);
                        let mut map_segments = shared.map_segments.lock().unwrap();
                        map_segments.entry(key)
                            .and_modify(|existing| {
                                if group_id < *existing {
                                    *existing = group_id;
                                }
                            })
                            .or_insert(group_id);
                        drop(map_segments);

                        // Update map_segments_terminators
                        if part.kmer1 != MISSING_KMER && part.kmer2 != MISSING_KMER {
                            let mut terminators = shared.map_segments_terminators.lock().unwrap();

                            // Add k2 to k1's terminator list
                            let list1 = terminators.entry(part.kmer1).or_insert_with(Vec::new);
                            list1.push(part.kmer2);
                            list1.sort_unstable();

                            // Add k1 to k2's terminator list (if different)
                            if part.kmer1 != part.kmer2 {
                                let list2 = terminators.entry(part.kmer2).or_insert_with(Vec::new);
                                list2.push(part.kmer1);
                                list2.sort_unstable();
                            }
                        }
                    }
                }

                // Add segment to group
                let segment = Arc::clone(&shared.segments.lock().unwrap()[group_id as usize].as_ref().unwrap());
                let in_group_id = if group_id < shared.no_raw_groups {
                    segment.lock().unwrap().add_raw(&part.seg_data, zstd_pool)
                } else {
                    segment.lock().unwrap().add(&part.seg_data, zstd_pool)
                };

                // Buffer collection metadata
                buffered_placements.push(SegmentPlacement {
                    sample_name: part.sample_name,
                    contig_name: part.contig_name,
                    seg_part_no: part.seg_part_no,
                    group_id,
                    in_group_id,
                    is_rev_comp: part.is_rev_comp,
                    seg_size: part.seg_data.len() as u32,
                });

                // Flush if batch full
                if buffered_placements.len() == MAX_BUFF_SIZE {
                    shared.collection.lock().unwrap()
                        .add_segments_placed(&buffered_placements);
                    buffered_placements.clear();
                }
            }
        }
    }

    // Final flush
    if !buffered_placements.is_empty() {
        shared.collection.lock().unwrap()
            .add_segments_placed(&buffered_placements);
    }
}
```

### Worker Thread Registration Handler

```rust
ContigProcessingStage::Registration => {
    // Barrier 1: All workers arrive
    barrier.wait();

    // Leader-only: register segments
    if worker_id == 0 {
        register_segments(&shared, num_workers);
        process_fallback_minimizers(&shared);
    }

    // Barrier 2: Wait for registration
    barrier.wait();

    // All workers: store segments
    store_segments(worker_id, &shared, &zstd_pool);

    // Barrier 3: Wait for storage
    barrier.wait();

    // Thread 0 & 1: cleanup and progress
    if worker_id == 0 || worker_id == 1 {
        shared.buffered_segments.clear(num_workers.saturating_sub(1).max(1));

        if num_workers == 1 {
            // Update processed_samples
            // Flush archive buffers
        }

        if shared.config.adaptive_mode {
            // Switch queue (if needed)
        }
    }

    // Barrier 4: Ready to continue
    barrier.wait();
}
```

---

## 6. Key Implementation Notes

### Parallel Safety

1. **get_vec_id()**: Atomic fetch_add ensures no two workers get same block
2. **get_part()**: Uses virt_begin (atomic increment) for lock-free reading
3. **map_segments**: Mutex-protected, check-and-update pattern
4. **map_segments_terminators**: Mutex-protected, append-and-sort pattern

### Memory Considerations

1. **Lazy CSegment creation**: Only create when first segment arrives
2. **Batch collection updates**: Reduces mutex contention
3. **Reuse local variables**: Avoid allocations in hot loop

### C++ AGC Quirks to Match

1. **BOTH thread 0 and 1 do cleanup**: Redundant but must match exactly
2. **Sort after every push_back**: Inefficient but exact behavior
3. **Keep lower group_id**: When (k1, k2) already exists, prefer earlier sample

### Testing Strategy

1. **Unit tests**:
   - process_new() with various NEW segment patterns
   - Concurrent store_segments() calls
   - map_segments update with conflicts

2. **Integration tests**:
   - Full registration → storage cycle
   - Verify map_segments_terminators accumulation
   - Verify collection metadata matches C++ AGC

---

## 7. Summary

**register_segments** is primarily organizational:
- Sort existing segments
- Assign group IDs to NEW segments
- Register archive streams
- Prepare for parallel storage

**store_segments** does the heavy lifting:
- Parallel atomic iteration
- Lazy segment creation
- **Actual semantic "registration"** (map_segments + terminators)
- Write segments to archive
- Track collection metadata

**Critical for correctness**:
- Terminator accumulation enables split detection in later samples
- Lower group_id preference ensures earlier samples are reference
- Exact barrier synchronization prevents race conditions
