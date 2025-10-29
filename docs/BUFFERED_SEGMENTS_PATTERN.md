# Buffered Segment Storage Pattern in C++ AGC

**Purpose**: Document C++ AGC's CBufferedSegPart structure for Rust implementation

**C++ Reference**:
- Definition: `agc_compressor.h` lines 27-536 (CBufferedSegPart)
- Usage: `agc_compressor.cpp` lines 954-971 (register_segments), 974-1093 (store_segments)

---

## C++ AGC CBufferedSegPart Structure

### Overview

CBufferedSegPart is a thread-safe buffer that:
1. **Collects segments** during contig compression (from multiple worker threads)
2. **Separates NEW vs KNOWN** segments
3. **Assigns group IDs** to NEW segments at registration barriers
4. **Distributes segments** to worker threads for storage

**Key insight**: This is the central data structure that enables streaming - segments are buffered during compression, then registered and stored at synchronization barriers.

---

## Data Structures

### seg_part_t (Private, lines 29-120)

```cpp
struct seg_part_t {
    uint64_t kmer1;           // First k-mer of segment
    uint64_t kmer2;           // Last k-mer of segment
    string sample_name;
    string contig_name;
    contig_t seg_data;        // Compressed segment data
    bool is_rev_comp;         // Is segment reverse complemented?
    uint32_t seg_part_no;     // Part number within contig

    // Ordering: by (sample_name, contig_name, seg_part_no)
    bool operator<(const seg_part_t& x) const;
};
```

**Purpose**: Internal storage for segments within `list_seg_part_t`

### kk_seg_part_t (Public, lines 124-165)

```cpp
struct kk_seg_part_t {
    uint64_t kmer1;
    uint64_t kmer2;
    string sample_name;
    string contig_name;
    contig_t seg_data;
    bool is_rev_comp;
    uint32_t seg_part_no;

    // Same fields as seg_part_t, used for NEW segments
    bool operator<(const kk_seg_part_t& x) const;
};
```

**Purpose**: Public-facing type for NEW segments in `s_seg_part` set

### list_seg_part_t (Private, lines 169-296)

```cpp
struct list_seg_part_t {
    mutex mtx;                    // Thread-safe access
    vector<seg_part_t> l_seg_part;
    size_t virt_begin = 0;        // Virtual start for pop operations

    void emplace(uint64_t kmer1, uint64_t kmer2, ...);  // Thread-safe add
    void sort();                  // Sort by (sample, contig, part_no)
    bool pop(seg_part_t& seg_part);  // Pop from virt_begin
    void clear();                 // Reset virt_begin and clear vector
};
```

**Purpose**: Thread-safe list of segments for a single group_id

---

## Main Data Members (lines 298-303)

```cpp
vector<list_seg_part_t> vl_seg_part;  // KNOWN segments, indexed by group_id
set<kk_seg_part_t> s_seg_part;        // NEW segments (no group_id yet)
mutex mtx;                             // Protects s_seg_part and vl_seg_part resizing
atomic<int32_t> a_v_part_id;          // Atomic counter for reading segments
```

**Architecture**:
- `vl_seg_part[group_id]` contains all segments for that group
- `s_seg_part` contains NEW segments (kmer pairs not seen before)
- During `register_segments()`, NEW segments are assigned group IDs and moved to `vl_seg_part`

---

## Key Methods

### Constructor (lines 308-311)

```cpp
CBufferedSegPart(uint32_t no_raw_groups) {
    vl_seg_part.resize(no_raw_groups);
}
```

**no_raw_groups**: Initial number of group IDs (expands as new segments found)

### add_known (lines 320-324)

```cpp
void add_known(uint32_t group_id, uint64_t kmer1, uint64_t kmer2,
               const string& sample_name, const string& contig_name,
               contig_t&& seg_data, bool is_rev_comp, uint32_t seg_part_no) {
    vl_seg_part[group_id].emplace(...);  // Thread-safe (internal mutex)
}
```

**Usage**: Add segment to KNOWN group during compression
**Thread-safety**: `list_seg_part_t::emplace()` has internal mutex

### add_new (lines 326-331)

```cpp
void add_new(uint64_t kmer1, uint64_t kmer2,
             const string& sample_name, const string& contig_name,
             contig_t& seg_data, bool is_rev_comp, uint32_t seg_part_no) {
    lock_guard<mutex> lck(mtx);
    s_seg_part.emplace(kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no);
}
```

**Usage**: Add NEW segment (kmer pair not yet assigned a group_id)
**Thread-safety**: Locks `mtx` to protect `s_seg_part` set

### sort_known (lines 333-377)

```cpp
void sort_known(uint32_t nt) {
    lock_guard<mutex> lck(mtx);

    // Parallel sort of all vl_seg_part lists
    // Each worker thread processes job_step groups
    auto job = [&] {
        while (true) {
            uint64_t j_from = seg_part_id.fetch_add(job_step);
            if (j_from >= n_seg_part) break;

            for (uint64_t j = j_from; j < j_to; ++j)
                vl_seg_part[j].sort();
        }
    };

    // Launch nt threads
    for (uint64_t i = 0; i < nt - 1; ++i)
        v_fut.emplace_back(async(job));
    job();  // Main thread also works

    for (auto& f : v_fut)
        f.wait();
}
```

**Usage**: Called during `register_segments()` **before** `process_new()`
**Purpose**: Sort all KNOWN segments by (sample, contig, part_no) for deterministic storage
**Parallelism**: Divides work among `nt` threads

### process_new (lines 384-415)

```cpp
uint32_t process_new() {
    lock_guard<mutex> lck(mtx);

    map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
    uint32_t group_id = (uint32_t)vl_seg_part.size();

    // Assign group IDs to NEW kmer pairs
    for (const auto& x : s_seg_part) {
        auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));
        if (p == m_kmers.end())
            m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
    }

    uint32_t no_new = group_id - (uint32_t)vl_seg_part.size();

    // Resize vl_seg_part to accommodate new groups
    vl_seg_part.resize(group_id);

    // Move NEW segments to KNOWN groups
    for (auto& x : s_seg_part) {
        add_known(m_kmers[make_pair(x.kmer1, x.kmer2)],
                  x.kmer1, x.kmer2, x.sample_name, x.contig_name,
                  const_cast<contig_t&&>(x.seg_data), x.is_rev_comp, x.seg_part_no);
    }

    s_seg_part.clear();

    return no_new;  // Number of NEW groups created
}
```

**Usage**: Called during `register_segments()` **after** `sort_known()`
**Returns**: Number of NEW segment groups created
**Side effects**:
- Assigns group IDs to all NEW kmer pairs
- Resizes `vl_seg_part` to add new groups
- Moves all segments from `s_seg_part` to `vl_seg_part`
- Clears `s_seg_part`

### distribute_segments (lines 417-435)

```cpp
void distribute_segments(uint32_t src_id, uint32_t dest_id_from, uint32_t dest_id_to) {
    uint32_t no_in_src = vl_seg_part[src_id].size();
    uint32_t dest_id_curr = dest_id_from;
    seg_part_t seg_part;

    for (uint32_t i = 0; i < no_in_src; ++i) {
        if (dest_id_curr != src_id) {
            vl_seg_part[src_id].pop(seg_part);
            vl_seg_part[dest_id_curr].append_no_lock(seg_part);
        }

        if (++dest_id_curr == dest_id_to)
            dest_id_curr = dest_id_from;
    }
}
```

**Usage**: Called during `register_segments()` **after** `process_new()`
**Purpose**: Distribute segments from group 0 to worker threads
**Pattern**: Round-robin distribution

**Example**: `distribute_segments(0, 0, num_workers)`
- Distributes group 0 segments among workers 0 to num_workers-1
- Worker 0 keeps every num_workers-th segment
- Workers 1..N-1 receive remaining segments

### clear (lines 461-507)

```cpp
void clear(uint32_t nt) {
    lock_guard<mutex> lck(mtx);

    // Parallel clear of all vl_seg_part lists
    auto job = [&] {
        while (true) {
            uint64_t loc_idx = idx.fetch_add(job_step);
            if (loc_idx >= n_seg_part) break;

            for (; loc_idx < upp_idx; ++loc_idx)
                vl_seg_part[loc_idx].clear();
        }
    };

    // Launch nt threads
    for (uint64_t i = 0; i < nt - 1; ++i)
        v_fut.emplace_back(async(job));

    s_seg_part.clear();
    job();

    for (auto& f : v_fut)
        f.wait();
}
```

**Usage**: Called after `store_segments()` completes
**Purpose**: Clear all buffered segments for next round
**Parallelism**: Divides work among `nt` threads

### Reading Methods (lines 509-535)

```cpp
void restart_read_vec() {
    lock_guard<mutex> lck(mtx);
    a_v_part_id = (int32_t)vl_seg_part.size() - 1;
}

int get_vec_id() {
    return a_v_part_id.fetch_sub(part_id_step);  // Atomic decrement
}

bool is_empty_part(int group_id) {
    return group_id < 0 || vl_seg_part[group_id].empty();
}

bool get_part(int group_id, uint64_t& kmer1, uint64_t& kmer2,
              string& sample_name, string& contig_name,
              contig_t& seg_data, bool& is_rev_comp, uint32_t& seg_part_no) {
    return vl_seg_part[group_id].pop(kmer1, kmer2, sample_name, contig_name,
                                      seg_data, is_rev_comp, seg_part_no);
}
```

**Usage**: Used during `store_segments()` for parallel reading

**Pattern**:
1. `restart_read_vec()` - Set counter to max group_id
2. Worker threads call `get_vec_id()` atomically to claim group ranges
3. Worker processes all segments in group via `get_part()` loop

---

## Usage in register_segments (agc_compressor.cpp:954-971)

```cpp
void CAGCCompressor::register_segments(uint32_t n_t) {
    buffered_seg_part.sort_known(n_t);          // 1. Sort all KNOWN segments

    uint32_t no_new = buffered_seg_part.process_new();  // 2. Assign group IDs to NEW

    // 3. Register new segment streams with archive
    for (uint32_t i = 0; i < no_new; ++i)
        out_archive->RegisterStreams(
            ss_ref_name(archive_version, no_segments + i),
            ss_delta_name(archive_version, no_segments + i));

    no_segments += no_new;

    if (no_segments > v_segments.size())
        v_segments.resize(no_segments);

    // 4. Distribute group 0 to worker threads
    buffered_seg_part.distribute_segments(0, 0, no_raw_groups);

    // 5. Prepare for reading
    buffered_seg_part.restart_read_vec();
}
```

---

## Usage in store_segments (agc_compressor.cpp:990-1028)

```cpp
void CAGCCompressor::store_segments(ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx) {
    while (true) {
        int block_group_id = buffered_seg_part.get_vec_id();  // Atomic

        if (block_group_id < 0)
            break;

        // Process segments in reverse order (highest group_id first)
        for (int group_id = block_group_id; group_id > block_group_id - part_id_step; --group_id) {
            if (!buffered_seg_part.is_empty_part(group_id)) {
                while (buffered_seg_part.get_part(group_id, kmer1, kmer2,
                                                   sample_name, contig_name,
                                                   seg_data, is_rev_comp, seg_part_no)) {
                    // Create CSegment if needed
                    if (v_segments[group_id] == nullptr) {
                        v_segments[group_id] = make_shared<CSegment>(...);

                        // Update map_segments and map_segments_terminators
                        seg_map_mtx.lock();
                        map_segments[make_pair(kmer1, kmer2)] = group_id;
                        map_segments_terminators[kmer1].push_back(kmer2);
                        // ...
                        seg_map_mtx.unlock();
                    }

                    // Add segment to CSegment (compresses with ZSTD)
                    v_segments[group_id]->add(...);
                }
            }
        }
    }
}
```

**Key Pattern**:
- Workers atomically claim group ID ranges via `get_vec_id()`
- Each worker processes all segments in their claimed groups
- Segments are written to archive via `CSegment::add()`

---

## Rust Implementation Strategy

### Module Structure

```rust
// ragc-core/src/segment_buffer.rs

use std::collections::{HashMap, BTreeSet};
use std::sync::{Arc, Mutex, atomic::{AtomicI32, Ordering}};

/// Segment part data (matches seg_part_t)
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct SegmentPart {
    kmer1: u64,
    kmer2: u64,
    sample_name: String,
    contig_name: String,
    seg_data: Vec<u8>,
    is_rev_comp: bool,
    seg_part_no: u32,
}

/// Thread-safe list of segments for one group (matches list_seg_part_t)
struct SegmentPartList {
    parts: Mutex<Vec<SegmentPart>>,
    virt_begin: Mutex<usize>,
}

impl SegmentPartList {
    fn emplace(&self, part: SegmentPart);  // Thread-safe add
    fn sort(&self);                         // Sort by (sample, contig, part_no)
    fn pop(&self) -> Option<SegmentPart>;  // Pop from virt_begin
    fn clear(&self);
    fn is_empty(&self) -> bool;
}

/// Buffered segments storage (matches CBufferedSegPart)
pub struct BufferedSegments {
    /// KNOWN segments indexed by group_id
    vl_seg_part: Vec<SegmentPartList>,

    /// NEW segments (no group_id yet)
    s_seg_part: Mutex<BTreeSet<SegmentPart>>,

    /// Atomic counter for reading segments
    a_v_part_id: AtomicI32,

    /// Mutex for resizing vl_seg_part
    resize_mtx: Mutex<()>,
}

impl BufferedSegments {
    pub fn new(no_raw_groups: usize) -> Self;

    pub fn add_known(&self, group_id: u32, kmer1: u64, kmer2: u64, ...);
    pub fn add_new(&self, kmer1: u64, kmer2: u64, ...);

    pub fn sort_known(&self, num_threads: usize);
    pub fn process_new(&self) -> u32;
    pub fn distribute_segments(&self, src_id: u32, dest_from: u32, dest_to: u32);
    pub fn clear(&self, num_threads: usize);

    pub fn restart_read_vec(&self);
    pub fn get_vec_id(&self) -> i32;
    pub fn is_empty_part(&self, group_id: i32) -> bool;
    pub fn get_part(&self, group_id: i32) -> Option<(u64, u64, String, String, Vec<u8>, bool, u32)>;
}
```

---

## Testing Strategy

```rust
#[test]
fn test_buffered_segments_new_and_known() {
    let buf = BufferedSegments::new(10);

    // Add KNOWN segment to group 5
    buf.add_known(5, 100, 200, "sample1".into(), "chr1".into(),
                  vec![0, 1, 2, 3], false, 0);

    // Add NEW segment
    buf.add_new(300, 400, "sample1".into(), "chr1".into(),
                vec![4, 5, 6, 7], false, 1);

    // Process NEW segments
    let no_new = buf.process_new();
    assert_eq!(no_new, 1);  // One new group created
}

#[test]
fn test_buffered_segments_parallel_read() {
    let buf = Arc::new(BufferedSegments::new(100));

    // Add segments to multiple groups...

    buf.restart_read_vec();

    // Spawn worker threads
    let handles: Vec<_> = (0..4).map(|_| {
        let b = buf.clone();
        thread::spawn(move || {
            while let Some(group_id) = {
                let id = b.get_vec_id();
                if id >= 0 { Some(id) } else { None }
            } {
                while let Some(part) = b.get_part(group_id) {
                    // Process part...
                }
            }
        })
    }).collect();

    for h in handles {
        h.join().unwrap();
    }
}
```

---

## Implementation Status

- [C] **Studied**: C++ AGC CBufferedSegPart structure ✓
- [R] **Implementation**: Next step - rewrite segment_buffer.rs
- [✓] **Verified**: Will verify with unit tests

**Note**: This is a complex, performance-critical data structure. The Rust implementation must match the C++ threading patterns exactly.
