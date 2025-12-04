# C++ AGC Segment Ordering Mechanism

**Date**: 2025-12-04
**Purpose**: Document the exact ordering mechanism used by C++ AGC to determine `in_group_id` assignment

---

## Overview

C++ AGC assigns `in_group_id` values (the order of segments within a group) based on **lexicographic sorting** of segments by `(sample_name, contig_name, seg_part_no)` **within each group** before writing to the archive.

This document traces the complete flow from segment creation to archive writing.

---

## Data Structures

### 1. `kk_seg_part_t` (agc_compressor.h:124-165)

Represents a segment with metadata, used in the NEW segments set.

```cpp
struct kk_seg_part_t {
    uint64_t kmer1;          // Front k-mer
    uint64_t kmer2;          // Back k-mer
    string sample_name;
    string contig_name;
    contig_t seg_data;       // Actual sequence data
    bool is_rev_comp;
    uint32_t seg_part_no;    // Segment index within contig

    // Lexicographic ordering
    bool operator<(const struct kk_seg_part_t& x) const {
        if (sample_name != x.sample_name)
            return sample_name < x.sample_name;
        if (contig_name != x.contig_name)
            return contig_name < x.contig_name;
        return seg_part_no < x.seg_part_no;
    }
};
```

### 2. `seg_part_t` (agc_compressor.h:29-120)

Same structure as `kk_seg_part_t` but used in vectors for KNOWN segments.

```cpp
struct seg_part_t {
    uint64_t kmer1;
    uint64_t kmer2;
    string sample_name;
    string contig_name;
    contig_t seg_data;
    bool is_rev_comp;
    uint32_t seg_part_no;

    // Same lexicographic ordering as kk_seg_part_t
    bool operator<(const struct seg_part_t& x) const {
        if (sample_name != x.sample_name)
            return sample_name < x.sample_name;
        if (contig_name != x.contig_name)
            return contig_name < x.contig_name;
        return seg_part_no < x.seg_part_no;
    }
};
```

### 3. `CBufferedSegPart` (agc_compressor.h:27-536)

Manages segments during batching:

```cpp
class CBufferedSegPart {
    // For NEW segments (not yet assigned to a group)
    set<kk_seg_part_t> s_seg_part;  // Automatically sorted by operator<

    // For KNOWN segments (already have a group_id)
    vector<list_seg_part_t> vl_seg_part;  // One vector per group

    struct list_seg_part_t {
        vector<seg_part_t> l_seg_part;  // Segments in this group

        void sort() {
            std::sort(l_seg_part.begin(), l_seg_part.end());
        }
    };
};
```

---

## Segment Flow

### Phase 1: Segment Collection

During contig processing, segments are classified as NEW or KNOWN:

**NEW segments** (first time seeing this k-mer pair):
- Added to `s_seg_part` set via `add_new()` (agc_compressor.h:326-330)
- Automatically sorted by set's `operator<` (lexicographic)

**KNOWN segments** (k-mer pair already seen):
- Added to `vl_seg_part[group_id]` vector via `add_known()` (agc_compressor.h:320-324)
- **NOT sorted yet** - appended in arrival order

### Phase 2: Batch Flush - `register_segments()`

Called at batch boundaries (agc_compressor.cpp:972-986):

```cpp
void CAGCCompressor::register_segments(uint32_t n_t)
{
    // STEP 1: Sort all known segments within each group
    buffered_seg_part.sort_known(n_t);          // LINE 974

    // STEP 2: Process new segments (assign group IDs)
    uint32_t no_new = buffered_seg_part.process_new();  // LINE 976

    // STEP 3: Register streams for new groups
    for (uint32_t i = 0; i < no_new; ++i)
        out_archive->RegisterStreams(...);

    // STEP 4: Distribute segments (for raw groups optimization)
    buffered_seg_part.distribute_segments(0, 0, no_raw_groups);
}
```

#### Step 1: `sort_known()` (agc_compressor.h:333-373)

**CRITICAL**: Sorts segments **within each group** before writing:

```cpp
void sort_known(uint32_t nt) {
    // Multi-threaded sorting
    auto job = [&] {
        while (true) {
            // Process batches of groups
            for (uint64_t j = j_from; j < j_to; ++j, ++p)
                p->sort();  // LINE 362: Sorts l_seg_part vector
        }
    };

    // Launch worker threads
    vector<thread> workers;
    for (uint32_t i = 0; i < nt; ++i)
        workers.emplace_back(job);

    for (auto& w : workers)
        w.join();
}
```

**Effect**: After this step, every `vl_seg_part[group_id].l_seg_part` vector is sorted by (sample_name, contig_name, seg_part_no).

#### Step 2: `process_new()` (agc_compressor.h:384-415)

Assigns group IDs to new segments:

```cpp
uint32_t process_new() {
    map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
    uint32_t group_id = (uint32_t)vl_seg_part.size();

    // Iterate through s_seg_part (set - already sorted!)
    for (const auto& x : s_seg_part) {  // LINE 392
        auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

        if (p == m_kmers.end())
            m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
    }

    // Resize vl_seg_part to accommodate new groups
    vl_seg_part.resize(group_id);

    // Move segments from s_seg_part to vl_seg_part
    for (auto& x : s_seg_part) {  // LINE 406
        add_known(m_kmers[make_pair(x.kmer1, x.kmer2)],
                  x.kmer1, x.kmer2, x.sample_name, x.contig_name,
                  const_cast<contig_t&&>(x.seg_data),
                  x.is_rev_comp, x.seg_part_no);
    }

    s_seg_part.clear();
    return no_new;
}
```

**Key observations:**
1. `s_seg_part` is a **set**, so iteration at line 392 and 406 is in **sorted order** (by operator<)
2. This means segments are added to `vl_seg_part` **in lexicographic order** when first processed
3. Group IDs are assigned to unique k-mer pairs in the order they appear in the sorted set

### Phase 3: Archive Writing - `store_segments()`

Called after `register_segments()` (agc_compressor.cpp:992-1068):

```cpp
void CAGCCompressor::store_segments(ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    // Iterate through groups in reverse order (high to low group_id)
    while (true) {
        int block_group_id = buffered_seg_part.get_vec_id();
        if (block_group_id < 0)
            break;

        for (int group_id = block_group_id;
             group_id > block_group_id - CBufferedSegPart::part_id_step;
             --group_id) {

            // Get segments from group (in sorted order!)
            while (buffered_seg_part.get_part(group_id, kmer1, kmer2,
                   sample_name, contig_name, seg_data, is_rev_comp, seg_part_no)) {

                // Create segment group if first segment
                if (v_segments[group_id] == nullptr) {
                    v_segments[group_id] = make_shared<CSegment>(...);
                    // Update map_segments...
                }

                // Add segment to archive - assigns in_group_id
                if (group_id < (int)no_raw_groups)
                    in_group_id = v_segments[group_id]->add_raw(seg_data, ...);
                else
                    in_group_id = v_segments[group_id]->add(seg_data, ...);

                // Record in collection descriptor
                collection_desc->add_segments_placed(..., in_group_id, ...);
            }
        }
    }
}
```

**Key observations:**
1. `get_part()` (agc_compressor.h:532-535) calls `vl_seg_part[group_id].pop()`
2. This retrieves segments **in the order they were sorted** by `sort_known()`
3. `in_group_id` is assigned sequentially as segments are added via `->add()` or `->add_raw()`
4. **Result**: `in_group_id` reflects the lexicographic order within each group

---

## Summary: in_group_id Assignment Order

**The `in_group_id` for segments in a group is determined by:**

1. **Primary**: Lexicographic order by `(sample_name, contig_name, seg_part_no)`
2. **Applied**: Per-group sorting via `sort_known()` before writing
3. **Assigned**: Sequentially during `store_segments()` as segments are written

**Example:**

Group 17 contains:
- AAA#0/chrI/seg0
- AAA#0/chrI/seg1
- BBT#0/chrI/seg4
- ZZZ#0/chrII/seg0

After `sort_known()`, they're ordered:
1. AAA#0/chrI/seg0 → in_group_id = 0
2. AAA#0/chrI/seg1 → in_group_id = 1
3. BBT#0/chrI/seg4 → in_group_id = 2
4. ZZZ#0/chrII/seg0 → in_group_id = 3

---

## Comparison with Initial RAGC Implementation

**RAGC (before fix needed):**
- Segments in `pending_batch_segments` are sorted before group assignment
- BUT: Within each group, segments may NOT be sorted before writing
- This could cause different `in_group_id` assignments

**To match C++ AGC, RAGC must:**
1. Collect all segments for a batch
2. Assign group IDs to segments (in any order)
3. **SORT segments within each group** by (sample_name, contig_name, place)
4. Write segments to archive in that sorted order

---

## Next Steps

1. Check if RAGC sorts segments within groups before writing
2. If not, add per-group sorting before `add_segment_to_pack()`
3. Verify that segments within each group are written in lexicographic order
