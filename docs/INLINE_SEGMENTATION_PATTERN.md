# Inline Segmentation and Split Detection in C++ AGC

**Purpose**: Document C++ AGC's inline segmentation during compress_contig for Rust implementation

**C++ Reference**:
- compress_contig: `agc_compressor.cpp` lines 2000-2054
- add_segment: `agc_compressor.cpp` lines 1275-1507
- find_cand_segment_with_missing_middle_splitter: lines 1510-1597

---

## Overview

**Key Insight**: C++ AGC does NOT write segments immediately. Instead:
1. **Segment** contig at splitters during `compress_contig()`
2. For each segment, **determine key** (k1, k2) in `add_segment()`
3. **Buffer segment** (NEW or KNOWN) in `buffered_seg_part`
4. **Write later** during registration barrier via `store_segments()`

This enables:
- Accumulation of terminators from earlier samples
- Split detection using accumulated terminators
- Deferred writing until all segments from sample are ready

---

## compress_contig Flow (lines 2000-2054)

### Main Loop: Scan for Splitters

```cpp
bool CAGCCompressor::compress_contig(
    contig_processing_stage_t contig_processing_stage,
    string sample_name,
    string id,
    contig_t& contig,
    ZSTD_CCtx* zstd_cctx,
    ZSTD_DCtx* zstd_dctx,
    uint32_t thread_id,
    my_barrier& bar
) {
    CKmer kmer(kmer_length, kmer_mode_t::canonical);

    uint64_t pos = 0;
    uint64_t split_pos = 0;
    CKmer split_kmer(kmer_length, kmer_mode_t::canonical);  // Empty initially
    uint32_t seg_part_no = 0;

    // Scan contig for splitters
    for (auto x : contig) {
        if (x >> 2)         // x > 3 (not A/C/G/T)
            kmer.Reset();
        else {
            kmer.insert_canonical(x);

            if (kmer.is_full()) {
                uint64_t d = kmer.data_canonical();

                // Check if this k-mer is a splitter
                if (bloom_splitters.check(d) && hs_splitters.check(d)) {
                    // Found splitter! Extract segment [split_pos, pos+1)
                    auto seg_id = add_segment(
                        sample_name, id, seg_part_no,
                        move(get_part(contig, split_pos, pos + 1 - split_pos)),
                        split_kmer,  // Front k-mer (empty for first segment)
                        kmer,        // Back k-mer (current splitter)
                        zstd_cctx, zstd_dctx, thread_id, bar
                    );

                    ++seg_part_no;

                    // If segment was split into 2 parts
                    if (seg_id.contains_second)
                        ++seg_part_no;

                    // Update for next segment
                    split_pos = pos + 1 - kmer_length;
                    split_kmer = kmer;
                    kmer.Reset();
                }
            }
        }

        ++pos;
    }

    // Check if contig failed to segment (adaptive mode)
    if (adaptive_compression &&
        contig_processing_stage == all_contigs &&
        split_kmer == CKmer(kmer_length, kmer_mode_t::canonical)) {  // No splitters found

        if (contig.size() >= segment_size)
            find_new_splitters(contig, thread_id);  // Add to vv_splitters

        return false;  // Signal failure - will be reprocessed
    }

    // Add final segment (from last splitter to end)
    if (split_pos < contig.size())
        add_segment(
            sample_name, id, seg_part_no,
            move(get_part(contig, split_pos, contig.size() - split_pos)),
            split_kmer,    // Front k-mer (last splitter)
            CKmer(kmer_length, kmer_mode_t::canonical),  // Back k-mer (empty - end of contig)
            zstd_cctx, zstd_dctx, thread_id, bar
        );

    return true;  // Success
}
```

**Key Points**:
- `split_kmer` tracks last splitter (empty for first segment)
- Segments include splitter k-mer at both ends (overlap of kmer_length)
- Returns `false` if no splitters found (adaptive mode) → contig reprocessed later
- Returns `true` on success

---

## add_segment Flow (lines 1275-1507)

### Step 1: Determine Segment Key (k1, k2)

```cpp
pair_segment_desc_t CAGCCompressor::add_segment(
    const string& sample_name,
    const string& contig_name,
    uint32_t seg_part_no,
    contig_t &&segment,
    CKmer kmer_front,  // Front terminal k-mer (may be empty)
    CKmer kmer_back,   // Back terminal k-mer (may be empty)
    ZSTD_CCtx* zstd_cctx,
    ZSTD_DCtx* zstd_dctx,
    uint32_t thread_id,
    my_barrier& bar
) {
    pair<uint64_t, uint64_t> pk, pk2(~0ull, ~0ull);
    contig_t segment_rc;
    bool store_rc = false;

    // Case 1: No terminal splitters (first or last segment of contig)
    if (!kmer_front.is_full() && !kmer_back.is_full()) {
        if (fallback_filter) {
            // Use fallback minimizers to find similar segment
            tie(pk, store_rc) = find_cand_segment_using_fallback_minimizers(segment, 1);
            if (pk != pk_empty && store_rc)
                reverse_complement_copy(segment, segment_rc);
        } else {
            pk = pk_empty;  // (~0ull, ~0ull)
        }
    }

    // Case 2: Both terminal splitters present
    else if (kmer_front.is_full() && kmer_back.is_full()) {
        // Canonical ordering: smaller k-mer first
        if (kmer_front.data() < kmer_back.data())
            pk = make_pair(kmer_front.data(), kmer_back.data());
        else {
            pk = make_pair(kmer_back.data(), kmer_front.data());
            reverse_complement_copy(segment, segment_rc);
            store_rc = true;
        }
    }

    // Case 3: Only front splitter present
    else if (kmer_front.is_full()) {
        reverse_complement_copy(segment, segment_rc);
        tie(pk, store_rc) = find_cand_segment_with_one_splitter(
            kmer_front, segment, segment_rc, zstd_dctx, bar);

        // Fallback if no match
        if (pk.first == ~0ull || pk.second == ~0ull) {
            tie(pk, store_rc) = find_cand_segment_using_fallback_minimizers(segment, 5);
        }
    }

    // Case 4: Only back splitter present
    else if (kmer_back.is_full()) {
        // Similar to Case 3, with direction swapping
        // ...
    }

    // At this point: pk = (k1, k2) identifying this segment

    auto p = map_segments.find(pk);

    // ... (split detection logic - see below) ...
}
```

**Key Points**:
- Four cases based on which terminal k-mers are present
- Canonical ordering: always `k1 < k2`, reverse complement if needed
- Fallback minimizers for segments without clear terminal k-mers
- `store_rc` tracks whether to store reverse complement

### Step 2: Split Detection (lines 1366-1467)

**Critical Section**: If pk not in map_segments, try to split

```cpp
// Check if segment can be split
if (!concatenated_genomes &&
    p == map_segments.end() &&
    pk.first != ~0ull && pk.second != ~0ull &&
    map_segments_terminators.count(pk.first) &&   // k1 exists as terminator
    map_segments_terminators.count(pk.second))    // k2 exists as terminator
{
    static size_t split_attempts = 0;
    split_attempts++;
    cerr << "CPP_SPLIT_ATTEMPT " << split_attempts
         << ": (" << pk.first << "," << pk.second << ")"
         << " sample=" << sample_name
         << " contig=" << contig_name
         << " size=" << segment.size() << endl;

    // Prepare reverse complement if needed
    if (segment_rc.empty())
        reverse_complement_copy(segment, segment_rc);

    // Find middle splitter shared by both terminators
    auto split_match = find_cand_segment_with_missing_middle_splitter(
        kmer1, kmer2,
        use_rc ? segment_rc : segment,
        use_rc ? segment : segment_rc,
        zstd_dctx, bar
    );

    if (split_match.first != ~0ull) {
        // Split successful!
        static size_t splits_executed = 0;
        splits_executed++;
        cerr << "CPP_SPLIT_EXEC " << splits_executed
             << ": middle=" << split_match.first
             << " split_pos=" << split_match.second << endl;

        uint32_t left_size = split_match.second;
        uint32_t right_size = segment.size() - split_match.second;

        if (left_size == 0) {
            // Entire segment belongs to right group (k_middle, k2)
            pk = minmax(split_match.first, kmer2.data());
            // ... adjust store_rc ...
        }
        else if (right_size == 0) {
            // Entire segment belongs to left group (k1, k_middle)
            pk = minmax(kmer1.data(), split_match.first);
            // ... adjust store_rc ...
        }
        else {
            // Split into 2 segments with overlap
            uint32_t seg2_start_pos = left_size - kmer_length / 2;
            segment2.assign(segment.begin() + seg2_start_pos, segment.end());
            segment.resize(seg2_start_pos + kmer_length);

            // Segment 1: (k1, k_middle)
            pk = make_pair(kmer_front.data(), split_match.first);
            segment_id = map_segments.at(pk);  // Must exist

            // Segment 2: (k_middle, k2)
            pk2 = make_pair(split_match.first, kmer_back.data());
            segment_id2 = map_segments.at(pk2);  // Must exist
        }
    }

    p = map_segments.find(pk);  // Re-check after potential split
}
```

**Split Conditions**:
1. NOT concatenated_genomes mode
2. Key (k1, k2) NOT in map_segments (NEW segment)
3. Both k1 and k2 are valid terminators (~0ull means invalid)
4. Both k1 and k2 exist in map_segments_terminators (seen before)

**Split Outcomes**:
- **No split**: left_size == 0 or right_size == 0 → entire segment to one group
- **Split**: left_size > 0 and right_size > 0 → two overlapping segments

### Step 3: Buffer Segment (lines 1491-1504)

```cpp
if (p == map_segments.end()) {
    // NEW segment - no existing group
    buffered_seg_part.add_new(
        pk.first, pk.second,
        sample_name, contig_name,
        store_rc ? segment_rc : segment,
        store_rc, seg_part_no
    );
} else {
    // KNOWN segment - group exists
    if (segment_id2 == -1)
        segment_id = p->second;  // Use found group_id

    buffered_seg_part.add_known(
        segment_id,  // group_id
        ~0ull, ~0ull,  // k1, k2 (not needed for KNOWN)
        sample_name, contig_name,
        store_rc ? move(segment_rc) : move(segment),
        store_rc, seg_part_no
    );

    // If split into 2 segments, add second segment
    if (segment_id2 >= 0)
        buffered_seg_part.add_known(
            segment_id2,
            ~0ull, ~0ull,
            sample_name, contig_name,
            store2_rc ? move(segment2_rc) : move(segment2),
            store2_rc, seg_part_no + 1
        );
}

return pair_segment_desc_t(
    segment_desc_t(segment_id, 0, store_rc, segment_size),
    segment_desc_t(segment_id2, 0, store2_rc, segment2_size),
    segment_id2 >= 0  // contains_second
);
```

---

## find_cand_segment_with_missing_middle_splitter (lines 1510-1597)

### Purpose: Find shared splitter k_middle for split

```cpp
pair<uint64_t, uint32_t> find_cand_segment_with_missing_middle_splitter(
    CKmer kmer_front,  // k1
    CKmer kmer_back,   // k2
    contig_t& segment_dir,
    contig_t& segment_rc,
    ZSTD_DCtx* zstd_dctx,
    my_barrier& bar
) {
    auto p_front = map_segments_terminators.find(kmer_front.data());
    auto p_back = map_segments_terminators.find(kmer_back.data());

    if (p_front == map_segments_terminators.end() ||
        p_back == map_segments_terminators.end())
        return make_pair(~0ull, 0);

    // Find shared terminators between k1 and k2
    vector<uint64_t> shared_splitters;
    set_intersection(
        p_front->second.begin(), p_front->second.end(),
        p_back->second.begin(), p_back->second.end(),
        back_inserter(shared_splitters)
    );

    if (shared_splitters.empty())
        return make_pair(~0ull, 0);

    // Try each shared splitter to find split position
    for (auto middle_kmer : shared_splitters) {
        // Check if groups (k1, k_middle) and (k_middle, k2) exist
        auto pk_left = minmax(kmer_front.data(), middle_kmer);
        auto pk_right = minmax(middle_kmer, kmer_back.data());

        auto p_left = map_segments.find(pk_left);
        auto p_right = map_segments.find(pk_right);

        if (p_left == map_segments.end() || p_right == map_segments.end())
            continue;  // Groups don't exist, try next

        // Find k_middle position in segment
        uint32_t split_pos = find_middle_splitter(
            middle_kmer, segment_dir, zstd_dctx, bar,
            v_segments[p_left->second], v_segments[p_right->second]
        );

        if (split_pos != ~0u)
            return make_pair(middle_kmer, split_pos);
    }

    return make_pair(~0ull, 0);  // No valid split found
}
```

**Key Steps**:
1. Get terminator lists for k1 and k2
2. Find shared terminators (intersection)
3. For each shared k_middle:
   - Check if groups (k1, k_middle) and (k_middle, k2) exist
   - Try to find k_middle position in segment
   - Return first match

---

## Rust Implementation Strategy

### Module Structure

```rust
// ragc-core/src/inline_compression.rs

use crate::segment_buffer::BufferedSegments;
use crate::task::Task;
use std::sync::{Arc, Barrier, Mutex};
use std::collections::HashMap;

/// Compress contig with inline segmentation
///
/// Returns: true if successful, false if needs adaptive splitters
pub fn compress_contig(
    task: &Task,
    worker_id: usize,
    zstd_ctx: &mut ZstdContext,
    barrier: &Arc<Barrier>,
    shared: &Arc<SharedCompressorState>,
) -> bool {
    // 1. Scan contig for splitters
    let mut segments = Vec::new();
    let mut split_pos = 0;
    let mut split_kmer: Option<u64> = None;
    let mut seg_part_no = 0;

    let kmer_scanner = KmerScanner::new(kmer_length);

    for (pos, base) in task.sequence.iter().enumerate() {
        if let Some(kmer) = kmer_scanner.insert(*base) {
            if is_splitter(kmer, &shared.bloom_splitters, &shared.hs_splitters) {
                // Extract segment [split_pos, pos+1)
                let segment = task.sequence[split_pos..pos+1].to_vec();

                let seg_ids = add_segment(
                    &task.sample_name,
                    &task.contig_name,
                    seg_part_no,
                    segment,
                    split_kmer,
                    Some(kmer),
                    worker_id,
                    shared,
                );

                seg_part_no += 1;
                if seg_ids.1.is_some() {
                    seg_part_no += 1;  // Split into 2 segments
                }

                split_pos = pos + 1 - kmer_length;
                split_kmer = Some(kmer);
            }
        }
    }

    // Check if no splitters found (adaptive mode)
    if shared.adaptive_mode && task.stage == AllContigs && split_kmer.is_none() {
        if task.sequence.len() >= shared.segment_size {
            find_new_splitters(&task.sequence, worker_id, shared);
        }
        return false;  // Needs reprocessing
    }

    // Add final segment
    if split_pos < task.sequence.len() {
        add_segment(
            &task.sample_name,
            &task.contig_name,
            seg_part_no,
            task.sequence[split_pos..].to_vec(),
            split_kmer,
            None,
            worker_id,
            shared,
        );
    }

    true
}

/// Add segment to buffered_seg_part
///
/// Returns: (segment_id, segment_id2) where segment_id2 is Some if split occurred
fn add_segment(
    sample_name: &str,
    contig_name: &str,
    seg_part_no: u32,
    segment: Vec<u8>,
    kmer_front: Option<u64>,
    kmer_back: Option<u64>,
    worker_id: usize,
    shared: &Arc<SharedCompressorState>,
) -> (Option<u32>, Option<u32>) {
    // 1. Determine key (k1, k2)
    let (key, segment_to_store, store_rc) = determine_segment_key(
        &segment, kmer_front, kmer_back, shared
    );

    // 2. Check if group exists
    let map_segments = shared.map_segments.lock().unwrap();
    let group_id = map_segments.get(&key).copied();
    drop(map_segments);

    // 3. Try split if NEW segment
    let (segment1, segment2, group_id1, group_id2) = if group_id.is_none() {
        try_split_segment(
            key, &segment_to_store, store_rc,
            sample_name, contig_name, worker_id, shared
        )
    } else {
        (segment_to_store, None, group_id, None)
    };

    // 4. Buffer segment(s)
    if let Some(gid) = group_id1 {
        shared.buffered_seg_part.add_known(
            gid, key.0, key.1,
            sample_name.to_string(), contig_name.to_string(),
            segment1, store_rc, seg_part_no
        );
    } else {
        shared.buffered_seg_part.add_new(
            key.0, key.1,
            sample_name.to_string(), contig_name.to_string(),
            segment1, store_rc, seg_part_no
        );
    }

    if let Some(segment2_data) = segment2 {
        shared.buffered_seg_part.add_known(
            group_id2.unwrap(),
            // ...
            seg_part_no + 1
        );
    }

    (group_id1, group_id2)
}

/// Try to split NEW segment using shared terminators
fn try_split_segment(
    key: (u64, u64),
    segment: &[u8],
    store_rc: bool,
    sample_name: &str,
    contig_name: &str,
    worker_id: usize,
    shared: &Arc<SharedCompressorState>,
) -> (Vec<u8>, Option<Vec<u8>>, Option<u32>, Option<u32>) {
    // Check split eligibility
    let terminators = shared.map_segments_terminators.lock().unwrap();

    if !terminators.contains_key(&key.0) || !terminators.contains_key(&key.1) {
        return (segment.to_vec(), None, None, None);
    }

    // Find shared terminators
    let shared_splitters = find_shared_terminators(key.0, key.1, &terminators);

    if shared_splitters.is_empty() {
        return (segment.to_vec(), None, None, None);
    }

    // Try each shared splitter
    for middle_kmer in shared_splitters {
        let pk_left = (key.0.min(middle_kmer), key.0.max(middle_kmer));
        let pk_right = (middle_kmer.min(key.1), middle_kmer.max(key.1));

        let map_segments = shared.map_segments.lock().unwrap();
        let gid_left = map_segments.get(&pk_left).copied();
        let gid_right = map_segments.get(&pk_right).copied();
        drop(map_segments);

        if gid_left.is_none() || gid_right.is_none() {
            continue;
        }

        // Find middle splitter position in segment
        if let Some(split_pos) = find_middle_splitter_position(
            middle_kmer, segment, shared
        ) {
            // Split segment with overlap
            let seg1 = segment[..split_pos + kmer_length].to_vec();
            let seg2 = segment[split_pos - kmer_length/2..].to_vec();

            return (seg1, Some(seg2), gid_left, gid_right);
        }
    }

    // No split found
    (segment.to_vec(), None, None, None)
}
```

---

## Testing Strategy

```rust
#[test]
fn test_compress_contig_with_splitters() {
    let shared = create_test_shared_state();

    // Add some known splitters
    shared.hs_splitters.insert(SPLITTER_KMER_1);
    shared.hs_splitters.insert(SPLITTER_KMER_2);

    let task = Task::new_contig(
        "sample1".into(),
        "chr1".into(),
        create_contig_with_splitters(SPLITTER_KMER_1, SPLITTER_KMER_2),
        ContigProcessingStage::AllContigs
    );

    let success = compress_contig(&task, 0, &mut zstd_ctx, &barrier, &shared);

    assert!(success);
    // Check that segments were buffered
    assert!(shared.buffered_seg_part.get_no_parts() > 0);
}

#[test]
fn test_split_detection() {
    let shared = create_test_shared_state();

    // Setup: k1 and k2 both exist as terminators
    shared.map_segments_terminators.insert(K1, vec![K_MID, K_OTHER]);
    shared.map_segments_terminators.insert(K2, vec![K_MID, K_OTHER2]);

    // Groups (K1, K_MID) and (K_MID, K2) exist
    shared.map_segments.insert((K1, K_MID), 10);
    shared.map_segments.insert((K_MID, K2), 11);

    // Try to add NEW segment with key (K1, K2)
    let (seg1, seg2, gid1, gid2) = try_split_segment(
        (K1, K2), &segment, false, "sample1", "chr1", 0, &shared
    );

    // Should split into 2 segments
    assert!(seg2.is_some());
    assert_eq!(gid1, Some(10));
    assert_eq!(gid2, Some(11));
}
```

---

## Implementation Status

- [C] **Studied**: C++ AGC compress_contig and add_segment ✓
- [R] **Implementation**: Next step - create inline_compression.rs
- [✓] **Verified**: Will verify with unit tests

**Note**: This is the most complex part of the streaming architecture. Split detection enables efficient encoding by reusing existing segment groups.
