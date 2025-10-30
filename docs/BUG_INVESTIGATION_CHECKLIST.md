# Bug Investigation Checklist: Data Corruption in Sample 2+

**Date**: 2025-10-30
**Bug**: Samples after first are corrupted (30% data missing)
**Minimal Reproduction**: 2 samples (AAA#0 perfect, AAB#0 corrupted)

---

## Test Results

| Sample | C++ AGC Size | RAGC Size | Data Loss | Status |
|--------|--------------|-----------|-----------|--------|
| AAA#0 | 12,309,305 | 12,309,305 | 0% | ✓ PERFECT |
| AAB#0 | 12,299,861 | 8,584,313 | 30.2% | ✗ CORRUPTED |

**Key Statistics from 2-sample test**:
- RAGC: 22,446 segments collected
- RAGC: 12,910 NEW segments, 9,536 KNOWN segments
- RAGC: 24 splits executed in Phase 4 (KNOWN segments)
- C++ AGC: 31 splits executed total
- RAGC: 12,180 groups created

---

## Investigation Plan

### Phase 1: Segment Accounting ⬜
**Goal**: Verify every segment from input reaches the archive

**Checks**:
- [ ] Count segments extracted from input FASTAs
- [ ] Count segments in Phase 1 buffer
- [ ] Count segments in Phase 3 (NEW)
- [ ] Count segments in Phase 4 (KNOWN)
- [ ] Count segments written to archive
- [ ] Verify: input_count == buffer_count == NEW_count + KNOWN_count == written_count

**Expected**: All segments accounted for at every stage

---

### Phase 2: Phase 1 Verification (Segment Collection) ⬜
**Goal**: Verify Phase 1 creates correct segments with correct keys

**C++ AGC Logic**:
```cpp
compress_contig() {
    for (each k-mer in contig) {
        if (is_splitter(kmer)) {
            add_segment(segment_data, kmer_front, kmer_back);
        }
    }
    add_segment(final_segment);
}
```

**RAGC Implementation**: `compressor_streaming.rs:1248-1287`

**Checks**:
- [ ] Count contigs processed: AAA#0 = 17 contigs, AAB#0 = 17 contigs
- [ ] Count segments per contig (spot check a few)
- [ ] Verify k-mer extraction at segment boundaries
- [ ] Verify key normalization (canonical form)
- [ ] Compare first few segments from AAB#0 between RAGC and C++ AGC

**Test**: Add logging for first 10 segments of AAB#0

---

### Phase 3: Phase 2 Verification (Sorting & Classification) ⬜
**Goal**: Verify segments are correctly classified as NEW vs KNOWN

**C++ AGC Logic**:
- NEW segments: `map_segments.find(pk) == map_segments.end()`
- KNOWN segments: `map_segments.find(pk) != map_segments.end()`

**RAGC Implementation**: `compressor_streaming.rs:1266-1269`
```rust
let is_new = !known_groups.contains(&key);
```

**Checks**:
- [ ] After AAA#0: How many unique keys in `known_groups`?
- [ ] For AAB#0 segments: How many match existing keys?
- [ ] Expected: AAA#0 creates ~12K groups, AAB#0 reuses ~9.5K of them
- [ ] Verify: NEW + KNOWN = total segments

**Test**: Log classification for first 20 segments of AAB#0

---

### Phase 4: Phase 3 Verification (NEW Segment Processing) ⬜
**Goal**: Verify NEW segments create groups correctly

**C++ AGC Logic** (`store_segments()` lines 995-1049):
```cpp
if (kmer1 != ~0ull && kmer2 != ~0ull) {
    map_segments_terminators[kmer1].push_back(kmer2);
    if (kmer1 != kmer2)
        map_segments_terminators[kmer2].push_back(kmer1);
}
```

**RAGC Implementation**: `compressor_streaming.rs:1559-1582`

**Checks**:
- [ ] Count groups created in Phase 3: Expected ~12,910 (matches NEW segment count)
- [ ] Verify terminators added for each group
- [ ] Count total terminator connections: Expected matches log output
- [ ] Verify bidirectional connections (k1→k2 AND k2→k1)

**Test**: After Phase 3, dump terminator map size and verify it matches expectations

---

### Phase 5: Phase 4 Verification (KNOWN Segment Processing) ✓
**Goal**: Identify where KNOWN segments are lost

**RESULT**: ✓ ALL SEGMENTS ARE ADDED - NO SEGMENTS LOST

Test results:
- Known segments processed: 9,536
- Segments reaching normal processing: 9,512
- Segments that split: 24 → 48 segments
- Total segments added: 9,560 = 9,512 + 48
- Expected: 9,560
- **Accounting is PERFECT**

**Conclusion**: Bug is NOT in Phase 4. All KNOWN segments are correctly added to groups.

**C++ AGC Logic** (`add_segment()` lines 1275-1504):
```cpp
// If group exists
auto p = map_segments.find(pk);
if (p != map_segments.end()) {
    buffered_seg_part.add_known(segment_id, ...);
}
// If group doesn't exist, try to split
else if (can_split) {
    // Split into two segments
    buffered_seg_part.add_known(segment_id1, ...);
    buffered_seg_part.add_known(segment_id2, ...);
}
// Otherwise, buffer as new
else {
    buffered_seg_part.add_new(pk.first, pk.second, ...);
}
```

**RAGC Implementation**: `compressor_streaming.rs:1632-1869`

**Critical Question**: What happens to segments that DON'T split?

**Checks**:
- [ ] Count KNOWN segments: Expected 9,536
- [ ] Count segments that enter split logic: Expected ? (need to check)
- [ ] Count segments that successfully split: 24 (from logs)
- [ ] Count segments that fail split checks: ?
- [ ] **Count segments added after split logic**: Should be 9,536 - 24*2 + 24 = 9,512
  - (9,536 original - 48 segments that became splits + 24 non-split segments)

**Wait, this math is wrong. Let me recalculate:**
- Start with 9,536 KNOWN segments
- 24 segments split into 48 segments
- So we should write: 9,536 - 24 + 48 = 9,560 segments

**Checks continued**:
- [ ] Line 1644: Count segments where `!groups.contains_key(&buffered_seg.key)` is TRUE
- [ ] Line 1800: Verify `continue` is hit for successful splits
- [ ] Line 1819: Count segments that reach "Normal processing"
- [ ] Verify: Every KNOWN segment is either split OR added normally

**Test**: Add counters for each path through Phase 4

---

### Phase 6: Split Logic Deep Dive ⬜
**Goal**: Understand why only 24 splits execute vs 31 in C++ AGC

**From logs**:
- RAGC: "Both target groups exist: 36"
- RAGC: "Split position calculated: 36"
- RAGC: "Splits executed: 24"
- Difference: 12 splits calculated but not executed

**Possible causes**:
- [ ] Degenerate splits (left_size == 0 or right_size == 0)
- [ ] Splits that update key but don't `continue`
- [ ] Segments lost in split processing

**Check line 1695-1801**:
```rust
if left_size == 0 {
    buffered_seg.key = ...;  // Update key
    // NO continue - falls through to line 1819
} else if right_size == 0 {
    buffered_seg.key = ...;  // Update key
    // NO continue - falls through to line 1819
} else {
    // Split into two
    stats_split_executed += 1;
    continue;  // Skip normal processing
}
```

**Checks**:
- [ ] Count degenerate splits (left_size == 0): ?
- [ ] Count degenerate splits (right_size == 0): ?
- [ ] Count full splits (both > 0): Should be 24
- [ ] Verify degenerate splits are added with updated key
- [ ] Total: degenerate + full should equal 36

**Test**: Add logging for left_size and right_size

---

### Phase 7: Group Addition Logic ⬜
**Goal**: Verify segments are added to groups

**RAGC Code** (lines 1820-1868):
```rust
let mut group_entry = groups.entry(buffered_seg.key.clone()).or_insert_with(|| {
    // Add terminators
    if buffered_seg.key.kmer_front != MISSING_KMER && buffered_seg.key.kmer_back != MISSING_KMER {
        group_terminators.entry(...).push(...);
    }
    // Create group
    (gid, GroupWriter::new(...))
});

if let Some(pack) = group_entry.1.add_segment(buffered_seg.segment, &self.config)? {
    pack_tx.send(compressed_pack)?;
}
```

**Checks**:
- [ ] Does `groups.entry()` ever fail to insert?
- [ ] Does `add_segment()` ever return None without error?
- [ ] Does `pack_tx.send()` ever fail?
- [ ] Are there any early returns or continues that skip addition?

**Test**: Add assertions that every segment reaches `add_segment()`

---

### Phase 8: Writer Thread Verification ⬜
**Goal**: Verify writer thread receives and writes all packs

**RAGC Code** (lines 1307-1350):
```rust
let writer_handle = {
    while let Ok(pack) = pack_rx.recv() {
        // Write pack
        packs_written += 1;
    }
};
```

**Checks**:
- [ ] Count packs sent: Should be logged
- [ ] Count packs received by writer: Should be logged as "23079 packs written"
- [ ] Verify: packs_sent == packs_received
- [ ] Check for any dropped packs

**Test**: Compare pack counts

---

### Phase 9: Collection Metadata ⬜
**Goal**: Verify collection correctly tracks segments per contig

**C++ AGC**: Collection stores segment count per contig for decompression

**Checks**:
- [ ] Count segments registered for AAA#0 contigs
- [ ] Count segments registered for AAB#0 contigs
- [ ] Verify: registered_count == segments_in_archive
- [ ] Check if collection has correct segment counts

**Test**: Extract collection metadata and compare

---

### Phase 10: Decompression Path ⬜
**Goal**: Verify decompression correctly reconstructs genome

**Checks**:
- [ ] Does decompressor read all segments?
- [ ] Does decompressor assemble segments in correct order?
- [ ] Are there any segments that fail to decompress?

**Test**: Add logging to decompression

---

## Hypotheses

### Hypothesis 1: Segments not added in Phase 4 ❓
**Theory**: Line 1644 check `!groups.contains_key(&buffered_seg.key)` is FALSE for most KNOWN segments, so they skip split logic AND normal addition

**Why this could happen**:
- In Phase 3, we create groups for NEW segments
- In Phase 4, KNOWN segments have keys that already exist
- Line 1644 checks if group DOESN'T exist - but it DOES exist (created in Phase 3)
- So most KNOWN segments skip the entire `if !groups.contains_key()` block
- They don't enter split logic
- They don't reach normal processing at line 1819
- **They are never added!**

**Test**: Count how many KNOWN segments have `groups.contains_key(&buffered_seg.key) == true`

---

### Hypothesis 2: Degenerate splits not counted ❓
**Theory**: Degenerate splits (left_size == 0 or right_size == 0) update key but don't add segment

**Why this could happen**:
- Line 1695-1700: Update key but DON'T continue
- Fall through to line 1819: Normal processing
- But at line 1820, we use `buffered_seg.key` (the UPDATED key)
- If this key already exists in groups, `or_insert_with` won't create new group
- Segment gets added to existing group
- But maybe with wrong metadata?

**Test**: Check if degenerate splits are added correctly

---

### Hypothesis 3: Phase 3 groups have wrong keys ❓
**Theory**: NEW segments create groups with incorrect keys, causing KNOWN segments to fail lookup

**Why this could happen**:
- Key normalization bug
- RC orientation bug
- K-mer extraction bug

**Test**: Compare keys between RAGC and C++ AGC for same segments

---

## Current Status

✅ **BUG FOUND!**

**Root Cause**: When segments split in Phase 4, they create segments with duplicate `seg_part_no` values.

**Example** (chrXV):
- Original segments from Phase 1: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...]
- Split at part_no=5 → creates seg(5) and seg(6)
- Split at part_no=8 → creates seg(8) and seg(9)
- Original segments 6 and 9 still exist!
- Result: TWO segments with part_no=6, TWO with part_no=9

When `add_segment_placed()` registers these segments, the second registration **overwrites** the first, causing segments to be lost from the collection metadata.

**Evidence**:
```
[BUG_CHECK] Splitting segment AAB#0:AAB#0#chrXV part_no=5 into part_no=5 and part_no=6
[BUG_CHECK] Splitting segment AAB#0:AAB#0#chrXV part_no=8 into part_no=8 and part_no=9
```

**Fix Required**: Renumber all subsequent segments when a split occurs, or track cumulative offsets per contig.

**Problem with Current Fix Attempt**:
The renumbering logic cannot distinguish between:
- Segment created BY split at position N (part_no = N or N+1, should keep its number)
- Original segment from Phase 1 (part_no = N+1, should be renumbered to N+2)

Both have the same `seg_part_no` value!

**Real Solution**:
C++ AGC processes segments ONE AT A TIME during initial segmentation. When it calls `add_segment()` with seg_part_no=N and the segment splits, it immediately creates two segments and writes them. There's no "Phase 1" that pre-numbers all segments.

RAGC's architecture is different:
- Phase 1: Create ALL segments, number them 0, 1, 2, ...
- Phase 4: Some segments split

This architectural difference causes the duplicate part_no problem.

**Possible Solutions**:
1. Renumber segments DURING Phase 4 as splits occur (complex - need to track per-contig state)
2. Use temporary numbers for split segments, then renumber everything at the end
3. Change architecture to match C++ AGC (process contigs sequentially, not in phases)
