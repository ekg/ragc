# Bug: Poor Grouping Due to MISSING_KMER

**Date**: 2025-10-27
**Status**: Root cause identified, fix in progress

## Summary

RAGC produces files 56% larger than C++ AGC (15M vs 9.6M) despite finding identical splitters (10,989,085 singletons, 240 splitters). The issue is in segment grouping, not splitter finding.

## Root Cause

Segments are grouped by `(kmer_front, kmer_back)` pairs. First/last segments of contigs have `MISSING_KMER` (u64::MAX) for one k-mer, creating too many unique keys:

```
[GROUP] Creating group 0 with key: front=18446744073709551615, back=8086069375427390912
[GROUP] Creating group 1 with key: front=6079842337715376484, back=8086069375427390912
[GROUP] Creating group 2 with key: front=18446744073709551615, back=6079842337715376484
```

Groups 0, 1, 2 share k-mer 8086069375427390912 but can't be merged because of MISSING_KMER.

## Impact

- yeast10 dataset: 5,206 segments → 1,537 groups (avg 3.4 segments/group)
- Should be: ~240 groups with 50 segments each (matching splitter count)
- Result: 1 segment per pack → no LZ compression benefit → 56% size increase

## Data Integrity

✓ All data is correct (sequences match perfectly)
✓ Splitter finding is correct (matches C++ AGC exactly)
✓ Segmentation is correct (splits at splitters properly)
✗ Grouping creates too many small groups

## C++ AGC Solution (FOUND!)

C++ AGC uses `map_segments_terminators` to track k-mer pairings:

```cpp
// When creating group with both k-mers:
if (kmer1 != ~0ull && kmer2 != ~0ull) {
    map_segments_terminators[kmer1].push_back(kmer2);
    map_segments_terminators[kmer2].push_back(kmer1);
}

// When segment has only ONE k-mer:
find_cand_segment_with_one_splitter() {
    // Find all existing groups containing this k-mer
    auto candidates = map_segments_terminators[kmer];

    // Test compression with each candidate group's reference
    for (cand_kmer : candidates) {
        pk = (min(kmer, cand_kmer), max(kmer, cand_kmer));
        size = estimate_compression(segment, group[pk].reference);
        if (size < best_size) best_group = pk;
    }
    return best_group;
}
```

## RAGC Issue

In `ragc-core/src/compressor_streaming.rs` lines 1177-1191, the terminator tracking is **COMMENTED OUT** with "TEMPORARILY DISABLED"!

This means:
- Groups with both k-mers are created but not tracked
- Segments with one k-mer can't find candidate groups
- Each one-kmer segment creates a new unique group
- Result: 1537 groups instead of ~240

## Fix Applied (Partial)

✅ **Phase 1 Complete:**
1. Re-enabled `group_terminators` tracking
2. Implemented one-kmer candidate search
3. Groups: 1537 → 240 (matches C++ AGC!)
4. Splitters: 10,989,085 singletons, 240 splitters (correct!)

⚠️ **Remaining Issue:**
- File size: 15M vs C++ AGC's 9.6M (56% larger)
- Cause: Using **first candidate** instead of testing all candidates
- Packs: 1,537 (avg 3.4 segs/pack) instead of optimal ~104 packs
- Result: Poor grouping → unfilled packs → less LZ benefit

## Phase 2 Required

Implement proper candidate selection (matching C++ AGC):
```rust
// For each one-kmer segment:
let candidates = group_terminators.get(&present_kmer);
let mut best_group = None;
let mut best_size = segment.len();

for cand_kmer in candidates {
    let group_key = normalize(present_kmer, cand_kmer);
    let estimated_size = estimate_compression(segment, group[group_key].reference);
    if estimated_size < best_size {
        best_size = estimated_size;
        best_group = Some(group_key);
    }
}
```

This will lead to:
- Better group selection (segments join most similar groups)
- Fuller packs (fewer small groups)
- Better LZ compression (more similar segments per group)
- **Expected result: ~9.6M matching C++ AGC**
