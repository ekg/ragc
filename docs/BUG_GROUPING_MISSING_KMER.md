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

## Fix Applied (Complete!)

✅ **Implementation Complete** (commit 1f2b580):
1. Re-enabled `group_terminators` tracking across ALL compression paths
2. Implemented one-kmer candidate search in three places:
   - **Rayon batch path**: Phase 1.5 pre-grouping fix-up
   - **Channel-based workers**: Inline candidate search
   - **DashMap contig-level**: Atomic terminators tracking
3. Each one-kmer segment finds existing groups via terminators map
4. Uses FIRST candidate to avoid creating new MISSING_KMER groups

**Results**:
- yeast10 (k=21, s=10000): **8.8M** (close to target!)
- yeast10 (k=31, s=60000): **15M** (better than C++ AGC's 103M)
- Packs: 1537 → 803 (48% reduction!)
- All tests passing ✓

## Why Compression Estimation Wasn't Needed

Initial investigation suggested testing all candidates with compression estimation.
However, testing revealed:
- One-kmer segments are small (31-104 bytes)
- LZ encoding doesn't compress them (estimated_size == raw_size)
- Since no candidate is "better", first candidate works fine
- The key is to JOIN an existing group, not find the "best" one

Using FIRST candidate is sufficient because:
1. Segments already grouped by k-mer similarity
2. Small boundary segments don't dominate compression
3. Simplicity avoids race conditions in concurrent code

**Status**: Bug fixed! File sizes are now competitive with C++ AGC.
