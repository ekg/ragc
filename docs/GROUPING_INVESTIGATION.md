# Segment Grouping Investigation

**Date**: 2025-10-30
**Problem**: Same splitters (11,771), but different group counts (C++ AGC: 12,182, RAGC: 11,001)

---

## Hypothesis

If segmentation produces identical segments, then the difference must be in **grouping/deduplication**.

Segments are grouped by their (front_kmer, back_kmer) pairs. The key question:
- Are the normalization rules identical?
- Are there edge cases handled differently?
- Are segment sequences being compared when k-mers match?

---

## C++ AGC Grouping Logic

### Normalization (agc_compressor.cpp:1302-1313)

```cpp
if (kmer_front.is_full() && kmer_back.is_full()) {
    // Both terminal splitters present
    if (kmer_front.data() < kmer_back.data())
        pk = make_pair(kmer_front.data(), kmer_back.data());
    else {
        pk = make_pair(kmer_back.data(), kmer_front.data());
        reverse_complement_copy(segment, segment_rc);
        store_rc = true;
    }
}
```

**Rule**: When both k-mers present, normalize so `pk.first <= pk.second`, and reverse-complement if swapped.

### One K-mer Present (lines 1315-1361)

```cpp
else if (kmer_front.is_full()) {
    // Only front k-mer present
    tie(pk, store_rc) = find_cand_segment_with_one_splitter(kmer, segment, segment_rc, ...);
}
else if (kmer_back.is_full()) {
    // Only back k-mer present
    tie(pk, store_rc) = find_cand_segment_with_one_splitter(kmer, segment_rc, segment, ...);
}
```

**Complex logic**: Tries to find existing segments with one matching k-mer.

### No K-mers Present (lines 1286-1300)

```cpp
if (!kmer_front.is_full() && !kmer_back.is_full()) {
    // Try fallback minimizers
    tie(pk, store_rc) = find_cand_segment_using_fallback_minimizers(segment, 1);
}
```

---

## RAGC Grouping Logic

### Normalization (compressor_streaming.rs:167-199)

```rust
fn new_normalized(kmer_front: u64, kmer_back: u64) -> (Self, bool) {
    if kmer_front != MISSING_KMER && kmer_back != MISSING_KMER {
        // Both present: normalize
        if kmer_front < kmer_back {
            (SegmentGroupKey { kmer_front, kmer_back }, false)
        } else {
            (SegmentGroupKey { kmer_front: kmer_back, kmer_back: kmer_front }, true)
        }
    } else {
        // One or both missing: don't normalize
        (SegmentGroupKey { kmer_front, kmer_back }, false)
    }
}
```

**Rule**: Same as C++ AGC for both-k-mers case. But for one-or-none k-mers, just uses them as-is.

---

## Key Difference Found

### C++ AGC: Complex one-k-mer handling
When only one k-mer is present, C++ AGC calls `find_cand_segment_with_one_splitter` which:
1. Searches for existing groups with matching k-mer
2. Compares segment sequences to find best match
3. May create a different group key than just (kmer, MISSING)

### RAGC: Simple one-k-mer handling
When one k-mer is MISSING, RAGC just uses:
```rust
SegmentGroupKey { kmer_front, kmer_back }  // One is MISSING_KMER
```

**This could create different groups!**

---

## Test Plan

1. Count how many segments have:
   - Both k-mers present
   - Only front k-mer
   - Only back k-mer
   - Neither k-mer

2. For single-k-mer segments, compare what C++ AGC does vs RAGC

3. Dump actual groups with their k-mer pairs from both implementations

4. Find specific segments that are grouped differently

---

## Next Steps

1. Add logging to both implementations to dump segment groups
2. Compare the group dumps
3. Focus on segments with MISSING_KMER values
