## Root Cause Found: Excessive Group Creation

**Problem:** 803 groups created, should be ~240 (matching splitter count)
**Result:** 6.5 segments/pack instead of 50 → poor LZ compression

**Analysis:**
- 5206 segments ÷ 50 (PACK_CARDINALITY) = ~104 ideal packs
- Actual: 803 packs → 7.7x more than ideal
- This means ~803 groups with 6.5 segments each on average

**One-kmer fix helped but not enough:**
- Before: 1537 groups
- After: 803 groups (48% reduction)
- Optimal: 240 groups (matching splitters)
- Still 3.3x too many groups!

## Root Cause: Missing Segment Splitting Feature

**FOUND!** C++ AGC analysis reveals the missing piece (see `CPP_AGC_GROUPING_ALGORITHM.md`).

### What C++ AGC Does (RAGC Doesn't)

When a segment has both k-mers but no group exists for (k1, k2):

1. **Find shared splitters** using set intersection:
   - k1's terminators: [kA, kB, k_middle, kC]
   - k2's terminators: [kD, k_middle, kE]
   - Shared: [k_middle]

2. **Split the segment** at the optimal position:
   - Part 1: (k1, k_middle) → joins existing group A
   - Part 2: (k_middle, k2) → joins existing group B
   - **No new group created!**

3. **Calculate optimal split position**:
   - Test LZ compression cost at each position
   - Find position that minimizes total compression

### Example Scenario

Existing groups:
- Group A: (k1, k_middle) with 50 segments
- Group B: (k_middle, k2) with 50 segments

New segment arrives with (k1, k2):

**C++ AGC**:
- Finds k_middle as shared splitter
- Splits segment → Part 1 to Group A, Part 2 to Group B
- **Result**: Still 2 groups ✓

**RAGC**:
- No group for (k1, k2) exists
- Creates **new Group C**: (k1, k2)
- **Result**: Now 3 groups ✗

**Impact on yeast10**:
- 240 splitters create a network of k-mer pairs
- Segments that could split create new groups instead
- **803 groups instead of ~240** (3.3x explosion)

### Why One-Kmer Fix Wasn't Enough

The one-kmer fix (commit 1f2b580) addressed **boundary segments** (segments with only one k-mer).

But the majority of segments have **both k-mers**!

These segments create new unique (k1, k2) pairs that:
- ✓ Should be split and joined to existing groups (C++ AGC)
- ✗ Create new groups instead (RAGC)

**This is why**: 1537 → 803 groups (48% reduction) but still 3.3x too many.

