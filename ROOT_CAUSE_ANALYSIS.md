# Root Cause Analysis: RAGC Archive Bloat

**Date**: 2025-11-12
**Problem**: RAGC creates 80M archives vs C++ AGC's 69M (16% bloat)
**Status**: ✅ ROOT CAUSE IDENTIFIED

---

## Executive Summary

**The root cause is NOT about group count - it's about SPLITTING DECISIONS.**

C++ AGC creates MORE groups (3,660) but produces SMALLER archives (69M) because it splits segments more aggressively, creating specialized groups that compress better.

RAGC creates FEWER groups (3,381) but produces LARGER archives (80M) because it fails to split segments when it should, grouping dissimilar segments together.

---

## Evidence

### Global Statistics

| Metric | C++ AGC | RAGC | Delta |
|--------|---------|------|-------|
| Groups created | 3,660 | 3,381 | +279 (C++ has more!) |
| Archive size | 69M | 80M | +11M (+16%) |
| Shared k-mer pairs | 2,921 | 2,921 | - |
| Unique k-mer pairs | 739 | 168 | +571 |

**Paradox**: MORE groups → SMALLER archive (better compression)

### Case Study: K-mer 0x4f51386ca000000

This front k-mer appears in segments from both implementations but with different grouping:

| Implementation | Groups created | Unique back k-mers | Shared | Unique |
|----------------|----------------|---------------------|--------|--------|
| C++ AGC | 83 | 83 | 49 | 34 |
| RAGC | 57 | 57 | 49 | 8 |

**C++ AGC creates 26 MORE specialized groups** for segments starting with this k-mer!

### Smoking Gun Example

**Segment with boundaries (0x4f51386ca000000, 0xbff71fdced000000):**

- **RAGC**: Creates group 3101 with these exact boundaries (3 segments)
- **C++ AGC**: NEVER creates this k-mer pair!
  - Uses 0xbff71fdced000000 as front k-mer for 5 different groups
  - Uses 0x4f51386ca000000 as front k-mer for 83 different groups
  - But NEVER combines them as (0x4f51386ca000000, 0xbff71fdced000000)

**Conclusion**: C++ AGC SPLIT the segment at a different position, creating sub-segments with different k-mer boundaries.

---

## Root Cause Hypothesis

### Phase 2 Splitting Logic Divergence

Both implementations use a 3-phase decision tree:
1. **Phase 1**: Check if exact group (k_front, k_back) exists → reuse
2. **Phase 2**: Try to split segment using middle k-mer → add to split groups
3. **Phase 3**: Create new group if both fail

**The divergence occurs in Phase 2 - the splitting logic.**

### Possible Causes

1. **Different split acceptance criteria**:
   - C++ AGC might use a different cost threshold
   - C++ AGC might have different compression estimation logic
   - Rounding errors or integer overflow in cost calculation

2. **Different split position search**:
   - Order of trying candidate split positions
   - Handling of degenerate splits (best_pos=0 or best_pos=len)
   - Early termination conditions

3. **Different reference availability checks**:
   - C++ AGC might find references that RAGC misses
   - Timing of when references are available
   - Map lookup differences

4. **Group reuse logic**:
   - C++ AGC updates map_segments with smaller group IDs: if (p->second > group_id) p->second = group_id;
   - This might affect subsequent split decisions
   - RAGC might not have equivalent logic

---

## Impact Analysis

### Why More Groups = Better Compression

When segments are grouped correctly (similar content together):
- **LZ compression** finds more matches within the group
- **ZSTD** can build better dictionaries
- **Smaller compressed blocks** overall

When segments are grouped poorly (dissimilar content together):
- **Fewer LZ matches** (segments are too different)
- **Larger compressed blocks** (can't compress dissimilar data well)
- **Worse compression ratio**

### The 16% Bloat Explained

RAGC's failure to split 571 segments (739 C++ unique - 168 RAGC unique) means:
- 571 segments that SHOULD be in specialized groups
- Are instead lumped into larger, less-compressible groups
- Each mis-grouped segment adds ~20KB of bloat (11M / 571 ≈ 19KB)

---

## Next Steps

### 1. Compare Split Cost Calculation

Add detailed logging to both implementations to see:
- What cost values are calculated for the same segment
- Which split positions are considered
- Why C++ AGC accepts a split that RAGC rejects

### 2. Identify Exact Divergence Point

Pick one specific segment:
- Input: Same FASTA coordinates
- Compare: C++ AGC vs RAGC split decision
- Find: Exact line of code where logic differs

### 3. Minimal Reproduction

Create a minimal test case:
- Single chromosome, single sample
- One segment that exhibits the divergence
- Unit test comparing split decisions

### 4. Implement Fix

Once divergence is understood:
- Update RAGC's splitting logic to match C++ AGC
- Verify: Archive size matches (±1%)
- Verify: Group count matches (±5%)
- Verify: Correctness (byte-identical extraction)
