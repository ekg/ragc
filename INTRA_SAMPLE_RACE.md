# Threading Impact Analysis (2025-11-11)

## Results Summary

### Baseline (Before Any Fixes)
- RAGC multi-threaded: 81M, **3,395 groups**
- C++ AGC multi-threaded: 69M, ~3,220 groups
- Gap: +12M (+17%), **+175 extra groups**

### After Barriers Implementation
- RAGC multi-threaded: 81M, **3,387 groups** (-8 groups)
- Gap: +12M (+17%), **+167 extra groups**
- **Conclusion**: Barriers had minimal effect

### After Atomic Group ID Allocation
- RAGC multi-threaded: 81M, **3,387 groups** (no change)
- Gap: +12M (+17%), **+167 extra groups**
- **Conclusion**: Atomic allocation had zero effect

### Single-Threaded Test
- RAGC single-threaded: 80M, **3,381 groups** (-6 groups from multi-threaded)
- C++ AGC single-threaded: 69M, ~3,220 groups (same as multi-threaded)
- Gap: +11M (+16%), **+161 extra groups**
- **Conclusion**: Threading contributes only 6 extra groups (negligible)

## Key Findings

**Total reduction from all threading fixes: -14 groups** (3,395 → 3,381)
- Barriers: -8 groups
- Single-threaded vs multi-threaded: -6 groups
- Atomic allocation: 0 groups

**Remaining bloat: +161 groups** (3,381 vs 3,220)

## Root Cause Identified

The 17% archive bloat is **NOT caused by threading or race conditions** but by **fundamental algorithmic differences** in how RAGC groups segments compared to C++ AGC.

### Evidence
1. ✅ Single-threaded RAGC still has 161 extra groups (same bloat as multi-threaded)
2. ✅ Barriers and atomic allocation had minimal/zero effect
3. ✅ C++ AGC produces identical results single-threaded vs multi-threaded (barriers work)
4. ✅ RAGC's grouping algorithm creates more groups even without any concurrency

## Next Steps

Need to compare the core segment grouping/splitting logic between RAGC and C++ AGC:

1. **Segment splitting logic**: How do we decide when to split a segment?
2. **Middle k-mer selection**: How do we find the middle k-mer for splitting?
3. **Group creation heuristics**: When do we create a new group vs reuse existing?
4. **Segment classification**: How do we classify segments as "new" vs "known"?

The bug is in the **algorithmic logic**, not the **concurrency model**.
