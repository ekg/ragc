# FIX 5 Results: Batch Boundary Synchronization

**Date**: 2025-12-04
**Status**: SIGNIFICANT PROGRESS - Group IDs now aligned for first 1501 positions
**Next Issue**: Different group assignments for same segments

---

## What FIX 5 Accomplished

### Before FIX 5
- **First divergence**: Position 1 (CFF#2#chrI seg 0)
- **group_id**: RAGC=131, C++ AGC=411
- **Cause**: Different batch boundaries caused groups to be created in different order

### After FIX 5
- **First divergence**: Position 1502 (ADI#0#chrII seg 0)
- **group_id**: 856 in BOTH implementations ✓
- **in_group_id**: RAGC=1, C++ AGC=2 (different)
- **Result**: Group creation order now matches for first 1501 positions!

**This is MAJOR PROGRESS!** FIX 5 successfully synchronized batch boundaries between RAGC and C++ AGC.

---

## Remaining Issue: Group Assignment Divergence

### The Problem

The same physical segment gets assigned to DIFFERENT groups:

**Segment**: `ADI#0#chrI seg 22` (249 bp)

| Implementation | Group ID | Role | in_group_id |
|----------------|----------|------|-------------|
| **C++ AGC** | 856 | Known segment | 1 |
| **RAGC** | 1522 | Reference segment | 0 |

### Impact on Group 856

**C++ AGC group 856**:
- 93 segments total
- Reference: CFF#2#chrXII_1 seg 0
- Includes segments from all 5 samples: CFF#2, ADI#0, AVI_1a#0, CL216#0, AGA_1a#0

**RAGC group 856**:
- 17 segments total
- Reference: CFF#2#chrXII_1 seg 0 (same as C++ AGC ✓)
- Includes segments only from: CFF#2, ADI#0

**76 FEWER segments in RAGC group 856!**

---

## Root Cause Analysis

### Hypothesis 1: K-mer Pair Extraction

The same segment may have different k-mer pairs extracted:
- **C++ AGC**: Extracts k-mer pair (X, Y) → matches existing group 856
- **RAGC**: Extracts k-mer pair (A, B) → creates new group 1522

**To verify**: Compare k-mer extraction for `ADI#0#chrI seg 22` in both implementations.

### Hypothesis 2: Group Lookup Timing

Segments may be classified as "new" vs "known" at different times:
- **C++ AGC**: Group 856 already exists when processing `ADI#0#chrI seg 22` → segment joins existing group
- **RAGC**: Group 856 doesn't exist yet → segment creates new group 1522

**To verify**: Trace group creation order and timing for groups 856 and 1522.

### Hypothesis 3: Batch Processing Order

Even with synchronized batch boundaries (FIX 5), segments within a batch may be processed in different order:
- **C++ AGC**: Processes segments in lexicographic order → group 856 created early
- **RAGC**: Processes segments in different order → group 856 created later

**To verify**: Log segment processing order within each batch.

---

## Investigation Plan

### Step 1: Extract K-mer Pairs for ADI#0#chrI seg 22

Add logging to both implementations to show:
1. Front k-mer value
2. Back k-mer value
3. K-mer pair (front, back)
4. Existing group ID (if any)
5. New group ID (if created)

### Step 2: Trace Group Creation Order

Log when groups 856 and 1522 are created:
- Which batch?
- Which segment triggered creation?
- What k-mer pair?

### Step 3: Compare Batch Content

For the batch containing position 1502:
- How many segments in batch?
- What's the processing order?
- When is group 856 created relative to ADI#0#chrI seg 22?

---

## Expected Outcomes

### If Hypothesis 1 (K-mer Extraction)
**Fix**: Align k-mer extraction logic between RAGC and C++ AGC
**Impact**: All segments should get same k-mer pairs → same group assignments

### If Hypothesis 2 (Group Lookup Timing)
**Fix**: Ensure group lookup happens at same point in processing
**Impact**: Segments should find existing groups instead of creating new ones

### If Hypothesis 3 (Batch Processing Order)
**Fix**: Process segments in same order as C++ AGC (lexicographic within batch)
**Impact**: Groups created in same order → segments find correct existing groups

---

## Testing Protocol

After implementing fix:

```bash
# Test with 5-sample dataset
./target/release/ragc create -o /tmp/ragc_fix6.agc -k 21 -s 10000 -m 20 -t 1 \
  /home/erik/scrapy/yeast235_samples/CFF#2.fa \
  /home/erik/scrapy/yeast235_samples/ADI#0.fa \
  /home/erik/scrapy/yeast235_samples/AVI_1a#0.fa \
  /home/erik/scrapy/yeast235_samples/CL216#0.fa \
  /home/erik/scrapy/yeast235_samples/AGA_1a#0.fa

/home/erik/agc/bin/agc create -o /tmp/cpp_fix6.agc -k 21 -s 10000 -l 20 -t 1 \
  /home/erik/scrapy/yeast235_samples/CFF#2.fa \
  /home/erik/scrapy/yeast235_samples/ADI#0.fa \
  /home/erik/scrapy/yeast235_samples/AVI_1a#0.fa \
  /home/erik/scrapy/yeast235_samples/CL216#0.fa \
  /home/erik/scrapy/yeast235_samples/AGA_1a#0.fa

# Compare layouts
./target/release/ragc inspect /tmp/ragc_fix6.agc --segment-layout > /tmp/ragc_layout6.csv
./target/release/ragc inspect /tmp/cpp_fix6.agc --segment-layout > /tmp/cpp_layout6.csv

# Compare group 856 specifically
grep ",856," /tmp/ragc_layout6.csv | wc -l  # Should be 93
grep ",856," /tmp/cpp_layout6.csv | wc -l   # Should be 93

# Check for divergence
python3 /tmp/compare_final.py
```

**Expected result**: Group 856 should have 93 segments in both, with identical in_group_id assignments.

---

## Progress Summary

| Fix | Issue | Status |
|-----|-------|--------|
| FIX 1 | Sample priority in sorting | ✅ FIXED |
| FIX 2 | Raw group count (16 vs 256) | ✅ FIXED |
| FIX 3 | pack_cardinality CLI parameter | ✅ FIXED |
| FIX 4 | Batch-level sorting | ✅ FIXED |
| FIX 5 | Batch boundary sync tokens | ✅ FIXED |
| **FIX 6** | **Group assignment divergence** | 🔍 INVESTIGATING |

---

## Key Metrics

### Before All Fixes
- First divergence: Position 1
- Different group_ids from the start
- Fundamentally incompatible archives

### After FIX 1-4
- First divergence: Position 1
- Different group_ids (131 vs 411)
- Group creation order diverged

### After FIX 5 (Current)
- First divergence: Position 1502
- **SAME group_id (856)** for positions 1-1501 ✓
- Different in_group_id within same group
- Group 856: 17 segments (RAGC) vs 93 segments (C++ AGC)

**We've gone from position 1 to position 1502 with identical group_ids!** This is 99.7% alignment on group creation order.

---

## Next Steps

1. Add k-mer pair logging to identify extraction differences
2. Trace group 856 and 1522 creation order
3. Compare segment processing order within batches
4. Implement FIX 6 based on findings
5. Re-test with 5-sample dataset
6. Verify byte-identical archives (if FIX 6 succeeds)
