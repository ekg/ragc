# Segment Split Threshold Bug - Multi-Sample Archive Bloat

**Date**: 2025-10-31
**Status**: Root cause identified, fix pending

## Problem

Multi-sample archives (multiple genome files) are **46% larger** than C++ AGC and have **data corruption** (1.2% data loss).

## Investigation

Created tiny test case: chr5 from 3 samples (one file per sample)

### Results
- **RAGC**: 220KB, 30 groups
- **C++ AGC**: 209KB, 24 groups (+5.36% bloat)

RAGC created 6 EXTRA groups (39-45) that should have been split and merged into existing groups.

## Root Cause

**RAGC's segment splitting threshold is TOO RESTRICTIVE**

**Location**: `ragc-core/src/compressor_streaming.rs:1459`

```rust
let can_split = kmer_front != MISSING_KMER
             && kmer_back != MISSING_KMER
             && segment_data.len() >= 2 * self.config.segment_size  // BUG: TOO HIGH!
             && !self.config.concatenated_genomes;
```

### Current Behavior
- RAGC: Requires segments >= **20,000 bytes** to split (2 × segment_size of 10,000)
- C++ AGC: Splits segments as small as **1,505 bytes**

### Evidence from C++ AGC Output

```
CPP_SPLIT_ATTEMPT 1: ... size=1505  ✓ SPLIT
CPP_SPLIT_ATTEMPT 2: ... size=20037 ✓ SPLIT
CPP_SPLIT_ATTEMPT 3: ... size=14053 ✓ SPLIT
CPP_SPLIT_ATTEMPT 4: ... size=20004 ✓ SPLIT
CPP_SPLIT_ATTEMPT 5: ... size=231232 ✓ SPLIT
```

C++ AGC splits segments regardless of size, as long as a "middle k-mer" connection exists.

## The 6 Extra RAGC Groups

These groups should have been SPLIT:

| RAGC Group | K-mer Pair | Should Split To |
|------------|------------|-----------------|
| 39 | MISSING -> 574464... | Group 17 (existing) |
| 40 | 22732249... -> 90843236... | Groups 28+27 |
| 41 | 57446420... -> 38130178... | Groups 18+19 |
| 42 | 84714973... -> 10182693... | Groups 26+27 |
| 43 | 97362197... -> 25398870... | Groups 33+34 |
| 45 | MISSING -> 90843236... | Group 26 or 27 |

## Impact

1. **Archive Size**: +5-46% larger (varies by dataset)
2. **Data Corruption**: 1.2% data loss (143KB in 11.6MB test)
   - Likely related: segments go to wrong groups when splitting fails
3. **Performance**: Extra groups = more overhead

## Fix Required

Change split threshold to match C++ AGC. Options:

1. **Remove size check entirely** - split whenever middle k-mer exists
2. **Lower threshold** - use segment_size instead of 2 × segment_size
3. **Match C++ AGC exactly** - investigate their actual threshold

## Test Case

```bash
# Tiny test case for validation
cd /tmp/chr5_test
ls -lh  # 3 files: AEL#2_chr5.fa, AIF#2_chr5.fa, ALI#2_chr5.fa

# Create archives
ragc create -o /tmp/chr5_ragc.agc -k 21 -s 10000 -m 20 *.fa
agc create -k 21 -s 10000 -l 20 -o /tmp/chr5_agc.agc *.fa

# Compare
grep "^GROUP_KMER:" /tmp/ragc_chr5_full.log | wc -l  # 30 groups
grep "^GROUP_KMER:" /tmp/cpp_groups.txt | wc -l      # 24 groups
```

## References

- RAGC split logic: `compressor_streaming.rs:1453-1526`
- C++ AGC: `src/agc-create.cpp` (find_cand_segment_with_missing_middle_splitter)
- Test script: `scripts/compare_agc_ragc.sh`
