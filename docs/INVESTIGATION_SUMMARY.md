# Segment Grouping Investigation Summary

**Date**: 2025-10-30
**Problem**: Same splitters (11,771), different archive sizes (C++ AGC: 3.5 MB, RAGC: 3.9 MB = +11%)

---

## Findings So Far

### ✅ 1. Splitters Match Exactly
- Both implementations find **11,771 splitters** from reference genome
- Verified by running both with same parameters

### ✅ 2. Segmentation Algorithm Matches
- Created test comparing C++ AGC's `compress_contig` with RAGC's `split_at_splitters_with_size`
- Test Result: **Identical segment counts** on test contigs
- Both split at every splitter occurrence, create k-base overlaps, reset k-mers after splits

### ⚠️ 3. Group Counts Differ
- C++ AGC output mentions various segment counts
- RAGC reports creating 1000+ group_ids, but final stats show bugs
- Archive sizes differ by +11% (3.9 MB vs 3.5 MB)

### 🔍 4. Potential Grouping Logic Difference Found

**C++ AGC** (agc_compressor.cpp:1302-1361):
- Both k-mers present: Normalize to `(min, max)`, RC if needed
- One k-mer present: **Complex logic** - calls `find_cand_segment_with_one_splitter()`
  - Searches for existing groups with matching k-mer
  - Compares segment sequences
  - May assign different group than simple (kmer, MISSING)
- No k-mers present: Uses fallback minimizers

**RAGC** (compressor_streaming.rs:167-199):
- Both k-mers present: Normalize to `(min, max)`, RC if needed
- One k-mer present: **Simple logic** - just uses `(kmer, MISSING_KMER)`
  - No search for existing groups
  - No sequence comparison
  - Creates new group directly

**Hypothesis**: Segments with one k-mer are being grouped differently!

---

## What We Need To Compare

To prove this hypothesis, we need to dump from both archives:

1. **Group metadata**:
   - Group ID
   - Front k-mer value (or MISSING)
   - Back k-mer value (or MISSING)
   - Number of segments in this group

2. **Segment metadata** for each group:
   - Sample name
   - Contig name
   - Segment position in contig
   - Segment length
   - First 20 bytes of sequence (for verification)

3. **Statistics**:
   - How many groups have both k-mers?
   - How many groups have one k-mer?
   - How many groups have no k-mers?
   - Average segments per group

---

## Key Questions

1. **Do one-k-mer segments create different groups?**
   - C++ AGC may merge them with existing groups
   - RAGC creates new groups for each unique (kmer, MISSING) pair

2. **Are no-k-mer segments handled differently?**
   - C++ AGC uses fallback minimizers
   - RAGC may use (MISSING, MISSING) always?

3. **Is sequence deduplication happening in C++ AGC?**
   - When finding "best match" for one-k-mer segments
   - May group identical sequences under different k-mer pairs

---

## Next Steps

1. **Create dump tool** for RAGC archives
   - Read archive structure
   - Extract group metadata
   - Dump to JSON/text format

2. **Create equivalent tool** for C++ AGC
   - May need to modify C++ AGC source to add dump command
   - Or read archive format directly

3. **Compare dumps**:
   - Sort by (front_kmer, back_kmer)
   - Find groups that exist in one but not the other
   - Analyze the segments in those groups

4. **Identify specific bug**:
   - Is it the one-k-mer handling?
   - Is it the no-k-mer handling?
   - Is it something else?

5. **Fix RAGC** to match C++ AGC exactly

---

## Immediate Action

The fastest path forward:
1. Modify RAGC to log ALL group creations with k-mer pairs
2. Modify C++ AGC (or read source) to see what groups it creates
3. Compare the lists directly
4. Find the first divergence

This will pinpoint exactly where the logic differs.
