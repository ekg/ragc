# Segment Merging Bug in Multi-Sample Compression

## Status: PARTIALLY FIXED (Overlap), REMAINING ISSUE IDENTIFIED

### What We Fixed
- ✅ **K-mer overlap** (commit fb3ab43): Changed from k bytes to (k-1) bytes
  - Reduced data loss from 6.8% to ~2-3%
  - Single-sample archives now work perfectly (0% data loss)

### Remaining Issue: Segment Merging in Multi-Sample Archives

**Symptoms:**
- AAB#0 chrI alone: ~25 segments (correct)
- AAB#0 chrI with AAA#0 reference: ~17 segments (WRONG - segments merged!)
- Data loss: ~1.65-3% in subsequent samples

**Root Cause:**
In `streaming_compressor_queue.rs`, the `find_group_with_one_kmer()` function incorrectly matches segments to existing groups when:
1. AAB#0's segments have k-mers that don't exist in AAA#0's reference groups
2. The function searches for "connected" k-mers in the `map_segments_terminators`
3. It incorrectly matches AAB#0's segments to AAA#0's groups based on these connections
4. Multiple segments get assigned to the same group, causing them to be merged

**Evidence:**
```
AAB#0 chrI in multi-sample:
Segment 0: len=50,890 bytes (5x too large! Should be ~10,000)
Segment 2: len=19,368 bytes (2x too large)
Segment 9: len=19,999 bytes (2x too large)
```

**Test to Reproduce:**
```bash
# Single-sample: Perfect
./target/release/ragc create -o /tmp/aab_alone.agc -k 21 -s 10000 /tmp/check_AAB#0.fa
./target/release/ragc getset -o /tmp/extracted.fa /tmp/aab_alone.agc "AAB#0"
diff /tmp/check_AAB#0.fa /tmp/extracted.fa  # ✓ IDENTICAL

# Multi-sample: Bug!
./target/release/ragc create -o /tmp/multi.agc -k 21 -s 10000 /tmp/check_AAA#0.fa /tmp/check_AAB#0.fa
./target/release/ragc getset -o /tmp/extracted.fa /tmp/multi.agc "AAB#0"
diff /tmp/check_AAB#0.fa /tmp/extracted.fa  # ✗ DIFFERENT (~2% data loss)
```

**Location:**
- File: `ragc-core/src/streaming_compressor_queue.rs`
- Function: `find_group_with_one_kmer()` (line 758)
- Issue: Lines 807-836 - candidate group matching logic

**Fix Needed:**
The `find_group_with_one_kmer` function needs to be more conservative about matching segments to existing groups. It should only match when there's a direct k-mer correspondence, not through transitive connections in the terminators map.

Possible approaches:
1. Check that the matched group's reference segment actually shares sequence similarity
2. Create new groups more aggressively when k-mers don't match directly
3. Review C++ AGC's `find_cand_segment_with_one_splitter` logic for correct matching

**Impact:**
- Reference sample (first sample): ✓ Perfect
- Single-segment contigs (chrMT): ✓ Perfect in all samples
- Multi-segment contigs in subsequent samples: ✗ ~2-3% data loss

**Next Steps:**
1. Study C++ AGC's `find_cand_segment_with_one_splitter` implementation
2. Add debug output to track why segments are being matched to wrong groups
3. Implement stricter matching criteria in `find_group_with_one_kmer`
4. Add test case for multi-sample archives with divergent sequences
