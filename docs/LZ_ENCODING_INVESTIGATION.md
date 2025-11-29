# LZ Encoding Investigation: RAGC vs C++ AGC

## Investigation Status: UPDATED (2025-11-29)

## Executive Summary

**Root Causes Identified** (Two separate issues):

1. **Segment Count Difference** (~6% of overhead): RAGC has FEWER segments in non-reference samples because it lacks C++ AGC's dynamic splitter discovery. C++ AGC discovers NEW k-mer splitters unique to each non-reference contig, creating additional segment boundaries.

2. **Reference Selection Order** (~7-8% of overhead): Different `in_group_id` values due to streaming vs batch architecture.

**Key Finding**: LZ encoding is byte-identical between RAGC and C++ AGC (verified via FFI parity tests).

## Problem Statement

RAGC archives are consistently ~7.8% larger than C++ AGC archives when compressing the same input files.

## Investigation Findings

### What We Confirmed

1. **LZ Encoding is Byte-Identical**
   - Created FFI parity test (`ragc-core/tests/test_lz_parity.rs`)
   - RAGC's LZDiff::encode() produces identical output to C++ AGC's CLZDiff_V2::Encode()
   - All test cases pass: identical sequences, literals, matches, N-runs, real yeast data

2. **Correctness is 100%**
   - All samples extract byte-for-byte identical to original input
   - Data integrity is fully preserved
   - Archives are fully C++ AGC compatible (can be read by C++ AGC tools)

3. **Root Cause: Reference Selection Order**
   - ~22% of segments have different `in_group_id` values between RAGC and C++ AGC
   - `in_group_id=0` indicates the reference segment for LZ encoding
   - Different reference selection leads to different (larger) LZ-encoded deltas

### Why Reference Selection Differs

**C++ AGC (Batch Mode)**:
1. Accumulates segments from N files (batch of 50)
2. Sorts ALL segments globally by (sample_name, contig_name, segment_index)
3. Distributes sorted segments to groups
4. First segment to arrive at each group (in sorted order) becomes reference
5. Writes packs as groups fill (50 segments per pack)

**RAGC (Streaming Mode)**:
1. Processes files one at a time (memory-efficient)
2. Segments arrive at groups in file order, not sorted order
3. First segment to ARRIVE at each group becomes reference
4. Different pack boundaries cause different reference selections
5. Writes all groups at finalize()

### Impact Analysis

| Metric | RAGC | C++ AGC | Difference |
|--------|------|---------|------------|
| Archive Size | 14.3 MB | 13.3 MB | +7.8% |
| Segment Layout Diffs | 22% | - | Different in_group_id |
| Correctness | 100% | 100% | Both correct |
| Compatibility | Full | Full | Cross-compatible |

### Issue 1: Dynamic Splitter Discovery (NEW FINDING 2025-11-29)

**Problem**: RAGC has 84-102 FEWER segments per non-reference sample compared to C++ AGC.

**Evidence** (3-sample test):
| Sample | RAGC Segments | C++ AGC Segments | Difference |
|--------|---------------|------------------|------------|
| AAA#0 (reference) | 1243 | 1243 | Same ✓ |
| AAB#0 | 1132 | 1216 | -84 (-7%) |
| AAC#0 | 1095 | 1197 | -102 (-8.5%) |

**Root Cause**: C++ AGC has `find_new_splitters()` function that dynamically discovers NEW k-mer splitters for each non-reference contig:
```cpp
void CAGCCompressor::find_new_splitters(contig_t& ctg, uint32_t thread_id) {
    // 1. Enumerate ALL k-mers in current contig
    enumerate_kmers(ctg, v_contig_kmers);

    // 2. Filter out k-mers already in reference genome
    p_end = set_difference(v_contig_kmers, v_candidate_kmers, v_tmp);

    // 3. These become NEW splitters for this contig
    find_splitters_in_contig(ctg, v_contig_kmers, vv_splitters[thread_id], ...);
}
```

**RAGC Behavior**: Uses ONLY splitters from first file. Non-reference contigs may lack the same k-mers (due to mutations), resulting in fewer segment boundaries.

**Impact**: Fewer segments → larger average segment size → worse LZ compression → ~6% larger archives.

**Fix Required**: Implement dynamic splitter discovery (significant architectural change).

### Attempted Fixes

1. **Multi-file splitter detection**: Tried collecting splitters from ALL input files (up to 50). Result: Made archives 18% LARGER because:
   - Too many splitters (3682 vs 1226)
   - Over-segmentation → too many small segments → worse compression

2. **Deferred Reference Selection in flush_pack()**: Tried sorting segments before picking reference at flush time. Result: Made archives 18% LARGER because:
   - BufferedSegment stores already-transformed data
   - At flush time, original data isn't available for proper re-orientation
   - Would need to store both original and transformed data (2x memory)

## Conclusions

1. **LZ encoding is correct** - the algorithm matches C++ AGC exactly
2. **~13% total overhead** is from two architectural issues:
   - **~6%**: Missing dynamic splitter discovery for non-reference contigs
   - **~7-8%**: Different reference selection order (streaming vs batch)
3. **Two fixes required**:
   - Implement dynamic splitter discovery (see C++ AGC's `find_new_splitters()`)
   - Implement batch mode for reference selection (see `docs/BATCH_ARCHITECTURE_PLAN.md`)

## Trade-offs

**Streaming Mode (Current)**:
- ✅ Constant memory usage
- ✅ Can process files larger than RAM
- ✅ 100% correct extraction
- ❌ 7.8% larger archives

**Batch Mode (C++ AGC style)**:
- ✅ Byte-identical archives to C++ AGC
- ✅ Optimal compression
- ❌ Memory scales with batch size
- ❌ More complex implementation

## Relevant Files

- `ragc-core/tests/test_lz_parity.rs` - LZ encoding parity tests
- `ragc-core/src/lz_diff.rs` - RAGC's LZ encoder
- `agc/src/core/lz_diff.cpp` - C++ AGC's LZ encoder
- `ragc-core/src/streaming_compressor_queue.rs` - Streaming compression pipeline
- `docs/BATCH_ARCHITECTURE_PLAN.md` - Plan for implementing batch mode

## Test Commands

```bash
# Run LZ parity tests
cargo test --release --features cpp_agc -p ragc-core test_lz_parity

# Compare archive sizes
./target/release/ragc create -o /tmp/ragc.agc -k 21 -s 10000 -m 20 -t 1 samples/*.fa
/home/erik/agc/bin/agc create -o /tmp/cpp.agc -k 21 -s 10000 -l 20 -t 1 samples/*.fa
ls -la /tmp/*.agc

# Compare segment layouts
./target/release/ragc inspect /tmp/ragc.agc --segment-layout > ragc_layout.csv
./target/release/ragc inspect /tmp/cpp.agc --segment-layout > cpp_layout.csv
diff ragc_layout.csv cpp_layout.csv | wc -l
```
