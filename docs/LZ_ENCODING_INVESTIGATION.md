# LZ Encoding Investigation: RAGC vs C++ AGC

## Investigation Status: COMPLETE (2025-11-29)

## Executive Summary

**Root Cause Identified**: RAGC archives are ~7.8% larger than C++ AGC due to different reference selection order (`in_group_id` values), NOT LZ encoding differences.

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

### Attempted Fixes

1. **Deferred Reference Selection in flush_pack()**: Tried sorting segments before picking reference at flush time. Result: Made archives 18% LARGER because:
   - BufferedSegment stores already-transformed data
   - At flush time, original data isn't available for proper re-orientation
   - Would need to store both original and transformed data (2x memory)

## Conclusions

1. **LZ encoding is correct** - the algorithm matches C++ AGC exactly
2. **7.8% overhead is architectural** - streaming vs batch mode fundamentally differ in segment ordering
3. **Fix requires batch mode** - see `docs/BATCH_ARCHITECTURE_PLAN.md` for implementation plan

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
