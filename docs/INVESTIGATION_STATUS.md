# C++ AGC Compatibility Investigation - Status Report

**Date**: 2025-10-25
**Issue**: C++ AGC cannot read RAGC archives when they contain 8+ delta-encoded samples per group
**Status**: Root cause not yet identified, investigation ongoing

## Summary

RAGC archives work perfectly when read by RAGC itself, but fail when C++ AGC tries to extract samples 8+ (in_group_id ≥ 7) from RAGC-created archives.

## What Works

- ✅ RAGC can read its own archives (all samples extract correctly)
- ✅ Small archives (< 8 delta samples): C++ AGC extracts correctly
- ✅ Samples 1-7 in large archives: C++ AGC extracts correctly
- ✅ 2-sample archives: Both RAGC and C++ AGC work perfectly

## What Fails

- ❌ Sample 8+ in archives with 8+ delta-encoded samples
- ❌ C++ AGC error: "Corrupted archive!" or `std::length_error`
- ❌ Consistent failure at in_group_id=7 (8th sample in a group)

## Investigation Attempts

### 1. Tuple Packing Hypothesis ❌

**Theory**: C++ AGC uses tuple packing with threshold at value 6 (me < 6 uses 3/byte, me >= 6 uses 2/byte)

**Finding**: NOT the root cause
- C++ AGC only uses tuple packing for SOME reference segments based on repetitiveness
- Delta packs ALWAYS use plain ZSTD (segment.h:279)
- Collection metadata uses variable-width integer encoding, NOT tuple packing
- RAGC already uses plain ZSTD correctly

**Implementation**: Created complete tuple packing module (`ragc-core/src/tuple_packing.rs`) matching C++ AGC
- All unit tests pass
- Roundtrip compression/decompression works
- But applying tuple packing to all segments made corruption WORSE (all samples failed)

### 2. Segment Compression Format

**C++ AGC** (segment.h:172-215):
- `add_to_archive()`: Plain ZSTD, marker byte 0
- `add_to_archive_tuples()`: Tuple-packed ZSTD, marker byte 1
- Delta packs: Always use `add_to_archive()` with marker 0
- Reference segments: Choose based on repetitiveness

**RAGC** (`compressor_streaming.rs`):
- Uses plain ZSTD for all segments
- Marker byte 0
- Matches C++ AGC delta pack format ✅

### 3. Collection Metadata Encoding

**C++ AGC** (`collection_v3.cpp`):
- Uses `append()` function for variable-width integer encoding
- in_group_id, seg_part_no, is_rev_comp encoded as variable-width integers
- Thresholds: < 64, < 16384, < 4194304, >= 4194304

**RAGC** (`ragc-common/src/collection.rs`):
- Uses `CollectionVarInt::encode()` for variable-width integers
- Should match C++ AGC format
- Needs byte-by-byte verification ⏳

## Tested Scenarios

| Archive | Samples | Sample 1 | Sample 7 | Sample 8 | Sample 10 |
|---------|---------|----------|----------|----------|-----------|
| 2 samples | 2 | ✅ | N/A | N/A | N/A |
| 8 samples | 8 | ✅ | ❌ | ❌ | N/A |
| 10 samples | 10 | ✅ | ✅ | ❌ | ❌ |

**Observation**: 8-sample archive fails at sample 7, but 10-sample archive works at sample 7!
This suggests the bug depends on total sample count, not just in_group_id value.

## Current Hypothesis

The bug is likely in one of these areas:

1. **Collection metadata encoding**: Variable-width integer encoding may have subtle differences
2. **Group boundary handling**: How groups are split/merged when sample counts change
3. **Archive stream ordering**: Order of writing collection-details vs segment streams
4. **Metadata size fields**: Compressed vs uncompressed size handling

## Recommended Next Steps

### Immediate Actions

1. **Create minimal reproducible case**:
   ```bash
   # 2 tiny samples that trigger the bug
   head -20 sample1.fa > tiny1.fa
   head -20 sample2.fa > tiny2.fa
   # Create archives with both implementations
   # Compare byte-by-byte
   ```

2. **Compare collection metadata**:
   - Extract collection-details stream from both archives
   - Hex dump and compare byte-by-byte
   - Verify variable-width integer encoding matches exactly

3. **Inspect archive structure**:
   - Use `agc info` on both archives
   - Compare stream names, sizes, ordering
   - Check for missing or extra streams

4. **Debug with C++ AGC source**:
   - Add debug prints to C++ AGC decompressor
   - Find exactly where "Corrupted archive!" is triggered
   - Check what value causes the corruption

### Long-term Solutions

1. **Format specification**: Document exact AGC archive format
2. **Test suite**: Add compatibility tests that verify byte-identical output
3. **Fuzzing**: Generate random archives and cross-test with C++ AGC

## Files Modified in Investigation

- `ragc-core/src/tuple_packing.rs` - Created (for potential future use)
- `ragc-core/src/segment_compression.rs` - Reverted to plain ZSTD
- `ragc-core/src/lib.rs` - Added tuple_packing module
- `ragc-core/src/compressor_streaming.rs` - No changes (already correct)
- `docs/TUPLE_PACKING_BUG.md` - Created with findings
- `docs/INVESTIGATION_STATUS.md` - This file

## Test Artifacts

- `/tmp/test_10samples.agc` - 10-sample RAGC archive (samples 8-10 fail)
- `/tmp/test_8samples.agc` - 8-sample RAGC archive (samples 7-8 fail)
- `/tmp/ragc_2samples.agc` - 2-sample RAGC archive (all work)
- `/tmp/cpp_2samples.agc` - 2-sample C++ AGC archive (for comparison)

## Key Code References

### C++ AGC
- `/home/erik/agc/src/common/segment.h:172-215` - Segment compression
- `/home/erik/agc/src/common/segment.cpp:254-272` - Segment decompression
- `/home/erik/agc/src/common/collection_v3.cpp` - Collection metadata
- `/home/erik/agc/src/common/collection.h:126-160` - Variable-width encoding

### RAGC
- `ragc-core/src/segment_compression.rs` - Segment compression/decompression
- `ragc-core/src/compressor_streaming.rs` - Main compression pipeline
- `ragc-common/src/collection.rs` - Collection metadata encoding
- `ragc-common/src/types.rs` - Type definitions

## Conclusion

The bug is real and reproducible, but the root cause remains elusive. Tuple packing was a red herring. The issue likely lies in collection metadata encoding or archive structure, but requires byte-level comparison to identify.

The fact that RAGC can read its own archives perfectly suggests the format is internally consistent, but has subtle incompatibility with C++ AGC's expectations at sample count thresholds.
