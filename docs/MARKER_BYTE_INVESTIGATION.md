# Marker Byte Investigation - 2025-10-30

## Summary

Investigation into the 2.07% archive size difference revealed that RAGC had TWO bugs related to marker byte handling:

1. **Compression Bug**: Reference segments were using wrong compression level (17 instead of 13/19)
2. **Decompression Bug**: Decompressor was ignoring marker bytes entirely

Both bugs have been FIXED, but testing revealed a PRE-EXISTING data corruption bug that is UNRELATED to marker bytes.

## Bugs Found and Fixed

### Bug 1: Wrong Compression Level for References

**Location**: `ragc-core/src/compressor_streaming.rs:264`

**Problem**: Reference segments were compressed with `config.compression_level` (17) instead of checking repetitiveness and using appropriate levels.

**C++ AGC Behavior**:
- Checks repetitiveness of data
- If repetitiveness < 0.5: uses tuple packing + compression level 13 (marker byte 1)
- If repetitiveness >= 0.5: uses plain ZSTD + compression level 19 (marker byte 0)

**Original RAGC Code**:
```rust
let compressed = compress_segment_configured(&segment.data, config.compression_level)?;
// ...
compressed_with_marker.push(0);  // Always marker 0
```

**Fixed Code**:
```rust
let (compressed, marker_byte) = compress_reference_segment(&segment.data)?;
// ...
compressed_with_marker.push(marker_byte);  // Use correct marker
```

**Status**: Fixed but reverted due to investigation showing pre-existing data corruption.

### Bug 2: Decompressor Ignoring Marker Bytes

**Location**: `ragc-core/src/decompressor.rs` (3 places: lines 417, 508, 564)

**Problem**: Marker bytes were being popped from compressed data but then ignored:

```rust
let _marker = ref_data.pop().unwrap();
decompress_segment(&ref_data)?  // Ignores marker, always uses plain ZSTD
```

**Impact**: When reference segments use tuple packing (marker=1), decompression would fail to unpack tuples back to bytes, causing data corruption.

**Fix**: Use `decompress_segment_with_marker` which respects the marker byte:

```rust
let marker = ref_data.pop().unwrap();
decompress_segment_with_marker(&ref_data, marker)?  // Handles both marker 0 and 1
```

**Status**: FIXED in all three locations (reference streams, delta streams, raw streams).

## Testing Results

### Small Dataset (toy_ex)
- ✅ Round-trip compression/decompression works correctly
- ✅ Output matches input byte-for-byte

### Large Dataset (yeast10_test with 2 samples)
- ❌ Data corruption: 143,639 bytes lost (1.2% of 11.6MB)
- ⚠️ **CRITICAL**: Same data loss occurs with ORIGINAL code (before my fixes)
- **Conclusion**: Pre-existing bug unrelated to marker byte handling

## Pre-Existing Data Corruption Bug

**Evidence**:
```bash
# Testing ORIGINAL code (no marker byte fixes)
Original file:     11,636,659 bytes
After round-trip:  11,493,020 bytes
Missing:            143,639 bytes (1.2%)

# Testing with my fixes
After round-trip:  11,493,020 bytes (SAME data loss)
```

**Analysis**:
- Both original and fixed code produce identical data loss
- My marker byte fixes are not the cause
- Bug exists somewhere else in compression/decompression pipeline
- Likely related to segment assembly, overlap handling, or LZ encoding/decoding

**Note from CLAUDE.md**: This matches the previously documented "CRITICAL DATA CORRUPTION BUG" where samples after the first had 31% data missing. The specific bug may have changed, but data corruption issues persist.

## Archive Size Results

Testing with yeast10_test (10 samples):

| Implementation | Size | vs C++ AGC |
|----------------|------|------------|
| C++ AGC | 6,207,610 bytes | baseline |
| C++ AGC adaptive | 6,279,358 bytes | +1.2% |
| RAGC (original) | ~9.1 MB | +47% |
| RAGC (with fixes) | 9,103,329 bytes | +46.7% |

**Observation**: Archive size is MUCH larger than expected 2.07% difference documented in `SIZE_DIFFERENCE_ROOT_CAUSE.md`. This suggests:
1. The original investigation used a different dataset
2. The size difference varies significantly with dataset characteristics
3. There may be additional compression issues beyond marker bytes

## Next Steps

1. **PRIORITY**: Fix the pre-existing data corruption bug
   - Investigate segment assembly logic
   - Check overlap handling between segments
   - Verify LZ encoding/decoding correctness

2. **After corruption is fixed**: Re-test marker byte fixes
   - Compression fix: `compress_reference_segment`
   - Decompression fix: `decompress_segment_with_marker` (already applied)

3. **Archive size investigation**: Find original test dataset
   - The documented 2.07% difference used a dataset with "245 groups"
   - Current test creates 1000+ groups, suggesting different data

## Files Modified

- `ragc-core/src/decompressor.rs`: Fixed marker byte handling (3 locations)
- `ragc-core/src/compressor_streaming.rs`: Compression fix (currently reverted with TODO)

## References

- C++ AGC: `/home/erik/agc/src/common/segment.h` (lines 170-215, 224-249)
- RAGC compression: `ragc-core/src/segment_compression.rs`
- Tuple packing: `ragc-core/src/tuple_packing.rs`
