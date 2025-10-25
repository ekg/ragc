# Tuple Packing Implementation Notes

**Date**: 2025-10-25
**Status**: Module implemented but NOT integrated (requires repetitiveness check)

## Summary

Implemented complete tuple packing/unpacking matching C++ AGC format in `ragc-core/src/tuple_packing.rs`. All unit tests pass. However, **integration into compression pipeline requires repetitiveness check** that is not yet implemented.

## Why Not Integrated

C++ AGC does NOT always use tuple packing! It chooses based on repetitiveness:

```cpp
// segment.h:218-255
double best_frac = 0.0;
for (uint32_t i = 4; i < 32; ++i) {
    cnt += data[j] == data[(size_t) j + i];
    frac = (double)cnt / cur_size;
    if (frac > best_frac) best_frac = frac;
}

if (best_frac < 0.5)
    add_to_archive_tuples(...);  // Marker byte 1
else
    add_to_archive(...);         // Marker byte 0
```

**Translation**: If repetitiveness < 50%, use tuple packing. Otherwise use plain ZSTD.

## Current Status

- ✅ `tuple_packing.rs` - Complete and tested
- ❌ Repetitiveness check - Not implemented
- ❌ Integration - Cannot integrate without repetitiveness check

## What Happens If We Integrate Without Repetitiveness Check

Tested this - it breaks everything:
- Applied tuple packing to ALL reference segments
- Reference segments that C++ AGC expects plain ZSTD get tuple-packed data
- Both C++ AGC and RAGC fail to decompress

## Proper Integration Steps

1. Implement repetitiveness check (C++ AGC segment.h:224-249)
2. Add `compress_segment_auto()` that:
   - Checks repetitiveness
   - Chooses tuple packing if < 50%
   - Chooses plain ZSTD if >= 50%
   - Returns compressed data + marker byte (0 or 1)
3. Update reference segment compression to use `compress_segment_auto()`
4. Keep delta packs using plain ZSTD always

## Files

- `ragc-core/src/tuple_packing.rs` - Implemented, tested, ready to use
- `ragc-core/src/segment_compression.rs` - Currently uses plain ZSTD only
- `ragc-core/src/compressor_streaming.rs` - Reference segment compression (line ~3580)

## Recommendation

Do NOT integrate tuple packing until repetitiveness check is implemented. The current plain ZSTD approach is correct for delta packs and safe (if suboptimal) for reference segments.

The C++ AGC compatibility bug is NOT related to tuple packing - it's something else in the collection metadata or archive structure.
