# Known Issues

## ~~C++ AGC Compatibility Issue with Large Single-File Archives~~ [FIXED]

**Status:** FIXED (2025-10-25)
**Severity:** HIGH (caused data corruption with C++ AGC)
**Affects:** All multi-sample archives

### Symptoms (RESOLVED)

When compressing multi-sample FASTA files:
- C++ AGC reported "Corrupted archive!" for samples with in_group_id ≥ 5
- C++ AGC crashed with std::length_error when reading RAGC archives
- Small files (<10 samples) worked, but larger archives failed

### Root Cause

Two critical bugs were identified and fixed:

1. **Reference Segment Compression** (ragc-core/src/compressor_streaming.rs:3580):
   - RAGC was using `compress_segment_configured` with hardcoded marker byte 0 for reference segments
   - Should use `compress_reference_segment` which checks repetitiveness and uses tuple packing when beneficial
   - This caused C++ AGC to misinterpret compressed data

2. **Collection Metadata Tracking** (ragc-common/src/collection.rs:748-755):
   - Initial fix attempted to remove `&& seg.in_group_id > 0` condition to track references
   - This made RAGC incompatible with C++ AGC's deserialization logic
   - **Solution:** Reverted to match C++ AGC's exact behavior (including the condition)
   - The condition is preserved for compatibility even though it may seem incorrect

### Fixes Applied

```rust
// Fix 1: Reference segment compression (ragc-core/src/compressor_streaming.rs)
// BEFORE:
let compressed = compress_segment_configured(&ref_segment.data, self.config.compression_level)?;
compressed_with_marker.push(0); // Hardcoded marker

// AFTER:
let (compressed, marker) = compress_reference_segment(&ref_segment.data)?;
compressed_with_marker.push(marker); // Correct marker (0 = ZSTD, 1 = tuple-packed)

// Fix 2: Collection metadata tracking (ragc-common/src/collection.rs)
// Must match C++ AGC's exact condition (collection_v3.cpp:674-675):
if seg.in_group_id as i32 > prev_in_group_id && seg.in_group_id > 0 {
    self.set_in_group_id(seg.group_id as usize, seg.in_group_id as i32);
}
```

### Verification

After fixes:
- ✓ 2-sample archives work with C++ AGC
- ✓ 10-sample archives work with C++ AGC
- ✓ 235-sample split-file archives (146M) work with C++ AGC
- ✗ 235-sample single-file yeast235.fa.gz (87M) still crashes C++ AGC

**Note:** The single-file yeast235.fa.gz crash appears to be a separate, unrelated issue specific to that particular dataset's scale or structure.

### Related Files

- `ragc-core/src/compressor_streaming.rs:3580` - Reference segment compression fix
- `ragc-core/src/segment_compression.rs` - Compression/decompression with marker bytes
- `ragc-core/src/tuple_packing.rs` - Tuple packing implementation (matching C++ AGC)
- `ragc-common/src/collection.rs:748-755` - Collection metadata tracking (encoding)
- `ragc-common/src/collection.rs:907-912` - Collection metadata tracking (decoding)

### Last Updated

2025-10-25 - Fixes applied and verified

---

## Remaining Issue: yeast235.fa.gz Single-File Crash

**Status:** Under investigation
**Severity:** Medium (workaround available)
**Affects:** Large single-file multi-sample FASTA (235 samples, ~3GB)

### Symptoms

- ✓ RAGC can create and read the archive
- ✓ C++ AGC can list samples
- ✗ C++ AGC crashes (segfault) when extracting any sample
- ✓ Same samples in split-file format work fine with C++ AGC

### Workaround

Use split files instead of single-file format:
```bash
# Works: one file per sample
ragc create -o output.agc -k 21 -s 10000 -m 20 yeast_split/*.fa

# Crashes: all samples in one file
ragc create -o output.agc -k 21 -s 10000 -m 20 yeast235.fa.gz
```

### Next Steps

1. Investigate why C++ AGC crashes on large single-file archives
2. Compare segment/group structure between working and crashing archives
3. Check for buffer overflow or pointer issues in C++ AGC's deserialization
4. May be a C++ AGC bug rather than RAGC issue
