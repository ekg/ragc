# Archive Size Difference Root Cause Analysis

## Executive Summary

RAGC archives are 2.07% larger than C++ AGC (61,911 bytes for yeast10 dataset).

**Root Cause**: RAGC's reference pack compression produces ~260 bytes larger output per group compared to C++ AGC, despite using the same ZSTD compression level (19).

## Detailed Findings

### Archive Structure Comparison

**RAGC (yeast10)**:
- Total size: 3,051,829 bytes
- Footer: 3,762 bytes (SMALLER ✓)
- Data: 3,048,059 bytes (LARGER ✗)
- Streams: 252 (245 ref + 3 collection + 4 other)
- Parts: 252 (1 per stream)

**C++ AGC (yeast10)**:
- Total size: 2,989,918 bytes  
- Footer: 5,642 bytes (larger due to 261 empty delta streams)
- Data: 2,984,268 bytes
- Streams: 513 (245 ref + 261 delta + 3 collection + 4 other)
- Parts: 268

**Key Insight**: C++ AGC creates empty delta streams (261 total, 0 parts each) for groups that only have a reference segment. RAGC skips these empty streams, saving footer space. However, RAGC's actual compressed data is larger.

### Per-Group Size Comparison

Sample ref streams (bytes):

| Group | RAGC | C++ AGC | Difference |
|-------|------|---------|------------|
| xHr   | 14,647 | 14,469 | +178 |
| xIr   | 15,166 | 14,797 | +369 |
| xJr   | 15,063 | 14,696 | +367 |
| xKr   | 12,215 | 12,168 | +47 |
| xNr   | 14,855 | 14,621 | +234 |
| xOr   | 15,138 | 14,761 | +377 |

**Average overhead**: ~260 bytes per group
**Total for 245 groups**: 63,791 bytes

### Groups 0-15 Mystery

C++ AGC creates groups 0-15 (x0d through xFd) with 2 bytes each. RAGC doesn't create these groups at all. However, this accounts for only 32 bytes of the difference, not the 63KB.

### Compression Code Analysis

**RAGC** (`compressor_streaming.rs:263-273`):
```rust
let compressed = compress_segment_configured(&segment.data, config.compression_level)?;
let (compressed_data, uncompressed_size) =
    if compressed.len() + 1 < segment.data.len() {
        let mut compressed_with_marker = compressed;
        compressed_with_marker.push(0); // Marker byte
        (compressed_with_marker, segment.data.len() as u64)
    } else {
        (segment.data.clone(), 0)  // Uncompressed
    };
```

**Compression level**: 19 (confirmed in `segment_compression.rs:19`)

### Hypothesis

RAGC's compressed output is ~260 bytes larger per group despite:
1. Same compression level (19)
2. Same input data (segments are identical)
3. Same ZSTD algorithm

**Possible causes**:
1. ❓ Marker byte affecting compression (only 1 byte though)
2. ❓ Different ZSTD version or settings
3. ❓ Extra data being included in what gets compressed
4. ❓ Alignment or padding differences
5. ❓ Context/dictionary differences

### Next Steps

1. Compare RAGC vs C++ AGC's actual bytes being compressed for one segment
2. Check if ZSTD context/dictionary settings differ
3. Check if there's metadata being compressed along with data
4. Verify ZSTD versions match

## Files Modified

None yet - investigation phase.

## References

- Investigation started: 2025-10-30
- Test dataset: yeast10 (10 contigs, 245 groups)
- RAGC code: `ragc-core/src/compressor_streaming.rs:258-293`
- C++ AGC code: TBD
