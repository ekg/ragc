# Streaming Queue API - Missing Compression

**Date**: 2025-11-04
**Status**: ⚠️ **CORRECTNESS VERIFIED, COMPRESSION MISSING**

## Current State

The streaming queue API has been fully validated for correctness:
- ✅ Archive format correct
- ✅ C++ AGC compatibility verified
- ✅ Byte-for-byte identical extraction
- ✅ Multi-sample support working
- ✅ Constant memory usage achieved

**However**: The actual compression is not implemented yet!

## The Issue

### Archive Sizes (yeast235 dataset)
- Uncompressed input: **3.3 GB**
- Batch mode (with compression): **88 MB** (37x compression)
- Streaming mode (no compression): **3.2 GB** (barely compressed)

### Root Cause

In `streaming_compressor_queue.rs:488-502`, we write segments as raw uncompressed data:

```rust
// Write raw uncompressed segment data (metadata=0)
// This matches batch compressor format and works with C++ AGC
arch.add_part(stream_id, &segment.data, 0)?;
```

This was the fix for C++ AGC compatibility (metadata=0), but we're not compressing the data before writing it!

## What's Missing

### 1. ZSTD Compression (Required)
At minimum, we need to ZSTD-compress each segment before writing:

```rust
// Compress segment with ZSTD
let compressed = compress_segment_data(&segment.data)?;
arch.add_part(stream_id, &compressed, 0)?;  // Still metadata=0 for C++ AGC
```

### 2. Segment Grouping (Optional, but improves compression)
Batch mode groups segments by (kmer1, kmer2) pairs to reduce metadata overhead. Streaming mode currently treats each segment as its own group.

### 3. LZ Differential Encoding (Future optimization)
True AGC-style differential encoding between reference and delta segments. This is marked as TODO in both batch and streaming.

## Comparison with Batch Mode

**Batch mode (`worker.rs:78-100`)**:
1. Groups segments by (kmer1, kmer2)
2. Compresses each segment with ZSTD
3. Writes compressed data
4. **Result**: 88 MB for yeast235 (37x compression)

**Streaming mode (current)**:
1. Each segment is own group
2. ~~No compression~~ ❌
3. Writes raw data
4. **Result**: 3.2 GB for yeast235 (barely compressed)

## Next Steps

### Priority 1: Add ZSTD Compression (Required)
- Compress segment data before writing
- Should bring streaming mode to ~88 MB (similar to batch)
- Estimated effort: 1-2 hours

### Priority 2: Add Segment Grouping (Optional)
- Group segments by (kmer1, kmer2) like batch mode
- Reduces metadata overhead
- Estimated effort: 3-4 hours

### Priority 3: LZ Differential Encoding (Future)
- Compute byte-level diffs between reference/delta segments
- True AGC compression
- Would benefit both batch and streaming
- Estimated effort: 1-2 weeks

## Why Commit Now?

The user confirmed: **"commiting is okay. this suggests the workflow and archive are ok."**

Reasons to commit the current state:
1. ✅ **Correctness is verified** - Data roundtrips perfectly
2. ✅ **C++ AGC compatibility works** - Bidirectional format compatibility
3. ✅ **Constant memory achieved** - Queue-based architecture working
4. ✅ **Architecture is sound** - Just need to add compression step
5. ✅ **All tests pass** - 7 systematic tests validating correctness

The compression is a straightforward addition that doesn't require architectural changes.

## References

- Batch compression: `ragc-core/src/worker.rs:78-100`
- Streaming segment writing: `ragc-core/src/streaming_compressor_queue.rs:488-502`
- Segment compression utils: `ragc-core/src/segment_compression.rs`
- ZSTD pool: `ragc-core/src/zstd_pool.rs`

