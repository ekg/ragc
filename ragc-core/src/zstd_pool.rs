// ZSTD Compression Helpers
//
// Uses thread-local ZSTD contexts to avoid allocation overhead per-segment.
// This matches how C++ AGC uses ZSTD: one ZSTD_CCtx per thread, reused for all segments.
//
// Note: ZSTD_resetCCtx_internal (the memset) is still called on every compression,
// but we avoid the malloc/free overhead of creating/destroying contexts.

use anyhow::Result;
use ragc_common::types::{Contig, PackedBlock};
use std::cell::RefCell;

// Thread-local ZSTD compression context - matches C++ AGC's per-thread ZSTD_CCtx
thread_local! {
    static ZSTD_CCTX: RefCell<zstd_safe::CCtx<'static>> = RefCell::new(zstd_safe::CCtx::create());
}

/// Compress a segment using ZSTD at the specified level
///
/// Uses thread-local context reuse like C++ AGC's `ZSTD_compressCCtx()`.
/// This avoids per-segment context allocation/deallocation overhead.
///
/// C++ AGC pattern:
///   ZSTD_CCtx* zstd_cctx = ZSTD_createCCtx();  // once per thread
///   ZSTD_compressCCtx(zstd_cctx, ...);         // reused for all segments
///   ZSTD_freeCCtx(zstd_cctx);                  // at thread exit
///
/// This implementation follows the same pattern using thread_local.
pub fn compress_segment_pooled(data: &Contig, level: i32) -> Result<PackedBlock> {
    ZSTD_CCTX.with(|cctx| {
        let mut cctx = cctx.borrow_mut();

        // Pre-allocate output buffer (ZSTD_compressBound equivalent)
        let max_compressed_size = zstd_safe::compress_bound(data.len());
        let mut output = vec![0u8; max_compressed_size];

        // Use ZSTD_compressCCtx - same as C++ AGC
        match cctx.compress(&mut output, data, level) {
            Ok(compressed_size) => {
                output.truncate(compressed_size);
                Ok(output)
            }
            Err(code) => {
                let msg = zstd_safe::get_error_name(code);
                Err(anyhow::anyhow!("ZSTD compression failed: {}", msg))
            }
        }
    })
}

/// Decompress a segment using ZSTD
///
/// Note: Decompression is less critical for pooling since:
/// 1. Decompression contexts are smaller than compression contexts
/// 2. Decompression happens less frequently in the hot path
/// 3. Current implementation is adequate for now
///
/// Future optimization: Add thread-local decoder pool if profiling shows benefit
pub fn decompress_segment_pooled(compressed: &PackedBlock) -> Result<Contig> {
    // For now, use the existing decode_all
    // Decompression contexts are smaller and less critical to pool
    zstd::decode_all(compressed.as_slice())
        .map_err(|e| anyhow::anyhow!("Failed to decompress segment with ZSTD: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pooled_compress_decompress_roundtrip() {
        let original = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];

        let compressed = compress_segment_pooled(&original, 11).unwrap();
        let decompressed = decompress_segment_pooled(&compressed).unwrap();

        assert_eq!(original, decompressed);
    }

    #[test]
    fn test_pooled_multiple_compressions() {
        // Test that context reuse works correctly
        for _ in 0..10 {
            let original = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
            let compressed = compress_segment_pooled(&original, 11).unwrap();
            let decompressed = decompress_segment_pooled(&compressed).unwrap();
            assert_eq!(original, decompressed);
        }
    }

    #[test]
    fn test_pooled_different_levels() {
        let original = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];

        for level in [1, 3, 9, 17, 19].iter() {
            let compressed = compress_segment_pooled(&original, *level).unwrap();
            let decompressed = decompress_segment_pooled(&compressed).unwrap();
            assert_eq!(original, decompressed);
        }
    }

    #[test]
    fn test_pooled_large_data() {
        let mut original = Vec::new();
        for i in 0..10000 {
            original.push((i % 4) as u8);
        }

        let compressed = compress_segment_pooled(&original, 17).unwrap();
        let decompressed = decompress_segment_pooled(&compressed).unwrap();

        assert_eq!(original, decompressed);
        assert!(compressed.len() < original.len());
    }
}
