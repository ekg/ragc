// ZSTD Context Pooling - Matching C++ AGC Design
//
// C++ AGC creates ONE ZSTD_CCtx per thread and reuses it for ALL compressions
// with DIFFERENT levels by passing level as a parameter to ZSTD_compressCCtx().
//
// CRITICAL: Previous implementation used zstd::bulk::Compressor which LOCKS
// the compression level at creation, forcing context recreation when level changes.
// This doesn't match C++ AGC behavior and may affect compression efficiency!
//
// New implementation uses zstd_safe raw bindings to match C++ AGC exactly:
// - ONE thread-local ZSTD_CCtx reused for ALL compressions
// - Compression level passed as parameter to each call (not stored in context)
// - Matches C++ AGC: ZSTD_compressCCtx(ctx, ..., level)

use anyhow::Result;
use ragc_common::types::{Contig, PackedBlock};
use std::cell::RefCell;

thread_local! {
    /// Thread-local ZSTD compression context - EXACTLY like C++ AGC
    /// C++ AGC: ZSTD_CCtx* zstd_ctx = ZSTD_createCCtx() (reused per thread)
    /// Rust: zstd::zstd_safe::CCtx<'static> (same concept)
    ///
    /// CRITICAL: We do NOT store the compression level in the context!
    /// Level is passed as a parameter to each compression call, allowing
    /// the same context to handle levels 13, 17, and 19 without recreation.
    static ZSTD_CCTX: RefCell<Option<zstd::zstd_safe::CCtx<'static>>> = const { RefCell::new(None) };
}

/// Compress a segment using thread-local ZSTD context (MATCHING C++ AGC)
///
/// C++ AGC equivalent (segment.h:172-189):
/// ```cpp
/// size_t a_size = ZSTD_compressBound(data.size());
/// uint8_t *packed = new uint8_t[a_size+1u];
/// uint32_t packed_size = ZSTD_compressCCtx(zstd_ctx, packed, a_size, data.data(), data.size(), level);
/// vector<uint8_t> v_packed(packed, packed + packed_size + 1);
/// delete[] packed;
/// ```
///
/// Our Rust implementation does the same with thread-local context.
pub fn compress_segment_pooled(data: &Contig, level: i32) -> Result<PackedBlock> {
    ZSTD_CCTX.with(|ctx_cell| {
        let mut ctx_opt = ctx_cell.borrow_mut();

        // Get or create context (only created once per thread)
        let ctx = match ctx_opt.as_mut() {
            Some(c) => c,
            None => {
                let new_ctx = zstd::zstd_safe::CCtx::create();
                *ctx_opt = Some(new_ctx);
                ctx_opt.as_mut().unwrap()
            }
        };

        // Calculate maximum compressed size (matching C++ AGC ZSTD_compressBound)
        let max_size = zstd::zstd_safe::compress_bound(data.len());

        // Allocate output buffer (matching C++ AGC: new uint8_t[a_size])
        let mut compressed = vec![0u8; max_size];

        // Compress with level parameter (matching C++ AGC: ZSTD_compressCCtx)
        // CRITICAL: Level is passed HERE, not stored in context!
        let compressed_size = ctx
            .compress(&mut compressed, data, level)
            .map_err(|e| anyhow::anyhow!("ZSTD compression failed: {e}"))?;

        // Truncate to actual compressed size
        compressed.truncate(compressed_size);
        Ok(compressed)
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
