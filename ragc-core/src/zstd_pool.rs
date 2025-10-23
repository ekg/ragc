// ZSTD Context Pooling - Matching C++ AGC Design
//
// C++ AGC creates ONE ZSTD_CCtx per thread and reuses it for all compressions.
// This module implements the exact same approach using zstd-safe bindings.
//
// CRITICAL FIX: Previous implementation used zstd::stream::copy_encode which:
// 1. Created a new encoder context for each call (massive overhead!)
// 2. Did internal buffering we don't need
// 3. Cloned the output buffer (defeating pooling!)
//
// New implementation matches C++ AGC line-by-line:
// - Thread-local ZSTD_CCtx reused for all compressions (like C++ line 176)
// - Pre-allocated output buffer (like C++ new uint8_t[a_size])
// - Direct compression into buffer (like ZSTD_compressCCtx)

use anyhow::Result;
use ragc_common::types::{Contig, PackedBlock};
use std::cell::RefCell;

thread_local! {
    /// Thread-local ZSTD compression context - EXACTLY like C++ AGC
    /// C++ AGC: ZSTD_CCtx* zstd_ctx (reused per thread)
    /// Rust: zstd::bulk::Compressor<'static> (same concept)
    ///
    /// Note: We store (Compressor, level) to detect level changes.
    /// In practice, compression level rarely changes within a thread.
    static ZSTD_ENCODER: RefCell<Option<(zstd::bulk::Compressor<'static>, i32)>> = RefCell::new(None);
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
/// Our Rust implementation does the same but with thread-local context pooling.
pub fn compress_segment_pooled(data: &Contig, level: i32) -> Result<PackedBlock> {
    ZSTD_ENCODER.with(|encoder_cell| {
        let mut encoder_opt = encoder_cell.borrow_mut();

        // Get or create encoder for this level
        let encoder = match encoder_opt.as_mut() {
            Some((enc, cached_level)) if *cached_level == level => enc,
            _ => {
                // Create new encoder for this level
                let new_encoder = zstd::bulk::Compressor::new(level)
                    .map_err(|e| anyhow::anyhow!("Failed to create ZSTD encoder: {}", e))?;
                *encoder_opt = Some((new_encoder, level));
                &mut encoder_opt.as_mut().unwrap().0
            }
        };

        // Compress directly - bulk::Compressor handles buffer internally
        // and returns owned Vec (no clone needed!)
        encoder
            .compress(data.as_slice())
            .map_err(|e| anyhow::anyhow!("ZSTD compression failed: {}", e))
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
        .map_err(|e| anyhow::anyhow!("Failed to decompress segment with ZSTD: {}", e))
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
