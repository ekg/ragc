// Segment Compression
// ZSTD compression/decompression for segments

use anyhow::{Context, Result};
use ragc_common::types::{Contig, PackedBlock};

/// Default ZSTD compression level
/// Use level 17 to match C++ AGC's delta pack compression (segment.h:279)
/// C++ AGC uses levels 13-19 depending on segment type:
///   - Reference segments: level 13 (tuples) or 19 (plain) based on repetitiveness
///   - Delta packs: level 17
/// We use 17 as it matches C++ AGC's most common case and provides good balance
const DEFAULT_COMPRESSION_LEVEL: i32 = 17;

/// Compress a segment using ZSTD with default compression level
pub fn compress_segment(data: &Contig) -> Result<PackedBlock> {
    compress_segment_with_level(data, DEFAULT_COMPRESSION_LEVEL)
}

/// Compress a segment using ZSTD with compression level from config
pub fn compress_segment_configured(data: &Contig, level: i32) -> Result<PackedBlock> {
    compress_segment_with_level(data, level)
}

/// Compress a segment using ZSTD with specified compression level
pub fn compress_segment_with_level(data: &Contig, level: i32) -> Result<PackedBlock> {
    zstd::encode_all(data.as_slice(), level).context("Failed to compress segment with ZSTD")
}

/// Decompress a segment using ZSTD
pub fn decompress_segment(compressed: &PackedBlock) -> Result<Contig> {
    zstd::decode_all(compressed.as_slice()).context("Failed to decompress segment with ZSTD")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compress_decompress_roundtrip() {
        let original = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];

        let compressed = compress_segment(&original).unwrap();
        let decompressed = decompress_segment(&compressed).unwrap();

        assert_eq!(original, decompressed);
    }

    #[test]
    fn test_compress_empty() {
        let original = vec![];

        let compressed = compress_segment(&original).unwrap();
        let decompressed = decompress_segment(&compressed).unwrap();

        assert_eq!(original, decompressed);
    }

    #[test]
    fn test_compress_large() {
        // Create a large sequence with some repetition
        let mut original = Vec::new();
        for i in 0..1000 {
            original.push((i % 4) as u8);
        }

        let compressed = compress_segment(&original).unwrap();
        let decompressed = decompress_segment(&compressed).unwrap();

        assert_eq!(original, decompressed);

        // Compressed should be smaller than original
        assert!(compressed.len() < original.len());
    }

    #[test]
    fn test_different_compression_levels() {
        let original = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];

        for level in [1, 3, 9, 19].iter() {
            let compressed = compress_segment_with_level(&original, *level).unwrap();
            let decompressed = decompress_segment(&compressed).unwrap();
            assert_eq!(original, decompressed);
        }
    }
}
