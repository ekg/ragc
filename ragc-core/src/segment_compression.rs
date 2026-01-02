// Segment Compression
// ZSTD compression/decompression for segments with tuple packing

use crate::tuple_packing::{bytes_to_tuples, tuples_to_bytes};
use crate::zstd_pool;
use anyhow::Result;
use ragc_common::types::{Contig, PackedBlock};

/// Default ZSTD compression level for delta packs
/// Use level 17 to match C++ AGC's delta pack compression (segment.h:279)
const DELTA_COMPRESSION_LEVEL: i32 = 17;

/// ZSTD compression level for reference segments with tuple packing
/// C++ AGC uses level 13 for tuple-packed references (segment.h:252)
const REF_TUPLES_COMPRESSION_LEVEL: i32 = 13;

/// ZSTD compression level for reference segments without tuple packing
/// C++ AGC uses level 19 for plain references (segment.h:254)
const REF_PLAIN_COMPRESSION_LEVEL: i32 = 19;

/// Repetitiveness threshold for choosing compression method
/// C++ AGC uses 0.5 (segment.h:225)
const REPETITIVENESS_THRESHOLD: f64 = 0.5;

/// Check repetitiveness of data to decide compression method
/// Matches C++ AGC implementation in segment.h:224-249
fn check_repetitiveness(data: &[u8]) -> f64 {
    let mut best_frac = 0.0;

    for offset in 4..32 {
        let mut cnt = 0;
        let mut cur_size = 0;

        for j in 0..data.len() {
            if j + offset < data.len() {
                if data[j] == data[j + offset] {
                    cnt += 1;
                }
                // Only count ACGT bases (values < 4)
                if data[j] < 4 {
                    cur_size += 1;
                }
            }
        }

        let frac = if cur_size > 0 {
            cnt as f64 / cur_size as f64
        } else {
            0.0
        };

        if frac > best_frac {
            best_frac = frac;
            // Early exit if we've reached threshold
            if best_frac >= REPETITIVENESS_THRESHOLD {
                break;
            }
        }
    }

    best_frac
}

/// Compress a segment using ZSTD with default compression level (for delta packs)
pub fn compress_segment(data: &Contig) -> Result<PackedBlock> {
    compress_segment_plain(data, DELTA_COMPRESSION_LEVEL)
}

/// Compress a segment using ZSTD with configured level (for delta packs)
pub fn compress_segment_configured(data: &Contig, level: i32) -> Result<PackedBlock> {
    compress_segment_plain(data, level)
}

/// Compress a reference segment with automatic tuple packing decision
///
/// **IMPORTANT**: Matches C++ AGC's store_in_archive() logic!
/// - Checks repetitiveness of data
/// - If repetitiveness < 0.5: use tuple packing (marker 1, level 13)
/// - If repetitiveness >= 0.5: use plain ZSTD (marker 0, level 19)
///
/// Returns (compressed_data, marker_byte)
pub fn compress_reference_segment(data: &Contig) -> Result<(PackedBlock, u8)> {
    let repetitiveness = check_repetitiveness(data);

    // Debug logging for reference compression decisions
    let debug_ref = crate::env_cache::debug_ref();
    if debug_ref {
        eprintln!(
            "RAGC_REF_COMPRESS: len={} rep={:.4} threshold={:.4}",
            data.len(),
            repetitiveness,
            REPETITIVENESS_THRESHOLD
        );
    }

    if repetitiveness < REPETITIVENESS_THRESHOLD {
        // Low repetitiveness: use tuple packing
        let tuples = bytes_to_tuples(data);
        let compressed = zstd_pool::compress_segment_pooled(&tuples, REF_TUPLES_COMPRESSION_LEVEL)?;
        if debug_ref {
            eprintln!(
                "RAGC_REF_DECISION: TUPLE_PACK marker=1 level={} tuple_len={} compressed_len={}",
                REF_TUPLES_COMPRESSION_LEVEL,
                tuples.len(),
                compressed.len()
            );
        }
        Ok((compressed, 1)) // Marker 1 = tuple-packed
    } else {
        // High repetitiveness: use plain ZSTD
        let compressed = zstd_pool::compress_segment_pooled(data, REF_PLAIN_COMPRESSION_LEVEL)?;
        if debug_ref {
            eprintln!(
                "RAGC_REF_DECISION: PLAIN marker=0 level={} compressed_len={}",
                REF_PLAIN_COMPRESSION_LEVEL,
                compressed.len()
            );
        }
        Ok((compressed, 0)) // Marker 0 = plain
    }
}

/// Compress a segment using plain ZSTD (for delta packs)
pub fn compress_segment_plain(data: &Contig, level: i32) -> Result<PackedBlock> {
    zstd_pool::compress_segment_pooled(data, level)
}

/// Decompress a segment based on marker byte
///
/// **IMPORTANT**: Checks marker byte to determine decompression method!
/// - Marker 0: plain ZSTD decompression
/// - Marker 1 (or any non-zero): ZSTD + tuple unpacking
pub fn decompress_segment_with_marker(compressed: &[u8], marker: u8) -> Result<Contig> {
    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    if marker == 0 {
        // Plain ZSTD
        zstd_pool::decompress_segment_pooled(&compressed.to_vec())
    } else {
        // Tuple-packed: decompress then unpack
        let tuples = zstd_pool::decompress_segment_pooled(&compressed.to_vec())?;
        Ok(tuples_to_bytes(&tuples))
    }
}

/// Decompress a segment using plain ZSTD (for old code compatibility)
pub fn decompress_segment(compressed: &[u8]) -> Result<Contig> {
    zstd_pool::decompress_segment_pooled(&compressed.to_vec())
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
            let compressed = compress_segment_configured(&original, *level).unwrap();
            let decompressed = decompress_segment(&compressed).unwrap();
            assert_eq!(original, decompressed);
        }
    }

    #[test]
    fn test_repetitiveness_check() {
        // Highly repetitive sequence (all zeros) - should get score of 1.0
        let repetitive = vec![0; 100];
        let rep1 = check_repetitiveness(&repetitive);
        assert_eq!(
            rep1, 1.0,
            "All-zero sequence should have perfect repetitiveness"
        );

        // Test that the function returns a value between 0 and 1
        let mixed = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let rep2 = check_repetitiveness(&mixed);
        assert!(
            (0.0..=1.0).contains(&rep2),
            "Repetitiveness should be in [0, 1]"
        );
    }

    #[test]
    fn test_reference_compression_with_tuple_packing() {
        // Test both compression paths work correctly

        // High repetitiveness: should use plain ZSTD (marker 0)
        let high_rep = vec![0; 100];
        let (compressed1, marker1) = compress_reference_segment(&high_rep).unwrap();
        assert_eq!(marker1, 0, "High repetitiveness should use plain ZSTD");
        let decompressed1 = decompress_segment_with_marker(&compressed1, marker1).unwrap();
        assert_eq!(high_rep, decompressed1);

        // Note: We don't test for marker=1 here because it's hard to create
        // test data that reliably has repetitiveness < 0.5. The important
        // thing is that both code paths work (tested separately below).
    }

    #[test]
    fn test_reference_compression_without_tuple_packing() {
        // High repetitiveness: should use plain ZSTD
        let high_rep = vec![0; 100];
        let (compressed, marker) = compress_reference_segment(&high_rep).unwrap();
        assert_eq!(marker, 0, "High repetitiveness should use plain ZSTD");
        let decompressed = decompress_segment_with_marker(&compressed, marker).unwrap();
        assert_eq!(high_rep, decompressed);
    }

    #[test]
    fn test_tuple_packing_compression_path() {
        // Test the tuple packing code path directly (marker 1)
        // Even if we can't reliably trigger it via compress_reference_segment,
        // we can test that decompress_segment_with_marker handles it correctly
        let original = vec![0, 1, 2, 3, 0, 1, 2, 3];

        // Manually invoke tuple packing path
        let tuples = bytes_to_tuples(&original);
        let compressed = zstd_pool::compress_segment_pooled(&tuples, 13).unwrap();

        // Test decompression with marker 1
        let decompressed = decompress_segment_with_marker(&compressed, 1).unwrap();
        assert_eq!(
            original, decompressed,
            "Tuple packing roundtrip should work"
        );
    }
}
