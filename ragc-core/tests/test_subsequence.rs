#![allow(clippy::all)]
// Integration test for subsequence extraction
// Tests get_contig_length() and get_contig_range()

use ahash::AHashSet;
use ragc_core::{Decompressor, DecompressorConfig, StreamingQueueCompressor, StreamingQueueConfig};
use std::fs;

fn create_test_archive(path: &str, contigs: Vec<(&str, &str, Vec<u8>)>) {
    // Create compressor
    let config = StreamingQueueConfig {
        queue_capacity: 10 * 1024 * 1024,
        num_threads: 1,
        verbosity: 0,
        ..Default::default()
    };

    let splitters = AHashSet::new();
    let mut compressor = StreamingQueueCompressor::with_splitters(path, config, splitters)
        .expect("Failed to create compressor");

    // Push all contigs
    for (sample, contig, data) in contigs {
        compressor
            .push(sample.to_string(), contig.to_string(), data)
            .expect("Failed to push contig");
    }

    compressor.finalize().expect("Failed to finalize");
}

#[test]
fn test_get_contig_length() {
    let archive_path = "/tmp/test_contig_length.agc";

    // Create test data with known lengths
    let contig1_len = 5000;
    let contig2_len = 10000;

    let contigs = vec![
        ("sample1", "chr1", vec![0u8; contig1_len]), // A x 5000
        ("sample1", "chr2", vec![1u8; contig2_len]), // C x 10000
    ];

    create_test_archive(archive_path, contigs);

    // Test length retrieval
    let config = DecompressorConfig { verbosity: 0 };
    let mut dec = Decompressor::open(archive_path, config).expect("Failed to open archive");

    let len1 = dec
        .get_contig_length("sample1", "chr1")
        .expect("Failed to get length");
    let len2 = dec
        .get_contig_length("sample1", "chr2")
        .expect("Failed to get length");

    assert_eq!(len1, contig1_len, "chr1 length mismatch");
    assert_eq!(len2, contig2_len, "chr2 length mismatch");

    // Clean up
    fs::remove_file(archive_path).ok();
}

#[test]
fn test_get_contig_range_basic() {
    let archive_path = "/tmp/test_contig_range_basic.agc";

    // Create a sequence with different bases at different positions
    // Pattern: AAAAA...CCCCC...GGGGG...TTTTT... (1000 each)
    let mut seq = Vec::with_capacity(4000);
    seq.extend(vec![0u8; 1000]); // A
    seq.extend(vec![1u8; 1000]); // C
    seq.extend(vec![2u8; 1000]); // G
    seq.extend(vec![3u8; 1000]); // T

    let contigs = vec![("sample1", "chr1", seq.clone())];
    create_test_archive(archive_path, contigs);

    let config = DecompressorConfig { verbosity: 0 };
    let mut dec = Decompressor::open(archive_path, config).expect("Failed to open archive");

    // Get full contig for comparison
    let full_contig = dec
        .get_contig("sample1", "chr1")
        .expect("Failed to get full contig");

    // Test various ranges
    let test_cases = vec![
        (0, 100),     // Beginning
        (500, 600),   // Middle of first section
        (900, 1100),  // Across A/C boundary
        (2000, 3000), // Middle (all G)
        (3900, 4000), // End
        (0, 4000),    // Full range
        (1000, 1001), // Single base
    ];

    for (start, end) in test_cases {
        let range_result = dec
            .get_contig_range("sample1", "chr1", start, end)
            .expect(&format!("Failed to get range {}..{}", start, end));

        let expected = &full_contig[start..end];

        assert_eq!(
            range_result.len(),
            expected.len(),
            "Length mismatch for range {}..{}",
            start,
            end
        );
        assert_eq!(
            range_result, expected,
            "Content mismatch for range {}..{}",
            start, end
        );
    }

    // Clean up
    fs::remove_file(archive_path).ok();
}

#[test]
fn test_get_contig_range_edge_cases() {
    let archive_path = "/tmp/test_contig_range_edge.agc";

    let seq = vec![0u8; 1000];
    let contigs = vec![("sample1", "chr1", seq)];
    create_test_archive(archive_path, contigs);

    let config = DecompressorConfig { verbosity: 0 };
    let mut dec = Decompressor::open(archive_path, config).expect("Failed to open archive");

    // Empty range (start >= end)
    let empty1 = dec
        .get_contig_range("sample1", "chr1", 100, 100)
        .expect("Failed to get empty range");
    assert!(empty1.is_empty(), "start==end should return empty");

    let empty2 = dec
        .get_contig_range("sample1", "chr1", 200, 100)
        .expect("Failed to get reverse range");
    assert!(empty2.is_empty(), "start>end should return empty");

    // End beyond contig length (should clamp)
    let clamped = dec
        .get_contig_range("sample1", "chr1", 900, 2000)
        .expect("Failed to get clamped range");
    assert_eq!(clamped.len(), 100, "Should clamp to contig length");

    // Clean up
    fs::remove_file(archive_path).ok();
}

#[test]
fn test_get_contig_range_large_contig() {
    let archive_path = "/tmp/test_contig_range_large.agc";

    // Create a 500KB contig (multiple segments)
    let contig_size = 500_000;
    let mut seq = Vec::with_capacity(contig_size);
    for i in 0..contig_size {
        seq.push((i % 4) as u8);
    }

    let contigs = vec![("sample1", "chr1", seq.clone())];
    create_test_archive(archive_path, contigs);

    let config = DecompressorConfig { verbosity: 0 };
    let mut dec = Decompressor::open(archive_path, config).expect("Failed to open archive");

    // Get full contig for comparison
    let full_contig = dec
        .get_contig("sample1", "chr1")
        .expect("Failed to get full contig");
    assert_eq!(full_contig.len(), contig_size);

    // Test ranges across segment boundaries (~60KB segments)
    let test_cases = vec![
        (0, 1000),        // First segment only
        (59000, 61000),   // Across first segment boundary
        (119000, 121000), // Across second segment boundary
        (250000, 260000), // Middle
        (490000, 500000), // End
    ];

    for (start, end) in test_cases {
        let range_result = dec
            .get_contig_range("sample1", "chr1", start, end)
            .expect(&format!("Failed to get range {}..{}", start, end));

        let expected = &full_contig[start..end];

        assert_eq!(
            range_result.len(),
            expected.len(),
            "Length mismatch for range {}..{}: got {}, expected {}",
            start,
            end,
            range_result.len(),
            expected.len()
        );
        assert_eq!(
            range_result, expected,
            "Content mismatch for range {}..{}",
            start, end
        );
    }

    // Verify length calculation
    let reported_len = dec
        .get_contig_length("sample1", "chr1")
        .expect("Failed to get length");
    assert_eq!(reported_len, contig_size, "Length calculation mismatch");

    // Clean up
    fs::remove_file(archive_path).ok();
}
