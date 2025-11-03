// Integration test for StreamingQueueCompressor
// Verifies basic push() API and worker thread operation

use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
use std::collections::HashSet;

#[test]
fn test_streaming_queue_basic_flow() {
    // Create a small test dataset
    let test_data = vec![
        ("sample1".to_string(), "chr1".to_string(), vec![0u8; 1000]), // A x 1000
        ("sample1".to_string(), "chr2".to_string(), vec![1u8; 1000]), // C x 1000
        ("sample2".to_string(), "chr1".to_string(), vec![2u8; 1000]), // G x 1000
    ];

    // Create compressor with small queue for testing
    let config = StreamingQueueConfig {
        queue_capacity: 10 * 1024 * 1024, // 10 MB for test
        num_threads: 2,
        verbosity: 0, // Quiet for tests
        ..Default::default()
    };

    let splitters = HashSet::new(); // Will be determined from first contig
    let mut compressor = StreamingQueueCompressor::with_splitters(
        "/tmp/test_streaming_queue.agc",
        config,
        splitters,
    )
    .expect("Failed to create compressor");

    // Push all contigs
    for (sample, contig, data) in test_data {
        compressor
            .push(sample, contig, data)
            .expect("Failed to push contig");
    }

    // Finalize (waits for workers)
    compressor.finalize().expect("Failed to finalize");

    // If we get here, the basic pipeline worked!
    // Workers pulled contigs, split them, compressed them
    // (Archive writing is TODO, but the core flow works)
}

#[test]
fn test_streaming_queue_stats() {
    let config = StreamingQueueConfig {
        queue_capacity: 1024,
        num_threads: 1,
        verbosity: 0,
        ..Default::default()
    };

    let compressor = StreamingQueueCompressor::with_splitters(
        "/tmp/test_queue_stats.agc",
        config,
        HashSet::new(),
    )
    .expect("Failed to create compressor");

    let stats = compressor.queue_stats();
    assert_eq!(stats.current_size_bytes, 0);
    assert_eq!(stats.current_items, 0);
    assert_eq!(stats.capacity_bytes, 1024);
    assert!(!stats.is_closed);
}
