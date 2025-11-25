// Collection compatibility tests
// Read collections created by C++ implementation to validate format compatibility
//
// NOTE: These tests require pre-generated C++ collection files in ../test-data/
// They are ignored by default. See ragc-core/tests/cpp_compat.rs for active C++ compatibility tests.

use ragc_core::{Archive, CollectionV3};

#[test]
#[ignore] // Requires pre-generated test fixtures
fn test_read_cpp_simple_collection() {
    let mut archive = Archive::new_reader();
    archive
        .open("../test-data/collection_test1.agc")
        .expect("Failed to open C++ generated collection");

    let mut collection = CollectionV3::new();
    collection.set_config(1000, 21, None);
    collection
        .prepare_for_decompression(&archive)
        .expect("Failed to prepare for decompression");

    collection
        .load_batch_sample_names(&mut archive)
        .expect("Failed to load sample names");

    assert_eq!(collection.get_no_samples(), 2);

    let samples = collection.get_samples_list(false);
    assert_eq!(samples, vec!["sample1", "sample2"]);

    // Load contig batch
    collection
        .load_contig_batch(&mut archive, 0)
        .expect("Failed to load contig batch");

    // Check sample1 contigs
    assert_eq!(collection.get_no_contigs("sample1"), Some(2));
    let contigs = collection.get_contig_list("sample1").unwrap();
    assert_eq!(contigs, vec!["chr1", "chr2"]);

    // Check sample2 contigs
    assert_eq!(collection.get_no_contigs("sample2"), Some(1));
    let contigs = collection.get_contig_list("sample2").unwrap();
    assert_eq!(contigs, vec!["chr1"]);

    // Check segment descriptors for sample1
    let sample_desc = collection.get_sample_desc("sample1").unwrap();
    assert_eq!(sample_desc.len(), 2);

    // sample1/chr1 should have 2 segments
    let chr1_segments = &sample_desc[0].1;
    assert_eq!(chr1_segments.len(), 2);
    assert_eq!(chr1_segments[0].group_id, 0);
    assert_eq!(chr1_segments[0].in_group_id, 0);
    assert!(!chr1_segments[0].is_rev_comp);
    assert_eq!(chr1_segments[0].raw_length, 1021);

    assert_eq!(chr1_segments[1].group_id, 1);
    assert_eq!(chr1_segments[1].in_group_id, 0);
    assert!(!chr1_segments[1].is_rev_comp);
    assert_eq!(chr1_segments[1].raw_length, 1021);

    // sample1/chr2 should have 1 segment
    let chr2_segments = &sample_desc[1].1;
    assert_eq!(chr2_segments.len(), 1);
    assert_eq!(chr2_segments[0].group_id, 2);
    assert_eq!(chr2_segments[0].in_group_id, 0);
    assert!(chr2_segments[0].is_rev_comp);
    assert_eq!(chr2_segments[0].raw_length, 1021);

    // Check segment descriptors for sample2
    let sample_desc = collection.get_sample_desc("sample2").unwrap();
    assert_eq!(sample_desc.len(), 1);

    // sample2/chr1 should have 1 segment
    let chr1_segments = &sample_desc[0].1;
    assert_eq!(chr1_segments.len(), 1);
    assert_eq!(chr1_segments[0].group_id, 0);
    assert_eq!(chr1_segments[0].in_group_id, 1);
    assert!(!chr1_segments[0].is_rev_comp);
    assert_eq!(chr1_segments[0].raw_length, 1021);
}

#[test]
#[ignore] // Requires pre-generated test fixtures
fn test_read_cpp_collection_segments() {
    let mut archive = Archive::new_reader();
    archive
        .open("../test-data/collection_test2.agc")
        .expect("Failed to open C++ generated collection");

    let mut collection = CollectionV3::new();
    collection.set_config(1000, 21, None);
    collection
        .prepare_for_decompression(&archive)
        .expect("Failed to prepare for decompression");

    collection
        .load_batch_sample_names(&mut archive)
        .expect("Failed to load sample names");

    assert_eq!(collection.get_no_samples(), 1);

    let samples = collection.get_samples_list(false);
    assert_eq!(samples, vec!["sample1"]);

    // Load contig batch
    collection
        .load_contig_batch(&mut archive, 0)
        .expect("Failed to load contig batch");

    // Get segment descriptors
    let sample_desc = collection.get_sample_desc("sample1").unwrap();
    assert_eq!(sample_desc.len(), 1);

    let segments = &sample_desc[0].1;
    assert_eq!(segments.len(), 10);

    // Verify first few segments match C++ output
    assert_eq!(segments[0].group_id, 0);
    assert_eq!(segments[0].in_group_id, 0);
    assert!(segments[0].is_rev_comp);
    assert_eq!(segments[0].raw_length, 1021);

    assert_eq!(segments[1].group_id, 1);
    assert_eq!(segments[1].in_group_id, 1);
    assert!(!segments[1].is_rev_comp);
    assert_eq!(segments[1].raw_length, 1021);

    assert_eq!(segments[2].group_id, 2);
    assert_eq!(segments[2].in_group_id, 2);
    assert!(segments[2].is_rev_comp);
    assert_eq!(segments[2].raw_length, 1021);

    // Verify pattern continues
    #[allow(clippy::needless_range_loop)]
    for i in 0..10 {
        assert_eq!(segments[i].group_id, i as u32);
        assert_eq!(segments[i].in_group_id, (i % 3) as u32);
        assert_eq!(segments[i].is_rev_comp, i % 2 == 0);
        assert_eq!(segments[i].raw_length, 1021);
    }
}
