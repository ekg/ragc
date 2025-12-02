// Test that C++ AGC FFI archive creation works

use ragc_core::agc_index_ffi::{CppArchive, CppCollection, SegmentPlacement};

#[test]
fn test_cpp_archive_lifecycle() {
    let test_path = "/tmp/test_cpp_ffi_archive.agc";

    // Create and open archive
    let mut archive = CppArchive::new_writer();
    archive.open(test_path).expect("Failed to open archive");

    // Create collection
    let mut collection = CppCollection::new();
    collection
        .set_archives(&mut archive, 1, 1 << 20, 10000, 21)
        .expect("Failed to set archives");

    // CRITICAL: Register sample/contig BEFORE adding placements
    collection
        .register_sample_contig("test_sample", "test_contig")
        .expect("Failed to register sample/contig");

    // Add a test segment placement
    let placements = vec![SegmentPlacement {
        sample_name: "test_sample".to_string(),
        contig_name: "test_contig".to_string(),
        seg_part_no: 0,
        group_id: 16,
        in_group_id: 0,
        is_rev_comp: false,
        data_size: 1000,
    }];

    collection
        .add_segments_placed(&placements)
        .expect("Failed to add placements");

    // Finalize
    collection.complete_serialization();
    collection.store_contig_batch(0, 1);
    archive.close().expect("Failed to close archive");

    // Verify file exists
    assert!(std::path::Path::new(test_path).exists());

    // Cleanup
    std::fs::remove_file(test_path).ok();
}

#[test]
fn test_cpp_archive_multiple_samples() {
    let test_path = "/tmp/test_cpp_multi_sample.agc";

    // Create archive
    let mut archive = CppArchive::new_writer();
    archive.open(test_path).expect("Failed to open archive");

    let mut collection = CppCollection::new();
    collection
        .set_archives(&mut archive, 1, 1 << 20, 10000, 21)
        .expect("Failed to set archives");

    // Register multiple samples and contigs
    collection
        .register_sample_contig("sample1", "chr1")
        .expect("Failed to register sample1/chr1");
    collection
        .register_sample_contig("sample1", "chr2")
        .expect("Failed to register sample1/chr2");
    collection
        .register_sample_contig("sample2", "chr1")
        .expect("Failed to register sample2/chr1");

    // Add placements for each
    let placements = vec![
        SegmentPlacement {
            sample_name: "sample1".to_string(),
            contig_name: "chr1".to_string(),
            seg_part_no: 0,
            group_id: 16,
            in_group_id: 0,
            is_rev_comp: false,
            data_size: 1000,
        },
        SegmentPlacement {
            sample_name: "sample1".to_string(),
            contig_name: "chr2".to_string(),
            seg_part_no: 0,
            group_id: 17,
            in_group_id: 0,
            is_rev_comp: false,
            data_size: 2000,
        },
        SegmentPlacement {
            sample_name: "sample2".to_string(),
            contig_name: "chr1".to_string(),
            seg_part_no: 0,
            group_id: 16,
            in_group_id: 1,
            is_rev_comp: false,
            data_size: 1000,
        },
    ];

    collection
        .add_segments_placed(&placements)
        .expect("Failed to add placements");

    // Finalize
    collection.complete_serialization();
    collection.store_contig_batch(0, 2);
    archive.close().expect("Failed to close archive");

    // Verify file exists
    assert!(std::path::Path::new(test_path).exists());

    // Cleanup
    std::fs::remove_file(test_path).ok();
}
