// Archive compatibility tests
// Read archives created by C++ implementation to validate format compatibility
//
// NOTE: These tests require pre-generated C++ archive files in ../test-data/
// They are ignored by default. See ragc-core/tests/cpp_compat.rs for active C++ compatibility tests.

use ragc_core::Archive;

#[test]
#[ignore] // Requires pre-generated test fixtures
fn test_read_cpp_simple_archive() {
    let mut archive = Archive::new_reader();
    archive
        .open("../test-data/archive_test1.agc")
        .expect("Failed to open C++ generated archive");

    let stream_id = archive
        .get_stream_id("test_stream")
        .expect("Stream not found");

    assert_eq!(archive.get_num_parts(stream_id), 2);
    assert_eq!(archive.get_raw_size(stream_id), 100);

    let (data1, meta1) = archive.get_part(stream_id).unwrap().unwrap();
    assert_eq!(data1, b"Hello");
    assert_eq!(meta1, 42);

    let (data2, meta2) = archive.get_part(stream_id).unwrap().unwrap();
    assert_eq!(data2, b"World");
    assert_eq!(meta2, 99);

    assert!(archive.get_part(stream_id).unwrap().is_none());
}

#[test]
#[ignore] // Requires pre-generated test fixtures
fn test_read_cpp_multiple_streams() {
    let mut archive = Archive::new_reader();
    archive
        .open("../test-data/archive_test2.agc")
        .expect("Failed to open C++ generated archive");

    assert_eq!(archive.get_num_streams(), 2);

    let stream1 = archive.get_stream_id("stream1").expect("stream1 not found");
    let stream2 = archive.get_stream_id("stream2").expect("stream2 not found");

    assert_eq!(archive.get_num_parts(stream1), 2);
    assert_eq!(archive.get_num_parts(stream2), 1);

    let (data, meta) = archive.get_part_by_id(stream1, 0).unwrap();
    assert_eq!(data, b"Data1");
    assert_eq!(meta, 1);

    let (data, meta) = archive.get_part_by_id(stream2, 0).unwrap();
    assert_eq!(data, b"Data2");
    assert_eq!(meta, 2);

    let (data, meta) = archive.get_part_by_id(stream1, 1).unwrap();
    assert_eq!(data, b"Data3");
    assert_eq!(meta, 3);
}

#[test]
#[ignore] // Requires pre-generated test fixtures
fn test_read_cpp_large_metadata() {
    let mut archive = Archive::new_reader();
    archive
        .open("../test-data/archive_test3.agc")
        .expect("Failed to open C++ generated archive");

    let stream_id = archive
        .get_stream_id("meta_stream")
        .expect("Stream not found");

    let (_, meta0) = archive.get_part_by_id(stream_id, 0).unwrap();
    assert_eq!(meta0, 0);

    let (_, meta1) = archive.get_part_by_id(stream_id, 1).unwrap();
    assert_eq!(meta1, 255);

    let (_, meta2) = archive.get_part_by_id(stream_id, 2).unwrap();
    assert_eq!(meta2, 65536);
}
