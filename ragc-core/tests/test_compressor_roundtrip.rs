// Round-trip compression test
// Verifies that we can create an archive and read it back

use ragc_common::{Archive, CollectionV3};
use ragc_core::{Compressor, CompressorConfig};
use std::fs;
use std::path::Path;

#[test]
fn test_roundtrip_simple() {
    let archive_path = "/tmp/test_roundtrip.agc";
    let _ = fs::remove_file(archive_path);

    // Create archive
    {
        let config = CompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            verbosity: 0,
        };

        let mut compressor = Compressor::new(archive_path, config).unwrap();

        // Add multiple samples
        let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGTACGT
        let seq2 = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0]; // TGCATGCATGCA

        compressor.add_contig("sample1", "chr1", seq1.clone()).unwrap();
        compressor.add_contig("sample1", "chr2", seq2.clone()).unwrap();
        compressor.add_contig("sample2", "chr1", seq1.clone()).unwrap();

        compressor.finalize().unwrap();
    }

    // Read archive back
    {
        let mut archive = Archive::new_reader();
        archive.open(archive_path).unwrap();

        let mut collection = CollectionV3::new();
        collection.set_config(1000, 21, None);
        collection.prepare_for_decompression(&archive).unwrap();
        collection.load_batch_sample_names(&mut archive).unwrap();

        // Verify samples
        assert_eq!(collection.get_no_samples(), 2);
        let samples = collection.get_samples_list(false);
        assert_eq!(samples.len(), 2);
        assert!(samples.contains(&"sample1".to_string()));
        assert!(samples.contains(&"sample2".to_string()));

        // Load contigs
        collection.load_contig_batch(&mut archive, 0).unwrap();

        // Verify sample1 has 2 contigs
        assert_eq!(collection.get_no_contigs("sample1"), Some(2));
        let contigs = collection.get_contig_list("sample1").unwrap();
        assert_eq!(contigs.len(), 2);
        assert!(contigs.contains(&"chr1".to_string()));
        assert!(contigs.contains(&"chr2".to_string()));

        // Verify sample2 has 1 contig
        assert_eq!(collection.get_no_contigs("sample2"), Some(1));
        let contigs = collection.get_contig_list("sample2").unwrap();
        assert_eq!(contigs.len(), 1);
        assert!(contigs.contains(&"chr1".to_string()));

        // Verify segment descriptors
        let sample1_desc = collection.get_sample_desc("sample1").unwrap();
        assert_eq!(sample1_desc.len(), 2); // 2 contigs

        // Each contig should have 1 segment (no splitters in this test)
        assert_eq!(sample1_desc[0].1.len(), 1); // chr1 has 1 segment
        assert_eq!(sample1_desc[1].1.len(), 1); // chr2 has 1 segment

        archive.close().unwrap();
    }

    fs::remove_file(archive_path).unwrap();
}

#[test]
fn test_roundtrip_fasta() {
    let archive_path = "/tmp/test_roundtrip_fasta.agc";
    let _ = fs::remove_file(archive_path);

    // Create archive from FASTA
    {
        let config = CompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            verbosity: 0,
        };

        let mut compressor = Compressor::new(archive_path, config).unwrap();

        let fasta_path = Path::new("../test-data/test_simple.fasta");
        if fasta_path.exists() {
            compressor.add_fasta_file("test_sample", fasta_path).unwrap();
            compressor.finalize().unwrap();

            // Read back and verify
            let mut archive = Archive::new_reader();
            archive.open(archive_path).unwrap();

            let mut collection = CollectionV3::new();
            collection.set_config(1000, 21, None);
            collection.prepare_for_decompression(&archive).unwrap();
            collection.load_batch_sample_names(&mut archive).unwrap();

            assert_eq!(collection.get_no_samples(), 1);
            let samples = collection.get_samples_list(false);
            assert_eq!(samples, vec!["test_sample"]);

            collection.load_contig_batch(&mut archive, 0).unwrap();
            assert_eq!(collection.get_no_contigs("test_sample"), Some(2));

            let contigs = collection.get_contig_list("test_sample").unwrap();
            assert!(contigs.contains(&"chr1".to_string()));
            assert!(contigs.contains(&"chr2".to_string()));

            archive.close().unwrap();

            fs::remove_file(archive_path).unwrap();
        }
    }
}
