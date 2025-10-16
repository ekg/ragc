// Full round-trip test: FASTA → AGC → FASTA
// Verifies that we can compress and decompress without data loss

use ragc_core::{Compressor, CompressorConfig, Decompressor, DecompressorConfig};
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[test]
fn test_roundtrip_fasta_output() {
    let archive_path = "/tmp/test_roundtrip_output.agc";
    let output_fasta = "/tmp/test_roundtrip_output.fasta";
    let _ = fs::remove_file(archive_path);
    let _ = fs::remove_file(output_fasta);

    // Create archive with known data
    {
        let config = CompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            verbosity: 0,
        };

        let mut compressor = Compressor::new(archive_path, config).unwrap();

        // Add a simple test sequence: ACGTACGT repeated
        let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGTACGTACGT
        let seq2 = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0]; // TGCATGCATGCA

        compressor
            .add_contig("test_sample", "chr1", seq1.clone())
            .unwrap();
        compressor
            .add_contig("test_sample", "chr2", seq2.clone())
            .unwrap();

        compressor.finalize().unwrap();
    }

    // Decompress and write to FASTA
    {
        let config = DecompressorConfig { verbosity: 0 };
        let mut decompressor = Decompressor::open(archive_path, config).unwrap();

        decompressor
            .write_sample_fasta("test_sample", Path::new(output_fasta))
            .unwrap();
        decompressor.close().unwrap();
    }

    // Verify FASTA output
    {
        let file = fs::File::open(output_fasta).unwrap();
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();

        eprintln!("Output FASTA lines:");
        for (i, line) in lines.iter().enumerate() {
            eprintln!("  [{i}]: {line}");
        }

        assert_eq!(
            lines.len(),
            4,
            "Expected 4 lines, got {}: {:?}",
            lines.len(),
            lines
        ); // 2 headers + 2 sequences

        // Check chr1
        assert_eq!(lines[0], ">chr1");
        assert_eq!(lines[1], "ACGTACGTACGTACGT");

        // Check chr2
        assert_eq!(lines[2], ">chr2");
        assert_eq!(lines[3], "TGCATGCATGCA");
    }

    fs::remove_file(archive_path).unwrap();
    fs::remove_file(output_fasta).unwrap();
}

#[test]
fn test_roundtrip_real_fasta() {
    let test_fasta = Path::new("../test-data/test_simple.fasta");

    // Skip if test file doesn't exist
    if !test_fasta.exists() {
        eprintln!("Skipping test_roundtrip_real_fasta: test file not found");
        return;
    }

    let archive_path = "/tmp/test_roundtrip_real.agc";
    let output_fasta = "/tmp/test_roundtrip_real_output.fasta";
    let _ = fs::remove_file(archive_path);
    let _ = fs::remove_file(output_fasta);

    // Compress
    {
        let config = CompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            verbosity: 0,
        };

        let mut compressor = Compressor::new(archive_path, config).unwrap();
        compressor
            .add_fasta_file("test_sample", test_fasta)
            .unwrap();
        compressor.finalize().unwrap();
    }

    // Decompress
    {
        let config = DecompressorConfig { verbosity: 0 };
        let mut decompressor = Decompressor::open(archive_path, config).unwrap();
        decompressor
            .write_sample_fasta("test_sample", Path::new(output_fasta))
            .unwrap();
        decompressor.close().unwrap();
    }

    // Compare input and output (ignoring whitespace differences)
    {
        let original = fs::read_to_string(test_fasta).unwrap();
        let output = fs::read_to_string(output_fasta).unwrap();

        // Extract sequences (skip headers, concatenate all sequence lines)
        let original_seq: String = original
            .lines()
            .filter(|l| !l.starts_with('>'))
            .map(|l| l.trim())
            .collect();

        let output_seq: String = output
            .lines()
            .filter(|l| !l.starts_with('>'))
            .map(|l| l.trim())
            .collect();

        assert_eq!(
            original_seq.to_uppercase(),
            output_seq.to_uppercase(),
            "Sequence data should match after round-trip"
        );
    }

    fs::remove_file(archive_path).unwrap();
    fs::remove_file(output_fasta).unwrap();
}
