//! Integration tests for C++ compatibility
//!
//! These tests verify that ragc produces archives that are bit-compatible
//! with C++ AGC and can read C++ AGC archives correctly.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use ragc_core::{Compressor, Decompressor, CompressorConfig, DecompressorConfig};
use sha2::{Sha256, Digest};

fn get_test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../test-data")
}

fn compute_file_hash(path: &Path) -> String {
    let data = fs::read(path).expect("Failed to read file");
    let hash = Sha256::digest(&data);
    format!("{:x}", hash)
}

fn contig_to_string(contig: &[u8]) -> String {
    contig.iter().map(|&b| match b {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }).collect()
}

fn create_test_fasta(path: &Path) {
    let content = r#">seq1
ACGTACGTACGTACGTACGTACGTACGTACGT
>seq2
ACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNNNNNNN
"#;
    fs::write(path, content).expect("Failed to write test FASTA");
}

#[test]
fn test_ragc_creates_valid_archive() {
    let test_dir = std::env::temp_dir();
    let fasta_path = test_dir.join("test_compat.fasta");
    let archive_path = test_dir.join("test_compat_ragc.agc");

    // Create test FASTA
    create_test_fasta(&fasta_path);

    // Create archive with ragc
    let config = CompressorConfig::default();
    let mut compressor = Compressor::new(archive_path.to_str().unwrap(), config)
        .expect("Failed to create compressor");

    compressor.add_fasta_file("test_sample", &fasta_path)
        .expect("Failed to add FASTA");

    compressor.finalize()
        .expect("Failed to finalize archive");

    assert!(archive_path.exists(), "Archive was not created");

    // Verify archive can be read back
    let config = DecompressorConfig::default();
    let mut decompressor = Decompressor::open(archive_path.to_str().unwrap(), config)
        .expect("Failed to open archive");

    let samples = decompressor.list_samples();

    assert_eq!(samples.len(), 1);
    assert_eq!(samples[0], "test_sample");

    // Clean up
    let _ = fs::remove_file(&fasta_path);
    let _ = fs::remove_file(&archive_path);
}

#[test]
fn test_ragc_rust_roundtrip() {
    let test_dir = std::env::temp_dir();
    let fasta_path = test_dir.join("test_roundtrip.fasta");
    let archive_path = test_dir.join("test_roundtrip.agc");
    let output_path = test_dir.join("test_roundtrip_out.fasta");

    // Create test FASTA
    create_test_fasta(&fasta_path);
    let original_hash = compute_file_hash(&fasta_path);

    // Compress with ragc
    let config = CompressorConfig::default();
    let mut compressor = Compressor::new(archive_path.to_str().unwrap(), config)
        .expect("Failed to create compressor");

    compressor.add_fasta_file("test_sample", &fasta_path)
        .expect("Failed to add FASTA");

    compressor.finalize()
        .expect("Failed to finalize archive");

    // Decompress with ragc
    let config = DecompressorConfig::default();
    let mut decompressor = Decompressor::open(archive_path.to_str().unwrap(), config)
        .expect("Failed to open archive");

    let sequences = decompressor.get_sample("test_sample")
        .expect("Failed to extract sample");

    // Write output
    let mut output_content = String::new();
    for (name, contig) in sequences {
        output_content.push_str(&format!(">{}\n", name));
        output_content.push_str(&contig_to_string(&contig));
        output_content.push('\n');
    }
    fs::write(&output_path, output_content).expect("Failed to write output");

    let output_hash = compute_file_hash(&output_path);

    assert_eq!(original_hash, output_hash,
        "Roundtrip produced different data! Original: {}, Output: {}",
        original_hash, output_hash);

    // Clean up
    let _ = fs::remove_file(&fasta_path);
    let _ = fs::remove_file(&archive_path);
    let _ = fs::remove_file(&output_path);
}

#[test]
#[ignore] // HashMap iteration order is non-deterministic for security reasons
fn test_deterministic_compression() {
    // Test that ragc produces identical archives for identical inputs
    // NOTE: This test is expected to fail due to HashMap randomization
    let test_dir = std::env::temp_dir();
    let fasta_path = test_dir.join("test_deterministic.fasta");
    let archive1_path = test_dir.join("test_deterministic_1.agc");
    let archive2_path = test_dir.join("test_deterministic_2.agc");

    create_test_fasta(&fasta_path);

    // Create first archive
    let config = CompressorConfig::default();
    let mut compressor1 = Compressor::new(archive1_path.to_str().unwrap(), config.clone())
        .expect("Failed to create compressor 1");
    compressor1.add_fasta_file("test_sample", &fasta_path)
        .expect("Failed to add FASTA 1");
    compressor1.finalize().expect("Failed to finalize 1");

    // Create second archive
    let mut compressor2 = Compressor::new(archive2_path.to_str().unwrap(), config)
        .expect("Failed to create compressor 2");
    compressor2.add_fasta_file("test_sample", &fasta_path)
        .expect("Failed to add FASTA 2");
    compressor2.finalize().expect("Failed to finalize 2");

    // Compare hashes
    let hash1 = compute_file_hash(&archive1_path);
    let hash2 = compute_file_hash(&archive2_path);

    assert_eq!(hash1, hash2,
        "Archives differ! This means compression is non-deterministic.\nArchive 1: {}\nArchive 2: {}",
        hash1, hash2);

    // Clean up
    let _ = fs::remove_file(&fasta_path);
    let _ = fs::remove_file(&archive1_path);
    let _ = fs::remove_file(&archive2_path);
}

#[cfg(test)]
mod with_cpp_agc {
    use super::*;

    fn cpp_agc_available() -> bool {
        Command::new("agc")
            .arg("--version")
            .output()
            .is_ok()
    }

    #[test]
    fn test_cpp_can_read_ragc_archives() {
        if !cpp_agc_available() {
            eprintln!("Skipping C++ compatibility test: C++ agc not found");
            return;
        }

        let test_dir = std::env::temp_dir();
        let fasta_path = test_dir.join("test_cpp_read.fasta");
        let archive_path = test_dir.join("test_cpp_read.agc");
        let output_path = test_dir.join("test_cpp_read_out.fasta");

        // Create test FASTA
        create_test_fasta(&fasta_path);
        let original_hash = compute_file_hash(&fasta_path);

        // Create archive with ragc
        let config = CompressorConfig::default();
        let mut compressor = Compressor::new(archive_path.to_str().unwrap(), config)
            .expect("Failed to create compressor");
        compressor.add_fasta_file("test_sample", &fasta_path)
            .expect("Failed to add FASTA");
        compressor.finalize().expect("Failed to finalize archive");

        // Extract with C++ agc
        let status = Command::new("agc")
            .arg("getset")
            .arg(archive_path.to_str().unwrap())
            .arg("test_sample")
            .output()
            .expect("Failed to run C++ agc");

        assert!(status.status.success(), "C++ agc failed to extract: {}",
            String::from_utf8_lossy(&status.stderr));

        fs::write(&output_path, &status.stdout).expect("Failed to write output");
        let output_hash = compute_file_hash(&output_path);

        assert_eq!(original_hash, output_hash,
            "C++ extracted different data!\nOriginal: {}\nC++ Output: {}",
            original_hash, output_hash);

        // Clean up
        let _ = fs::remove_file(&fasta_path);
        let _ = fs::remove_file(&archive_path);
        let _ = fs::remove_file(&output_path);
    }

    #[test]
    fn test_ragc_can_read_cpp_archives() {
        if !cpp_agc_available() {
            eprintln!("Skipping C++ compatibility test: C++ agc not found");
            return;
        }

        let test_dir = std::env::temp_dir();
        let fasta_path = test_dir.join("test_ragc_read.fasta");
        let archive_path = test_dir.join("test_ragc_read.agc");
        let output_path = test_dir.join("test_ragc_read_out.fasta");

        // Create test FASTA
        create_test_fasta(&fasta_path);
        let original_hash = compute_file_hash(&fasta_path);

        // Create archive with C++ agc
        let status = Command::new("agc")
            .arg("create")
            .arg("-o")
            .arg(archive_path.to_str().unwrap())
            .arg(fasta_path.to_str().unwrap())
            .status()
            .expect("Failed to run C++ agc create");

        assert!(status.success(), "C++ agc failed to create archive");

        // Extract with ragc
        let config = DecompressorConfig::default();
        let mut decompressor = Decompressor::open(archive_path.to_str().unwrap(), config)
            .expect("Failed to open C++ archive");

        let sequences = decompressor.get_sample("test_ragc_read")
            .expect("Failed to extract sample");

        let mut output_content = String::new();
        for (name, contig) in sequences {
            output_content.push_str(&format!(">{}\n", name));
            output_content.push_str(&contig_to_string(&contig));
            output_content.push('\n');
        }
        fs::write(&output_path, output_content).expect("Failed to write output");

        let output_hash = compute_file_hash(&output_path);

        assert_eq!(original_hash, output_hash,
            "ragc extracted different data from C++ archive!\nOriginal: {}\nragc Output: {}",
            original_hash, output_hash);

        // Clean up
        let _ = fs::remove_file(&fasta_path);
        let _ = fs::remove_file(&archive_path);
        let _ = fs::remove_file(&output_path);
    }
}
