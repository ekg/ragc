// Full round-trip test: FASTA → AGC → FASTA
// Verifies that we can compress and decompress without data loss

use ragc_core::{Compressor, CompressorConfig, Decompressor, DecompressorConfig};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
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

/// Helper function to extract sequences from FASTA file, keyed by contig name
/// Returns a BTreeMap for deterministic ordering
fn extract_sequences_with_checksums(fasta_path: &Path) -> BTreeMap<String, String> {
    let content = fs::read_to_string(fasta_path).unwrap();
    let mut sequences = BTreeMap::new();
    let mut current_contig = None;
    let mut current_seq = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('>') {
            // Save previous contig
            if let Some(contig) = current_contig.take() {
                sequences.insert(contig, current_seq.clone());
                current_seq.clear();
            }
            // Start new contig - remove '>' and any sample prefix (e.g., "sample#0#chr1" -> "chr1")
            let header = &line[1..];
            let contig_name = header.split('#').last().unwrap_or(header).to_string();
            current_contig = Some(contig_name);
        } else {
            current_seq.push_str(line);
        }
    }

    // Save last contig
    if let Some(contig) = current_contig {
        sequences.insert(contig, current_seq);
    }

    sequences
}

/// Compute SHA256 checksum of all sequences concatenated in sorted order
fn compute_sequence_checksum(sequences: &BTreeMap<String, String>) -> String {
    let mut hasher = Sha256::new();

    // Concatenate all sequences in sorted order by contig name
    for (_contig, seq) in sequences.iter() {
        hasher.update(seq.as_bytes());
    }

    format!("{:x}", hasher.finalize())
}

#[test]
fn test_roundtrip_with_checksum_verification() {
    let test_fasta = Path::new("../test-data/test_simple.fasta");

    // Skip if test file doesn't exist
    if !test_fasta.exists() {
        eprintln!("Skipping test_roundtrip_with_checksum_verification: test file not found");
        return;
    }

    let archive_path = "/tmp/test_roundtrip_checksum.agc";
    let output_fasta = "/tmp/test_roundtrip_checksum_output.fasta";
    let _ = fs::remove_file(archive_path);
    let _ = fs::remove_file(output_fasta);

    // Extract original sequences and compute checksum
    let original_sequences = extract_sequences_with_checksums(test_fasta);
    let original_checksum = compute_sequence_checksum(&original_sequences);

    eprintln!("Original sequences:");
    for (contig, seq) in &original_sequences {
        eprintln!("  {}: {} bases", contig, seq.len());
    }
    eprintln!("Original checksum: {}", original_checksum);

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

    // Extract decompressed sequences and compute checksum
    let decompressed_sequences = extract_sequences_with_checksums(Path::new(output_fasta));
    let decompressed_checksum = compute_sequence_checksum(&decompressed_sequences);

    eprintln!("Decompressed sequences:");
    for (contig, seq) in &decompressed_sequences {
        eprintln!("  {}: {} bases", contig, seq.len());
    }
    eprintln!("Decompressed checksum: {}", decompressed_checksum);

    // Verify checksums match
    assert_eq!(
        original_checksum, decompressed_checksum,
        "SHA256 checksums must match for round-trip compression/decompression"
    );

    // Verify we have the same contigs
    assert_eq!(
        original_sequences.len(),
        decompressed_sequences.len(),
        "Number of contigs should match"
    );

    // Verify each contig individually
    for (contig, original_seq) in &original_sequences {
        let decompressed_seq = decompressed_sequences
            .get(contig)
            .expect(&format!("Contig {} missing in decompressed output", contig));

        assert_eq!(
            original_seq.to_uppercase(),
            decompressed_seq.to_uppercase(),
            "Sequence for contig {} should match exactly",
            contig
        );
    }

    eprintln!("✓ Round-trip test passed with cryptographic checksum verification!");

    fs::remove_file(archive_path).unwrap();
    fs::remove_file(output_fasta).unwrap();
}
