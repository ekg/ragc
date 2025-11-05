#![allow(clippy::all)]
// Minimal test to reproduce StreamingCompressor data corruption
use ragc_core::{Decompressor, DecompressorConfig};
use ragc_core::{StreamingCompressor, StreamingCompressorConfig};
use sha2::{Digest, Sha256};
use std::fs;

#[test]
fn test_streaming_compressor_extraction_correctness() {
    let archive_path = "/tmp/test_streaming_corruption.agc";
    let _ = fs::remove_file(archive_path);

    // Use yeast10 first sample
    let input_path = "/home/erik/scrapy/yeast10_test/AEL#2.fa";

    if !std::path::Path::new(input_path).exists() {
        eprintln!("Skipping test: {} not found", input_path);
        return;
    }

    // Read original sequence (normalize to uppercase for comparison)
    let original_content = fs::read_to_string(input_path).unwrap();
    let original_sequence: String = original_content
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<Vec<_>>()
        .join("")
        .to_uppercase()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect();

    let mut original_hasher = Sha256::new();
    original_hasher.update(original_sequence.as_bytes());
    let original_hash = format!("{:x}", original_hasher.finalize());

    println!("Original hash: {}", original_hash);

    // Compress with StreamingCompressor (the CLI path)
    {
        let config = StreamingCompressorConfig {
            kmer_length: 21,
            segment_size: 10000,
            min_match_len: 20,
            compression_level: 11,
            verbosity: 0, // Disable verbosity for cleaner output
            group_flush_threshold: 0,
            concatenated_genomes: false,
            periodic_flush_interval: 0,
            num_threads: 15, // Test multi-threaded determinism
            adaptive_mode: false,
        };

        let mut compressor = StreamingCompressor::new(archive_path, config).unwrap();
        compressor
            .add_fasta_files_with_splitters(&[(
                "AEL#2".to_string(),
                std::path::Path::new(input_path),
            )])
            .unwrap();
        compressor.finalize().unwrap();
    }

    // Extract and verify
    {
        let config = DecompressorConfig { verbosity: 0 };
        let mut decompressor = Decompressor::open(archive_path, config).unwrap();

        // Get the sample
        let samples = decompressor.list_samples();
        println!("Samples in archive: {:?}", samples);
        assert!(!samples.is_empty(), "No samples in archive");

        let sample_name = &samples[0];
        let contigs = decompressor.get_sample(sample_name).unwrap();

        // Extract sequence from result
        let extracted_sequence: String = contigs
            .iter()
            .flat_map(|(_, contig)| {
                contig.iter().map(|&b| match b {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                })
            })
            .collect();

        let mut extracted_hasher = Sha256::new();
        extracted_hasher.update(extracted_sequence.as_bytes());
        let extracted_hash = format!("{:x}", extracted_hasher.finalize());

        println!("Extracted hash: {}", extracted_hash);
        println!("Original length: {}", original_sequence.len());
        println!("Extracted length: {}", extracted_sequence.len());

        // Find where sequences diverge
        if original_hash != extracted_hash {
            let mut first_diff: Option<usize> = None;
            let mut diff_count = 0;

            for (i, (o, e)) in original_sequence
                .chars()
                .zip(extracted_sequence.chars())
                .enumerate()
            {
                if o != e {
                    if first_diff.is_none() {
                        first_diff = Some(i);
                        println!("\nFirst difference at position {}", i);
                        let start = i.saturating_sub(50);
                        let end = (i + 50).min(original_sequence.len());
                        println!(
                            "Original  [{}..{}]: {}",
                            start,
                            end,
                            &original_sequence[start..end]
                        );
                        println!(
                            "Extracted [{}..{}]: {}",
                            start,
                            end,
                            &extracted_sequence[start..end]
                        );
                    }
                    diff_count += 1;
                    if diff_count >= 100 {
                        break;
                    }
                }
            }

            println!("\nTotal differences found: {} (stopped at 100)", diff_count);
        }

        assert_eq!(
            original_hash, extracted_hash,
            "Data corruption detected! Sequences don't match.\nOriginal:  {}\nExtracted: {}",
            original_hash, extracted_hash
        );
    }
}
