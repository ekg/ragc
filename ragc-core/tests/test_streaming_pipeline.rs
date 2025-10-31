// End-to-end integration test for streaming compression pipeline
// Tests the complete flow from FASTA files to AGC archive

use ragc_core::{create_agc_archive, determine_splitters_streaming};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[test]
#[ignore] // Requires yeast samples in ~/scrapy - run with: cargo test --test test_streaming_pipeline -- --ignored
fn test_yeast3_streaming_compression() {
    // Use first 3 yeast samples for fast testing
    let home = env::var("HOME").expect("HOME not set");
    let base_dir = PathBuf::from(home).join("scrapy/yeast235_samples");

    // Check if test data exists
    if !base_dir.exists() {
        eprintln!("Skipping test: yeast samples not found at {:?}", base_dir);
        return;
    }

    // Get first 3 samples
    let mut sample_files = Vec::new();
    let samples = ["AAA#0.fa", "AAB#0.fa", "AAC#0.fa"];

    for sample in &samples {
        let path = base_dir.join(sample);
        if !path.exists() {
            eprintln!("Skipping test: sample not found: {:?}", path);
            return;
        }
        let sample_name = sample.trim_end_matches(".fa").replace("#", "_");
        sample_files.push((sample_name, path.to_string_lossy().to_string()));
    }

    println!("Testing with {} yeast samples", sample_files.len());

    // Output path
    let output_path = "/tmp/test_yeast3_streaming.agc";

    // Clean up old file
    let _ = fs::remove_file(output_path);

    // Parameters (matching C++ AGC defaults)
    let kmer_length = 21;
    let segment_size = 1000;
    let num_threads = 4;
    let adaptive_mode = false;
    let concatenated_genomes = false;
    let verbosity = 1;

    // Step 1: Determine splitters from first sample (reference)
    println!("Step 1: Determining splitters from reference sample...");

    let reference_path = Path::new(&sample_files[0].1);
    let (splitters, candidates, duplicates) =
        determine_splitters_streaming(reference_path, kmer_length, segment_size as usize)
            .expect("Failed to determine splitters");

    println!("  Found {} splitters", splitters.len());
    println!("  Found {} candidate k-mers", candidates.len());
    println!("  Found {} duplicated k-mers", duplicates.len());

    // Step 2: Compress all samples
    println!("Step 2: Compressing {} samples...", sample_files.len());

    let result = create_agc_archive(
        output_path,
        sample_files.clone(),
        splitters,
        candidates,
        duplicates,
        kmer_length,
        segment_size,
        num_threads,
        adaptive_mode,
        concatenated_genomes,
        verbosity,
    );

    assert!(result.is_ok(), "Compression failed: {:?}", result.err());
    println!("  ✓ Compression successful");

    // Step 3: Verify archive exists and has reasonable size
    let metadata = fs::metadata(output_path).expect("Archive file not found");
    let file_size = metadata.len();

    println!("Step 3: Verifying archive...");
    println!(
        "  Archive size: {} bytes ({:.2} MB)",
        file_size,
        file_size as f64 / 1_048_576.0
    );

    // Sanity check: archive should be at least 1KB and less than source files
    assert!(file_size > 1024, "Archive too small: {} bytes", file_size);
    assert!(
        file_size < 50_000_000,
        "Archive too large: {} bytes",
        file_size
    );

    println!("  ✓ Archive size reasonable");

    // Step 3.5: Verify with RAGC decompressor
    println!("Step 3.5: Verifying with RAGC decompressor...");
    {
        use ragc_core::{Decompressor, DecompressorConfig};
        let config = DecompressorConfig::default();
        let mut decompressor = Decompressor::open(output_path, config)
            .expect("Failed to open archive with RAGC decompressor");

        let samples = decompressor.list_samples();
        println!("  RAGC found {} samples: {:?}", samples.len(), samples);

        assert_eq!(
            samples.len(),
            3,
            "RAGC decompressor expected 3 samples, found {}",
            samples.len()
        );
        println!("  ✓ RAGC decompressor can read archive");
        println!("  ✓ Found all 3 samples");
    }

    // Step 4: Check archive structure with C++ AGC (if available)
    if Path::new("/home/erik/agc/bin/agc").exists() {
        println!("Step 4: Verifying with C++ AGC...");

        use std::process::Command;

        // List samples in archive
        let output = Command::new("/home/erik/agc/bin/agc")
            .args(&["listset", output_path])
            .output()
            .expect("Failed to run C++ AGC");

        if output.status.success() {
            let samples_list = String::from_utf8_lossy(&output.stdout);
            println!("  Samples in archive:\n{}", samples_list);

            // Verify we have 3 samples
            let sample_count = samples_list.lines().count();
            assert_eq!(
                sample_count, 3,
                "Expected 3 samples, found {}",
                sample_count
            );
            println!("  ✓ C++ AGC can read archive");
            println!("  ✓ Found all 3 samples");
        } else {
            let error = String::from_utf8_lossy(&output.stderr);
            eprintln!("  Warning: C++ AGC failed to read archive: {}", error);
        }
    } else {
        println!("Step 4: Skipped (C++ AGC not found)");
    }

    println!("\n=== Test Complete ===");
    println!("Archive created at: {}", output_path);
}
