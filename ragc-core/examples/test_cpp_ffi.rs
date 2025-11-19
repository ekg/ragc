// Test program for C++ AGC FFI compression
use ragc_core::agc_compress_ffi::compress_with_cpp_agc;
use std::path::Path;

fn main() -> anyhow::Result<()> {
    println!("Testing C++ AGC FFI compression...\n");

    // Use a simple test file
    let test_file = "/home/erik/ragc/test-data/test_simple.fasta";
    if !Path::new(test_file).exists() {
        eprintln!("Test file not found: {}", test_file);
        eprintln!("Please provide a valid FASTA file path as first argument");
        std::process::exit(1);
    }

    let output = "/tmp/test_cpp_ffi.agc";
    let sample_name = "test_sample";

    println!("Input: {}", test_file);
    println!("Output: {}", output);
    println!("Sample name: {}", sample_name);
    println!();

    // Prepare sample files: (sample_name, file_path)
    let sample_files = vec![(sample_name.to_string(), test_file.to_string())];

    // Compression parameters (matching C++ AGC defaults)
    let kmer_length = 21;
    let segment_size = 10000;
    let min_match_length = 20;
    let pack_cardinality = 50;
    let concatenated_genomes = false;
    let adaptive_compression = false;
    let verbosity = 2;
    let no_threads = 1; // Single-threaded for determinism
    let fallback_frac = 0.0;

    println!("Calling C++ AGC compression with parameters:");
    println!("  k-mer length: {}", kmer_length);
    println!("  segment size: {}", segment_size);
    println!("  min match: {}", min_match_length);
    println!("  threads: {}", no_threads);
    println!();

    // Call C++ AGC compression via FFI
    compress_with_cpp_agc(
        output,
        &sample_files,
        kmer_length,
        segment_size,
        min_match_length,
        pack_cardinality,
        concatenated_genomes,
        adaptive_compression,
        verbosity,
        no_threads,
        fallback_frac,
    )?;

    println!("\n✓ Compression succeeded!");
    println!("Output archive: {}", output);

    // Check file size
    if let Ok(metadata) = std::fs::metadata(output) {
        println!("Archive size: {} bytes", metadata.len());
    }

    Ok(())
}
