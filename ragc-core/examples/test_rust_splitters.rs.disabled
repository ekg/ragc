// Test Rust splitter integration with C++ AGC compression
// Verifies that using Rust-computed splitters produces byte-identical archives

use ragc_core::agc_compress_ffi::{compress_with_cpp_agc, compress_with_rust_splitters};
use std::process::Command;
use std::fs;

fn main() -> anyhow::Result<()> {
    println!("=== Testing Rust Splitter Integration ===\n");

    // Test files from previous verification
    let test_dir = "/tmp/ragc_rust_splitter_test";
    fs::create_dir_all(test_dir)?;

    let sample_files = vec![
        ("sample1".to_string(), "/tmp/test_sample.fa".to_string()),
    ];

    // Ensure test file exists
    if !std::path::Path::new(&sample_files[0].1).exists() {
        eprintln!("Creating test file: {}", sample_files[0].1);
        fs::write(&sample_files[0].1,
            ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n")?;
    }

    // Compression parameters (matching C++ AGC defaults)
    let kmer_length = 21;
    let segment_size = 10000;
    let min_match_length = 20;
    let pack_cardinality = 50;
    let concatenated_genomes = false;
    let adaptive_compression = false;
    let verbosity = 2;
    let no_threads = 1;
    let fallback_frac = 0.0;

    // Test 1: Native C++ AGC (baseline)
    println!("1. Creating archive with native C++ AGC...");
    let cpp_archive = format!("{}/native_cpp.agc", test_dir);
    compress_with_cpp_agc(
        &cpp_archive,
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

    let cpp_size = fs::metadata(&cpp_archive)?.len();
    let cpp_sha = get_sha256(&cpp_archive)?;
    println!("   Size: {} bytes", cpp_size);
    println!("   SHA256: {}", cpp_sha);

    // Test 2: C++ AGC with Rust-computed splitters
    println!("\n2. Creating archive with Rust splitters...");
    let rust_archive = format!("{}/rust_splitters.agc", test_dir);
    compress_with_rust_splitters(
        &rust_archive,
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

    let rust_size = fs::metadata(&rust_archive)?.len();
    let rust_sha = get_sha256(&rust_archive)?;
    println!("   Size: {} bytes", rust_size);
    println!("   SHA256: {}", rust_sha);

    // Compare results
    println!("\n=== Comparison ===");
    println!("Size difference: {} bytes", (rust_size as i64 - cpp_size as i64).abs());

    if cpp_sha == rust_sha {
        println!("✅ SUCCESS: Archives are BYTE-IDENTICAL!");
        println!("   Rust splitter replacement produces identical output to C++ AGC");
    } else {
        println!("⚠️  Archives differ:");
        println!("   Native C++:     {}", cpp_sha);
        println!("   Rust splitters: {}", rust_sha);
        println!("   Size diff: {} bytes", rust_size as i64 - cpp_size as i64);

        // Binary diff
        let diff_result = Command::new("cmp")
            .args(&["-l", &cpp_archive, &rust_archive])
            .output()?;

        if !diff_result.stdout.is_empty() {
            let diff_lines: Vec<&str> = std::str::from_utf8(&diff_result.stdout)?
                .lines()
                .take(10)
                .collect();
            println!("\n   First 10 byte differences:");
            for line in diff_lines {
                println!("      {}", line);
            }
        }
    }

    println!("\n=== Test archives created in {} ===", test_dir);
    Ok(())
}

fn get_sha256(path: &str) -> anyhow::Result<String> {
    let output = Command::new("sha256sum")
        .arg(path)
        .output()?;

    let stdout = String::from_utf8(output.stdout)?;
    let hash = stdout.split_whitespace().next()
        .ok_or_else(|| anyhow::anyhow!("Failed to parse sha256sum output"))?;

    Ok(hash.to_string())
}
