// Test Rust splitter integration with realistic yeast chr5 data
// Uses actual dataset with significant splitters to verify byte-identical output

use ragc_core::agc_compress_ffi::{compress_with_cpp_agc, compress_with_rust_splitters};
use std::process::Command;
use std::fs;

fn main() -> anyhow::Result<()> {
    println!("=== Testing Rust Splitters with Yeast chr5 ===\n");

    let test_dir = "/tmp/ragc_rust_splitter_chr5";
    fs::create_dir_all(test_dir)?;

    // Use existing yeast chr5 test files
    let sample_files = vec![
        ("AAA".to_string(), "/tmp/yeast_samples/sample_AAA_chr5.fa".to_string()),
        ("AAB".to_string(), "/tmp/yeast_samples/sample_AAB_chr5.fa".to_string()),
        ("AAC".to_string(), "/tmp/yeast_samples/sample_AAC_chr5.fa".to_string()),
    ];

    // Check if test files exist
    for (name, path) in &sample_files {
        if !std::path::Path::new(path).exists() {
            eprintln!("Error: Test file not found: {}", path);
            eprintln!("Run: scripts/prepare_yeast_chr5.sh");
            return Err(anyhow::anyhow!("Missing test files"));
        }
    }

    // Parameters from previous tests (matching C++ AGC)
    let kmer_length = 21;
    let segment_size = 10000;
    let min_match_length = 20;
    let pack_cardinality = 50;
    let concatenated_genomes = false;
    let adaptive_compression = false;
    let verbosity = 1;  // Reduce noise
    let no_threads = 1;  // Single-threaded for determinism
    let fallback_frac = 0.0;

    // Test 1: Native C++ AGC (baseline)
    println!("1. Creating archive with native C++ AGC splitters...");
    let cpp_archive = format!("{}/native_cpp.agc", test_dir);
    let start = std::time::Instant::now();
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
    let cpp_time = start.elapsed();

    let cpp_size = fs::metadata(&cpp_archive)?.len();
    let cpp_sha = get_sha256(&cpp_archive)?;
    println!("   Time: {:.2}s", cpp_time.as_secs_f64());
    println!("   Size: {} bytes ({:.1} KB)", cpp_size, cpp_size as f64 / 1024.0);
    println!("   SHA256: {}", cpp_sha);

    // Test 2: C++ AGC with Rust-computed splitters
    println!("\n2. Creating archive with Rust-computed splitters...");
    let rust_archive = format!("{}/rust_splitters.agc", test_dir);
    let start = std::time::Instant::now();
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
    let rust_time = start.elapsed();

    let rust_size = fs::metadata(&rust_archive)?.len();
    let rust_sha = get_sha256(&rust_archive)?;
    println!("   Time: {:.2}s", rust_time.as_secs_f64());
    println!("   Size: {} bytes ({:.1} KB)", rust_size, rust_size as f64 / 1024.0);
    println!("   SHA256: {}", rust_sha);

    // Compare results
    println!("\n=== Comparison ===");
    println!("Time difference: {:.2}s", (rust_time.as_secs_f64() - cpp_time.as_secs_f64()).abs());
    println!("Size difference: {} bytes", (rust_size as i64 - cpp_size as i64).abs());

    if cpp_sha == rust_sha {
        println!("\n✅ SUCCESS: Archives are BYTE-IDENTICAL!");
        println!("   First component replacement (splitters) verified!");
        println!("   Ready to replace next component.");
    } else {
        println!("\n❌ FAILURE: Archives differ");
        println!("   Native C++:     {}", cpp_sha);
        println!("   Rust splitters: {}", rust_sha);
        println!("   Size diff: {} bytes ({:.1}%)",
                 rust_size as i64 - cpp_size as i64,
                 ((rust_size as f64 - cpp_size as f64) / cpp_size as f64) * 100.0);

        // Show first differences
        println!("\n   Checking binary differences...");
        let diff_result = Command::new("cmp")
            .args(&["-l", &cpp_archive, &rust_archive])
            .output()?;

        if !diff_result.stdout.is_empty() {
            let diff_lines: Vec<&str> = std::str::from_utf8(&diff_result.stdout)?
                .lines()
                .take(5)
                .collect();
            println!("   First 5 byte differences (offset, cpp, rust):");
            for line in diff_lines {
                println!("      {}", line);
            }
        }
    }

    println!("\n=== Archives saved to {} ===", test_dir);
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
