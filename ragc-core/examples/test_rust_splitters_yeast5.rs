// Test Rust splitter integration with yeast5_chrV.fa (5 samples, PanSN format)
// Uses actual multi-sample dataset to verify byte-identical output

use ragc_core::agc_compress_ffi::{compress_with_cpp_agc, compress_with_rust_splitters};
use std::process::Command;
use std::fs;

fn main() -> anyhow::Result<()> {
    println!("=== Testing Rust Splitters with yeast5 chrV ===\n");

    let test_dir = "/tmp/ragc_rust_splitter_yeast5";
    fs::create_dir_all(test_dir)?;

    // Use existing test file (PanSN format: 5 samples in one file)
    let test_file = "./test_minimal/yeast5_chrV.fa";

    if !std::path::Path::new(test_file).exists() {
        eprintln!("Error: Test file not found: {}", test_file);
        return Err(anyhow::anyhow!("Missing test file"));
    }

    let sample_files = vec![
        ("yeast5".to_string(), test_file.to_string()),
    ];

    // Parameters from CLAUDE.md (matching verified tests)
    let kmer_length = 21;
    let segment_size = 10000;
    let min_match_length = 20;
    let pack_cardinality = 50;
    let concatenated_genomes = true;  // PanSN format
    let adaptive_compression = false;
    let verbosity = 1;  // Reduce noise
    let no_threads = 1;  // Single-threaded for determinism
    let fallback_frac = 0.0;

    println!("Input: {} ({:.1} MB)",
             test_file,
             fs::metadata(test_file)?.len() as f64 / (1024.0 * 1024.0));

    // Test 1: Native C++ AGC (baseline)
    println!("\n1. Creating archive with native C++ AGC splitters...");
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
    println!("   Size: {:.1} KB", cpp_size as f64 / 1024.0);
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
    println!("   Size: {:.1} KB", rust_size as f64 / 1024.0);
    println!("   SHA256: {}", rust_sha);

    // Compare results
    println!("\n=== Comparison ===");
    let time_diff = rust_time.as_secs_f64() - cpp_time.as_secs_f64();
    let size_diff = rust_size as i64 - cpp_size as i64;

    println!("Time difference: {:+.2}s ({:+.1}%)",
             time_diff,
             (time_diff / cpp_time.as_secs_f64()) * 100.0);
    println!("Size difference: {:+} bytes", size_diff);

    if cpp_sha == rust_sha {
        println!("\n✅ SUCCESS: Archives are BYTE-IDENTICAL!");
        println!("   ✓ Splitter detection replaced with Rust");
        println!("   ✓ Compression output unchanged");
        println!("   ✓ First component replacement verified!");
        println!("\n   Ready to replace next component.");
    } else {
        println!("\n❌ FAILURE: Archives differ");
        println!("   Native C++:     {}", cpp_sha);
        println!("   Rust splitters: {}", rust_sha);

        if size_diff != 0 {
            println!("   Size diff: {:+.1}%",
                     (size_diff as f64 / cpp_size as f64) * 100.0);
        }

        // Show first differences
        println!("\n   Binary diff (first 10 bytes):");
        let diff_result = Command::new("cmp")
            .args(&["-l", &cpp_archive, &rust_archive])
            .output()?;

        if !diff_result.stdout.is_empty() {
            let diff_output = std::str::from_utf8(&diff_result.stdout)?;
            for (i, line) in diff_output.lines().take(10).enumerate() {
                println!("      {} {}", i+1, line);
            }

            let total_diffs = diff_output.lines().count();
            if total_diffs > 10 {
                println!("      ... ({} more differences)", total_diffs - 10);
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
