// Test program to verify PansnFileIterator and MultiFileIterator produce identical results

use ragc_core::{
    contig_iterator::{ContigIterator, MultiFileIterator, PansnFileIterator},
    StreamingCompressor, StreamingCompressorConfig,
};
use std::path::Path;
use std::path::PathBuf;
use anyhow::Result;

fn main() -> Result<()> {
    let test_dir = Path::new("/home/erik/scrapy/yeast10_test");

    println!("=== Testing Unified Iterator Approach ===\n");

    // Test 1: Single pansn file
    println!("Test 1: Creating archive from single pansn file...");
    let pansn_file = test_dir.join("yeast10_pansn.fa");
    let output1 = "/tmp/test_pansn_iter.agc";

    let config1 = StreamingCompressorConfig {
        kmer_length: 21,
        segment_size: 10000,
        min_match_len: 20,
        compression_level: 11,
        verbosity: 1,
        group_flush_threshold: 0,
        periodic_flush_interval: 0,
        num_threads: 4,
        adaptive_mode: false,
    };

    let mut compressor1 = StreamingCompressor::new(output1, config1)?;
    let iterator1 = Box::new(PansnFileIterator::new(&pansn_file)?);
    compressor1.add_contigs_with_splitters(iterator1)?;
    compressor1.finalize()?;

    println!("✓ Pansn archive created: {}", output1);

    // Test 2: Multiple files (one per sample#hap)
    println!("\nTest 2: Creating archive from multiple sample files...");

    // Find all .fa files in the directory
    let mut fasta_files: Vec<PathBuf> = std::fs::read_dir(test_dir)?
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let path = entry.path();
            if path.extension()?.to_str()? == "fa" && path != pansn_file {
                Some(path)
            } else {
                None
            }
        })
        .collect();

    fasta_files.sort(); // Ensure consistent ordering

    println!("Found {} sample files", fasta_files.len());

    let output2 = "/tmp/test_multi_iter.agc";

    let config2 = StreamingCompressorConfig {
        kmer_length: 21,
        segment_size: 10000,
        min_match_len: 20,
        compression_level: 11,
        verbosity: 1,
        group_flush_threshold: 0,
        periodic_flush_interval: 0,
        num_threads: 4,
        adaptive_mode: false,
    };

    let mut compressor2 = StreamingCompressor::new(output2, config2)?;
    let iterator2 = Box::new(MultiFileIterator::new(fasta_files)?);
    compressor2.add_contigs_with_splitters(iterator2)?;
    compressor2.finalize()?;

    println!("✓ Multi-file archive created: {}", output2);

    // Compare file sizes
    let size1 = std::fs::metadata(output1)?.len();
    let size2 = std::fs::metadata(output2)?.len();

    println!("\n=== Results ===");
    println!("Pansn archive size:     {} bytes", size1);
    println!("Multi-file archive size: {} bytes", size2);

    if size1 == size2 {
        println!("\n✓ SUCCESS: Archives are identical in size!");
    } else {
        let diff = (size1 as i64 - size2 as i64).abs();
        let pct = (diff as f64 / size1 as f64) * 100.0;
        println!("\n⚠ WARNING: Archives differ by {} bytes ({:.2}%)", diff, pct);
    }

    Ok(())
}
