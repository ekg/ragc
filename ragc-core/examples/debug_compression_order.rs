// Debug program to check contig processing order during compression

use anyhow::Result;
use ragc_core::{
    contig_iterator::{MultiFileIterator, PansnFileIterator},
    StreamingCompressor, StreamingCompressorConfig,
};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

fn main() -> Result<()> {
    let test_dir = Path::new("/home/erik/scrapy/yeast10_test");
    let pansn_file = test_dir.join("yeast10_pansn.fa");

    // Find all sample files
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
    fasta_files.sort();

    let config = StreamingCompressorConfig {
        kmer_length: 21,
        segment_size: 10000,
        min_match_len: 20,
        compression_level: 11,
        verbosity: 0,
        group_flush_threshold: 0,
        periodic_flush_interval: 0,
        num_threads: 1, // Use single thread to ensure determinism
        adaptive_mode: false,
    };

    println!("=== Compressing with PansnFileIterator (single thread) ===");
    let pansn_archive = NamedTempFile::new()?;
    {
        let mut compressor =
            StreamingCompressor::new(pansn_archive.path().to_str().unwrap(), config.clone())?;
        let iterator = Box::new(PansnFileIterator::new(&pansn_file)?);
        compressor.add_contigs_with_splitters(iterator)?;
        compressor.finalize()?;
    }

    println!("\n=== Compressing with MultiFileIterator (single thread) ===");
    let multi_archive = NamedTempFile::new()?;
    {
        let mut compressor =
            StreamingCompressor::new(multi_archive.path().to_str().unwrap(), config.clone())?;
        let iterator = Box::new(MultiFileIterator::new(fasta_files)?);
        compressor.add_contigs_with_splitters(iterator)?;
        compressor.finalize()?;
    }

    println!("\n=== Comparing archive sizes ===");
    let pansn_size = std::fs::metadata(pansn_archive.path())?.len();
    let multi_size = std::fs::metadata(multi_archive.path())?.len();

    println!("Pansn: {pansn_size} bytes");
    println!("Multi: {multi_size} bytes");

    if pansn_size == multi_size {
        println!("✓ Sizes match!");
    } else {
        println!(
            "❌ Sizes differ by {} bytes",
            (pansn_size as i64 - multi_size as i64).abs()
        );
    }

    Ok(())
}
