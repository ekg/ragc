//! External layout parity test against a C++ AGC archive.
//!
//! This test is ignored by default. To run:
//!   RAGC_INPUTS="/path/A.fa /path/B.fa ..." \
//!   RAGC_CPP_ARCHIVE="/path/cpp_compare.agc" \
//!   cargo test -p ragc-core --test test_layout_parity_external -- --ignored --nocapture
//!
//! Notes:
//! - Uses the in-tree streaming compressor with C++ FFI cost enabled by default.
//! - Creates a temporary RAGC archive from the provided inputs (single-threaded, k=21, s=10000, m=20).
//! - Compares per-(sample,contig,index) raw segment lengths to the provided C++ archive.
//! - Group IDs may differ and are ignored; only lengths must match exactly.

use anyhow::{bail, Context, Result};
use ragc_core::{
    contig_iterator::ContigIterator,
    determine_splitters_streaming, Decompressor, DecompressorConfig, MultiFileIterator,
    StreamingQueueCompressor, StreamingQueueConfig,
};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn build_ragc_archive(out_path: &str, inputs: &[PathBuf]) -> Result<()> {
    if inputs.is_empty() { bail!("RAGC_INPUTS is empty"); }

    let k = 21usize;
    let segment_size = 10000usize;
    let min_match_len = 20usize;

    // Determine splitters from the first input (reference sample), streaming
    let (splitters, _, _) = determine_splitters_streaming(&inputs[0], k, segment_size)
        .with_context(|| format!("Failed to determine splitters from {:?}", inputs[0]))?;

    // Configure single-threaded, quiet
    let config = StreamingQueueConfig {
        k,
        segment_size,
        min_match_len,
        num_threads: 1,
        verbosity: 0,
        ..StreamingQueueConfig::default()
    };

    // Create compressor
    let mut compressor = StreamingQueueCompressor::with_splitters(out_path, config, splitters)
        .context("Failed to create streaming compressor")?;

    // Push all contigs from all inputs in the given order
    let mut it = MultiFileIterator::new(inputs.to_vec())?;
    while let Some((sample, contig, data)) = it.next_contig()? {
        compressor.push(sample, contig, data)?;
    }

    compressor.finalize()?;
    Ok(())
}

fn layout_map(archive_path: &str) -> Result<BTreeMap<(String, String), Vec<usize>>> {
    let mut dec = Decompressor::open(archive_path, DecompressorConfig { verbosity: 0 })
        .with_context(|| format!("Failed to open archive: {}", archive_path))?;
    let mut map: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (sample, contig, segs) in dec.get_all_segments()? {
        let lens = segs.iter().map(|s| s.raw_length as usize).collect::<Vec<_>>();
        map.insert((sample, contig), lens);
    }
    Ok(map)
}

#[test]
#[ignore]
fn layout_parity_against_cpp_archive() -> Result<()> {
    // Read environment config
    let cpp_archive = std::env::var("RAGC_CPP_ARCHIVE")
        .unwrap_or_else(|_| String::new());
    let inputs = std::env::var("RAGC_INPUTS")
        .unwrap_or_else(|_| String::new());

    if cpp_archive.is_empty() || inputs.is_empty() {
        eprintln!(
            "Skipping: set RAGC_CPP_ARCHIVE and RAGC_INPUTS to run this test"
        );
        return Ok(());
    }

    // Parse inputs
    let input_paths: Vec<PathBuf> = inputs
        .split_whitespace()
        .map(|s| PathBuf::from(s))
        .collect();

    // Build temp ragc archive
    let tmp_dir = tempfile::tempdir()?;
    let out_path = tmp_dir.path().join("ragc_test.agc");
    build_ragc_archive(out_path.to_str().unwrap(), &input_paths)?;

    // Load layouts
    let a = layout_map(out_path.to_str().unwrap())?; // RAGC
    let b = layout_map(&cpp_archive)?; // C++

    // Compare keys
    let keys_a: Vec<_> = a.keys().cloned().collect();
    let keys_b: Vec<_> = b.keys().cloned().collect();
    assert_eq!(keys_a, keys_b, "Sample/contig sets differ");

    // Compare lengths per segment index
    for (key, lens_a) in &a {
        let lens_b = &b[key];
        assert_eq!(lens_a.len(), lens_b.len(), "Row count differs for {:?}", key);
        for (i, (la, lb)) in lens_a.iter().zip(lens_b.iter()).enumerate() {
            assert_eq!(la, lb, "Length mismatch at {:?} index {}: {} vs {}", key, i, la, lb);
        }
    }

    Ok(())
}

