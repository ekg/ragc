// Integration test: Verify PansnFileIterator and MultiFileIterator produce equivalent archives
//
// This test ensures both input methods (single pansn file vs multiple sample files)
// produce archives with identical decompressed content.

use ragc_core::{
    contig_iterator::{MultiFileIterator, PansnFileIterator},
    Decompressor, DecompressorConfig, StreamingCompressor, StreamingCompressorConfig,
};
use std::collections::HashMap;
use std::path::PathBuf;
use tempfile::NamedTempFile;

#[test]
fn test_unified_iterator_equivalence() {
    // Skip if test data doesn't exist
    let test_dir = std::path::Path::new("/home/erik/scrapy/yeast10_test");
    if !test_dir.exists() {
        eprintln!("Skipping test: test data not found at {test_dir:?}");
        return;
    }

    let pansn_file = test_dir.join("yeast10_pansn.fa");
    if !pansn_file.exists() {
        eprintln!("Skipping test: pansn file not found");
        return;
    }

    // Find all sample files
    let mut fasta_files: Vec<PathBuf> = std::fs::read_dir(test_dir)
        .unwrap()
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

    if fasta_files.is_empty() {
        eprintln!("Skipping test: no sample files found");
        return;
    }

    fasta_files.sort();

    // Test config (matching the example)
    let config = StreamingCompressorConfig {
        kmer_length: 21,
        segment_size: 10000,
        min_match_len: 20,
        compression_level: 11,
        verbosity: 0, // Quiet for tests
        group_flush_threshold: 0,
        periodic_flush_interval: 0,
        num_threads: 1, // FIXME: Multi-threading has race conditions, use single thread for now
        adaptive_mode: false,
    };

    // Create temporary output files
    let pansn_archive = NamedTempFile::new().unwrap();
    let multi_archive = NamedTempFile::new().unwrap();

    // Compress using PansnFileIterator
    {
        let mut compressor =
            StreamingCompressor::new(pansn_archive.path().to_str().unwrap(), config.clone())
                .unwrap();
        let iterator = Box::new(PansnFileIterator::new(&pansn_file).unwrap());
        compressor.add_contigs_with_splitters(iterator).unwrap();
        compressor.finalize().unwrap();
    }

    // Compress using MultiFileIterator
    {
        let mut compressor =
            StreamingCompressor::new(multi_archive.path().to_str().unwrap(), config.clone())
                .unwrap();
        let iterator = Box::new(MultiFileIterator::new(fasta_files).unwrap());
        compressor.add_contigs_with_splitters(iterator).unwrap();
        compressor.finalize().unwrap();
    }

    // Verify archives have similar sizes (within 1%)
    let pansn_size = std::fs::metadata(pansn_archive.path()).unwrap().len();
    let multi_size = std::fs::metadata(multi_archive.path()).unwrap().len();

    let size_diff = (pansn_size as i64 - multi_size as i64).abs() as f64;
    let size_diff_pct = (size_diff / pansn_size as f64) * 100.0;

    assert!(
        size_diff_pct < 1.0,
        "Archive sizes differ by more than 1%: pansn={pansn_size} multi={multi_size} diff={size_diff_pct:.2}%"
    );

    // Decompress and verify content is identical (order-independent)
    let decompressor_config = DecompressorConfig::default();

    let pansn_samples = {
        let mut decompressor = Decompressor::open(
            pansn_archive.path().to_str().unwrap(),
            decompressor_config.clone(),
        )
        .unwrap();
        extract_all_samples(&mut decompressor)
    };

    let multi_samples = {
        let mut decompressor =
            Decompressor::open(multi_archive.path().to_str().unwrap(), decompressor_config)
                .unwrap();
        extract_all_samples(&mut decompressor)
    };

    // Verify same samples
    assert_eq!(
        pansn_samples.len(),
        multi_samples.len(),
        "Different number of samples"
    );

    for (sample_name, pansn_contigs) in &pansn_samples {
        let multi_contigs = multi_samples
            .get(sample_name)
            .unwrap_or_else(|| panic!("Sample {sample_name} missing from multi-file archive"));

        assert_eq!(
            pansn_contigs.len(),
            multi_contigs.len(),
            "Sample {sample_name} has different number of contigs"
        );

        // Compare contigs (order-independent)
        for (contig_name, pansn_seq) in pansn_contigs {
            let multi_seq = multi_contigs.get(contig_name).unwrap_or_else(|| {
                panic!(
                    "Contig {contig_name} missing from sample {sample_name} in multi-file archive"
                )
            });

            assert_eq!(
                pansn_seq, multi_seq,
                "Contig {contig_name} in sample {sample_name} has different sequence"
            );
        }
    }

    println!("✓ Both iterators produce equivalent archives");
    println!("  Pansn size: {pansn_size} bytes");
    println!("  Multi size: {multi_size} bytes ({size_diff_pct:.2}% difference)");
    println!("  Samples: {}", pansn_samples.len());
    println!(
        "  Total contigs: {}",
        pansn_samples.values().map(|c| c.len()).sum::<usize>()
    );
}

/// Extract all samples and contigs from an archive
fn extract_all_samples(
    decompressor: &mut Decompressor,
) -> HashMap<String, HashMap<String, Vec<u8>>> {
    let mut result = HashMap::new();

    for sample_name in decompressor.list_samples() {
        let contigs = decompressor.get_sample(&sample_name).unwrap();
        let contig_map: HashMap<String, Vec<u8>> = contigs.into_iter().collect();
        result.insert(sample_name, contig_map);
    }

    result
}
