// Minimal test to pinpoint memory corruption bug

use anyhow::Result;
use ragc_core::{
    contig_iterator::{ContigIterator, PansnFileIterator},
    Decompressor, DecompressorConfig, StreamingCompressor, StreamingCompressorConfig,
};
use std::path::Path;
use tempfile::NamedTempFile;

fn validate_bases(data: &[u8], label: &str) -> bool {
    for (i, &base) in data.iter().enumerate() {
        if base > 3 {
            println!("❌ {label} Invalid base at position {i}: {base} (0x{base:02X})");
            // Show context
            let start = i.saturating_sub(10);
            let end = (i + 10).min(data.len());
            println!("   Context [{}-{}]: {:?}", start, end, &data[start..end]);
            return false;
        }
    }
    true
}

fn main() -> Result<()> {
    let test_dir = Path::new("/home/erik/scrapy/yeast10_test");
    let pansn_file = test_dir.join("yeast10_pansn.fa");

    // Collect first 3 contigs from iterator with validation
    println!("=== Reading contigs from iterator ===");
    let mut iterator = PansnFileIterator::new(&pansn_file)?;
    let mut test_contigs = Vec::new();
    let mut count = 0;

    while let Some((sample_name, contig_name, sequence)) = iterator.next_contig()? {
        if !sequence.is_empty() {
            print!(
                "Contig {}: {} / {} ({} bases) ",
                count,
                sample_name,
                contig_name,
                sequence.len()
            );
            if validate_bases(&sequence, &format!("INPUT {count}")) {
                println!("✓");
            }
            test_contigs.push((sample_name.clone(), contig_name.clone(), sequence));
            count += 1;
            if count >= 3 {
                break;
            }
        }
    }

    // Create archive using add_contigs_with_splitters
    println!("\n=== Compressing with add_contigs_with_splitters ===");
    let archive = NamedTempFile::new()?;

    {
        let config = StreamingCompressorConfig {
            kmer_length: 21,
            segment_size: 10000,
            min_match_len: 20,
            compression_level: 11,
            verbosity: 2, // Verbose output
            group_flush_threshold: 0,
            concatenated_genomes: false,
            periodic_flush_interval: 0,
            num_threads: 1, // Single thread for determinism
            adaptive_mode: false,
        };

        let mut compressor = StreamingCompressor::new(archive.path().to_str().unwrap(), config)?;
        let iterator = Box::new(PansnFileIterator::new(&pansn_file)?);
        compressor.add_contigs_with_splitters(iterator)?;
        compressor.finalize()?;
    }

    // Decompress and validate
    println!("\n=== Decompressing and validating ===");
    let decompressor_config = DecompressorConfig::default();
    let mut decompressor =
        Decompressor::open(archive.path().to_str().unwrap(), decompressor_config)?;

    // Get unique sample names from test contigs
    let mut samples_to_check: Vec<String> =
        test_contigs.iter().map(|(s, _, _)| s.clone()).collect();
    samples_to_check.sort();
    samples_to_check.dedup();

    for sample_name in samples_to_check {
        println!("Sample: {sample_name}");
        let contigs = decompressor.get_sample(&sample_name)?;

        for (contig_name, decompressed) in contigs {
            // Find expected sequence
            if let Some((_, _, expected_seq)) = test_contigs
                .iter()
                .find(|(s, c, _)| s == &sample_name && c == &contig_name)
            {
                print!("  {contig_name} ... ");

                if validate_bases(
                    &decompressed,
                    &format!("OUTPUT {sample_name}/{contig_name}"),
                ) {
                    if decompressed == *expected_seq {
                        println!("✓ Matches input");
                    } else {
                        println!("❌ Mismatch with input!");
                        println!(
                            "   Expected len: {}, got len: {}",
                            expected_seq.len(),
                            decompressed.len()
                        );

                        // Find first difference
                        for (i, (exp, got)) in
                            expected_seq.iter().zip(decompressed.iter()).enumerate()
                        {
                            if exp != got {
                                println!(
                                    "   First diff at position {i}: expected {exp}, got {got}"
                                );
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(())
}
