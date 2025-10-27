// Prefix-based Sample Extraction Example
//
// This example demonstrates:
// - Listing samples by prefix
// - Extracting multiple samples matching a prefix
// - Useful for extracting specific genomes or haplotypes

use anyhow::Result;
use ragc_core::{Decompressor, DecompressorConfig};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <archive.agc> [prefix]", args[0]);
        eprintln!();
        eprintln!("Examples:");
        eprintln!(
            "  {} data.agc AAA         # Extract all samples starting with 'AAA'",
            args[0]
        );
        eprintln!(
            "  {} data.agc AAA#0       # Extract haplotype 0 of AAA",
            args[0]
        );
        eprintln!("  {} data.agc             # List all samples", args[0]);
        std::process::exit(1);
    }

    let archive_path = &args[1];
    let prefix = args.get(2).map(String::as_str);

    println!("Opening archive: {archive_path}");

    let mut decompressor = Decompressor::open(archive_path, DecompressorConfig::default())?;

    if let Some(prefix) = prefix {
        // Extract samples matching the prefix
        println!("\nSearching for samples with prefix: '{prefix}'");

        let matching_samples = decompressor.list_samples_with_prefix(prefix);

        if matching_samples.is_empty() {
            println!("No samples found matching prefix '{prefix}'");
            return Ok(());
        }

        println!("Found {} matching samples:", matching_samples.len());
        for sample in &matching_samples {
            println!("  - {sample}");
        }

        // Extract all matching samples
        println!("\nExtracting samples...");
        let samples_data = decompressor.get_samples_by_prefix(prefix)?;

        // Display statistics
        let mut total_contigs = 0;
        let mut total_bases = 0;

        for (sample_name, contigs) in &samples_data {
            let sample_bases: usize = contigs.iter().map(|(_, seq)| seq.len()).sum();
            println!(
                "  {}: {} contigs, {} bp",
                sample_name,
                contigs.len(),
                sample_bases
            );
            total_contigs += contigs.len();
            total_bases += sample_bases;
        }

        println!(
            "\nTotal: {} samples, {} contigs, {} bp",
            samples_data.len(),
            total_contigs,
            total_bases
        );

        // Optionally write to FASTA
        if args.len() >= 4 && args[3] == "--output" {
            if let Some(output_path) = args.get(4) {
                write_fasta(&samples_data, output_path)?;
                println!("\nWrote output to: {output_path}");
            }
        }
    } else {
        // Just list all samples
        let samples = decompressor.list_samples();
        println!("\nFound {} samples:", samples.len());
        for sample in &samples {
            println!("  {sample}");
        }
    }

    decompressor.close()?;

    Ok(())
}

fn write_fasta(
    samples_data: &std::collections::HashMap<String, Vec<(String, Vec<u8>)>>,
    output_path: &str,
) -> Result<()> {
    use std::io::Write;

    let mut file = std::fs::File::create(output_path)?;

    for contigs in samples_data.values() {
        for (contig_name, sequence) in contigs {
            writeln!(file, ">{contig_name}")?;
            file.write_all(sequence)?;
            writeln!(file)?;
        }
    }

    Ok(())
}
