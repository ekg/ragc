// Basic AGC Archive Extraction Example
//
// This example demonstrates:
// - Opening an AGC archive
// - Listing all samples
// - Extracting a single sample
// - Accessing contig names and sequences

use anyhow::Result;
use ragc_core::{Decompressor, DecompressorConfig};

fn main() -> Result<()> {
    // Open the AGC archive
    let archive_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "data.agc".to_string());

    println!("Opening archive: {}", archive_path);

    let mut decompressor = Decompressor::open(&archive_path, DecompressorConfig::default())?;

    // List all samples in the archive
    let samples = decompressor.list_samples();
    println!("\nFound {} samples:", samples.len());
    for (i, sample) in samples.iter().enumerate() {
        println!("  {}: {}", i + 1, sample);
        if i >= 9 {
            println!("  ... ({} more)", samples.len() - 10);
            break;
        }
    }

    // Extract the first sample (if available)
    if let Some(first_sample) = samples.first() {
        println!("\nExtracting sample: {}", first_sample);

        let contigs = decompressor.get_sample(first_sample)?;

        println!("  Found {} contigs:", contigs.len());
        for (contig_name, sequence) in &contigs {
            println!(
                "    {}: {} bp",
                contig_name,
                sequence.len()
            );
        }

        // Display first 100bp of first contig
        if let Some((contig_name, sequence)) = contigs.first() {
            let preview = String::from_utf8_lossy(
                &sequence[..sequence.len().min(100)]
            );
            println!("\nFirst 100bp of {}:", contig_name);
            println!("{}", preview);
        }
    }

    // Close the decompressor
    decompressor.close()?;

    Ok(())
}
