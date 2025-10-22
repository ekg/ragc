// Check if contig order differs between iterators

use anyhow::Result;
use ragc_core::contig_iterator::{ContigIterator, MultiFileIterator, PansnFileIterator};
use std::path::{Path, PathBuf};

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

    println!("=== Collecting contig order from PansnFileIterator ===");
    let mut pansn_order = Vec::new();
    let mut pansn_iter = PansnFileIterator::new(&pansn_file)?;

    while let Some((sample_name, contig_name, _sequence)) = pansn_iter.next_contig()? {
        pansn_order.push((sample_name, contig_name));
    }

    println!("Pansn iterator produced {} contigs", pansn_order.len());
    println!("First 10 contigs:");
    for (i, (sample, contig)) in pansn_order.iter().take(10).enumerate() {
        println!("  {i}: {sample} / {contig}");
    }

    println!("\n=== Collecting contig order from MultiFileIterator ===");
    let mut multi_order = Vec::new();
    let mut multi_iter = MultiFileIterator::new(fasta_files)?;

    while let Some((sample_name, contig_name, _sequence)) = multi_iter.next_contig()? {
        multi_order.push((sample_name, contig_name));
    }

    println!("Multi iterator produced {} contigs", multi_order.len());
    println!("First 10 contigs:");
    for (i, (sample, contig)) in multi_order.iter().take(10).enumerate() {
        println!("  {i}: {sample} / {contig}");
    }

    println!("\n=== Comparing Orders ===");
    let mut first_diff = None;
    for (i, (pansn_entry, multi_entry)) in pansn_order.iter().zip(multi_order.iter()).enumerate() {
        if pansn_entry != multi_entry {
            first_diff = Some(i);
            break;
        }
    }

    if let Some(i) = first_diff {
        println!("❌ First difference at position {i}");
        println!("   Pansn: {} / {}", pansn_order[i].0, pansn_order[i].1);
        println!("   Multi: {} / {}", multi_order[i].0, multi_order[i].1);

        println!(
            "\nContext (positions {} to {}):",
            i.saturating_sub(2),
            (i + 3).min(pansn_order.len())
        );
        for j in i.saturating_sub(2)..(i + 3).min(pansn_order.len()) {
            println!(
                "  Pansn[{}]: {} / {}",
                j, pansn_order[j].0, pansn_order[j].1
            );
            println!(
                "  Multi[{}]: {} / {}",
                j, multi_order[j].0, multi_order[j].1
            );
            println!();
        }
    } else {
        println!("✓ Orders are identical!");
    }

    Ok(())
}
