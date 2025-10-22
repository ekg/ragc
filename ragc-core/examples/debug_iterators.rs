// Debug program to compare what the two iterators produce

use anyhow::Result;
use ragc_core::contig_iterator::{ContigIterator, MultiFileIterator, PansnFileIterator};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
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

    println!("=== Testing PansnFileIterator ===");
    let mut pansn_data: HashMap<(String, String), Vec<u8>> = HashMap::new();
    let mut pansn_iter = PansnFileIterator::new(&pansn_file)?;
    let mut pansn_count = 0;

    while let Some((sample_name, contig_name, sequence)) = pansn_iter.next_contig()? {
        let key = (sample_name, contig_name);
        pansn_data.insert(key, sequence);
        pansn_count += 1;
    }
    println!("Pansn iterator produced {pansn_count} contigs");

    println!("\n=== Testing MultiFileIterator ===");
    let mut multi_data: HashMap<(String, String), Vec<u8>> = HashMap::new();
    let mut multi_iter = MultiFileIterator::new(fasta_files)?;
    let mut multi_count = 0;

    while let Some((sample_name, contig_name, sequence)) = multi_iter.next_contig()? {
        let key = (sample_name.clone(), contig_name.clone());
        multi_data.insert(key, sequence);
        multi_count += 1;
    }
    println!("Multi iterator produced {multi_count} contigs");

    println!("\n=== Comparing Results ===");

    // Check for missing contigs in multi
    for key in pansn_data.keys() {
        if !multi_data.contains_key(key) {
            println!("❌ Missing in multi: {} / {}", key.0, key.1);
        }
    }

    // Check for extra contigs in multi
    for key in multi_data.keys() {
        if !pansn_data.contains_key(key) {
            println!("❌ Extra in multi: {} / {}", key.0, key.1);
        }
    }

    // Compare sequences
    let mut differences = 0;
    for (key, pansn_seq) in &pansn_data {
        if let Some(multi_seq) = multi_data.get(key) {
            if pansn_seq != multi_seq {
                differences += 1;

                let mut pansn_hasher = Sha256::new();
                pansn_hasher.update(pansn_seq);
                let pansn_hash = format!("{:x}", pansn_hasher.finalize());

                let mut multi_hasher = Sha256::new();
                multi_hasher.update(multi_seq);
                let multi_hash = format!("{:x}", multi_hasher.finalize());

                println!(
                    "❌ DIFFERENT: {} / {} (pansn len: {}, multi len: {})",
                    key.0,
                    key.1,
                    pansn_seq.len(),
                    multi_seq.len()
                );
                println!("   Pansn SHA256: {pansn_hash}");
                println!("   Multi SHA256: {multi_hash}");

                if differences <= 3 {
                    // Show first few bytes
                    let show_bytes = 50.min(pansn_seq.len()).min(multi_seq.len());
                    println!(
                        "   Pansn first {} bytes: {:?}",
                        show_bytes,
                        &pansn_seq[..show_bytes]
                    );
                    println!(
                        "   Multi first {} bytes: {:?}",
                        show_bytes,
                        &multi_seq[..show_bytes]
                    );
                }
            }
        }
    }

    if differences == 0 {
        println!("✓ All sequences match!");
    } else {
        println!("\n❌ Found {differences} sequences with differences");
    }

    Ok(())
}
