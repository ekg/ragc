#![allow(unexpected_cfgs)]
#![allow(unused_imports)]
#![allow(clippy::len_zero)]

// Compare RAGC splitter detection vs C++ AGC splitter detection
use ragc_core::{determine_splitters, GenomeIO};
use std::collections::HashSet;
use std::fs::File;

fn main() -> anyhow::Result<()> {
    let test_file = "/tmp/yeast_test.fa";
    let k = 21;
    let segment_size = 10000;

    println!("Loading contigs from: {}", test_file);

    // Load contigs using RAGC's GenomeIO
    let file = File::open(test_file)?;
    let mut gio = GenomeIO::new(file);

    let mut contigs = Vec::new();
    while let Some((name, contig)) = gio.read_contig_converted()? {
        println!(
            "  Read contig '{}': {} bases (converted)",
            name,
            contig.len()
        );
        // Show first 20 bases in numeric encoding
        if contig.len() > 0 {
            print!("    First 20 bases (0=A,1=C,2=G,3=T): ");
            for &b in contig.iter().take(20) {
                print!("{} ", b);
            }
            println!();
        }
        contigs.push(contig);
    }

    println!("\nLoaded {} contigs total", contigs.len());
    for (i, contig) in contigs.iter().enumerate() {
        println!("  Contig {}: {} bases", i, contig.len());
    }

    // Run RAGC's splitter detection
    println!("\n=== RAGC Splitter Detection ===");
    let (splitters, singletons, duplicates) = determine_splitters(&contigs, k, segment_size);

    println!("RAGC Results:");
    println!("  Singletons: {}", singletons.len());
    println!("  Duplicates: {}", duplicates.len());
    println!("  Splitters:  {}", splitters.len());

    // Show first 10 splitters (sorted)
    let mut splitter_vec: Vec<_> = splitters.iter().cloned().collect();
    splitter_vec.sort();
    println!("\n  First 10 splitters:");
    for (i, sp) in splitter_vec.iter().take(10).enumerate() {
        println!("    {}: {:016x}", i, sp);
    }

    println!("\n=== C++ AGC Splitter Detection ===");
    println!("Create archive with C++ AGC to see its splitter count in output...");
    println!(
        "Run: /home/erik/agc/bin/agc create -o /tmp/cpp_test.agc -k 21 -s 10000 -l 20 -t 1 {}",
        test_file
    );
    println!("Look for 'No. of splitters:' in output");

    Ok(())
}
