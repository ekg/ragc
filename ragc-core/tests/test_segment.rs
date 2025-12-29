#![allow(clippy::all)]
// Integration test for segmentation
// Should produce output identical to C++ test_segment

use ahash::AHashSet;
use ragc_core::{split_at_splitters, Kmer, KmerMode};

fn print_segment(seg: &ragc_core::Segment, max_data_len: usize) {
    print!("front_kmer:{:x}", seg.front_kmer);
    print!(" back_kmer:{:x}", seg.back_kmer);
    print!(" len:{}", seg.len());
    print!(" data:");
    for (i, &byte) in seg.data.iter().enumerate() {
        if i >= max_data_len {
            break;
        }
        print!(" {byte}");
    }
    if seg.data.len() > max_data_len {
        print!(" ...");
    }
    println!();
}

fn main() {
    // Test 1: No splitters
    println!("# Test 1: No splitters (entire contig)");
    {
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let splitters = AHashSet::new();
        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }

    // Test 2: With one splitter
    println!("\n# Test 2: With one splitter");
    {
        let contig = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];

        // Get the first k-mer as a splitter
        let mut kmer = Kmer::new(3, KmerMode::Canonical);
        kmer.insert(0);
        kmer.insert(0);
        kmer.insert(0);
        let first_kmer = kmer.data();

        let mut splitters = AHashSet::new();
        splitters.insert(first_kmer);

        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }

    // Test 3: Multiple splitters
    println!("\n# Test 3: Multiple splitters");
    {
        let contig = vec![0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2];

        // Get AAA and CCC as splitters
        let mut kmer_aaa = Kmer::new(3, KmerMode::Canonical);
        kmer_aaa.insert(0);
        kmer_aaa.insert(0);
        kmer_aaa.insert(0);

        let mut kmer_ccc = Kmer::new(3, KmerMode::Canonical);
        kmer_ccc.insert(1);
        kmer_ccc.insert(1);
        kmer_ccc.insert(1);

        let mut splitters = AHashSet::new();
        splitters.insert(kmer_aaa.data());
        splitters.insert(kmer_ccc.data());

        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }

    // Test 4: Short contig (too short for k-mers)
    println!("\n# Test 4: Short contig");
    {
        let contig = vec![0, 1];
        let mut splitters = AHashSet::new();
        splitters.insert(12345);

        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }

    // Test 5: Splitter at start
    println!("\n# Test 5: Splitter at start");
    {
        let contig = vec![0, 0, 0, 1, 2, 3];

        let mut kmer = Kmer::new(3, KmerMode::Canonical);
        kmer.insert(0);
        kmer.insert(0);
        kmer.insert(0);

        let mut splitters = AHashSet::new();
        splitters.insert(kmer.data());

        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }

    // Test 6: Consecutive splitters
    println!("\n# Test 6: Consecutive splitters");
    {
        let contig = vec![0, 0, 0, 0, 0, 0];

        let mut kmer = Kmer::new(3, KmerMode::Canonical);
        kmer.insert(0);
        kmer.insert(0);
        kmer.insert(0);

        let mut splitters = AHashSet::new();
        splitters.insert(kmer.data());

        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }

    // Test 7: With N bases (non-ACGT)
    println!("\n# Test 7: With N bases");
    {
        let contig = vec![0, 0, 0, 4, 1, 1, 1];

        let mut splitters = AHashSet::new();
        splitters.insert(12345); // Some random splitter

        let segments = split_at_splitters(&contig, &splitters, 3);

        println!("num_segments:{}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            print!("seg{i} ");
            print_segment(seg, 20);
        }
    }
}
