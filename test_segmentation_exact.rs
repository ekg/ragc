#!/usr/bin/env rust-script
//! Test to compare C++ AGC's exact segmentation algorithm with RAGC's
//!
//! Purpose: Prove whether RAGC creates the same segment boundaries as C++ AGC

use std::collections::HashSet;

type Contig = Vec<u8>;

const MISSING_KMER: u64 = u64::MAX;

#[derive(Debug, Clone)]
struct Segment {
    data: Contig,
    front_kmer: u64,
    back_kmer: u64,
}

impl Segment {
    fn new(data: Contig, front_kmer: u64, back_kmer: u64) -> Self {
        Segment { data, front_kmer, back_kmer }
    }

    fn len(&self) -> usize {
        self.data.len()
    }
}

/// Simplified k-mer hasher for canonical k-mers
struct SimpleKmer {
    k: usize,
    mask: u64,
    data: u64,
    count: usize,
}

impl SimpleKmer {
    fn new(k: usize) -> Self {
        let mask = (1u64 << (2 * k)) - 1;
        SimpleKmer {
            k,
            mask,
            data: 0,
            count: 0,
        }
    }

    fn reset(&mut self) {
        self.data = 0;
        self.count = 0;
    }

    fn insert(&mut self, base: u8) {
        self.data = ((self.data << 2) | (base as u64)) & self.mask;
        if self.count < self.k {
            self.count += 1;
        }
    }

    fn is_full(&self) -> bool {
        self.count == self.k
    }

    fn value(&self) -> u64 {
        self.data
    }
}

/// C++ AGC's EXACT algorithm (direct port of compress_contig)
fn split_cpp_agc(contig: &Contig, splitters: &HashSet<u64>, k: usize) -> Vec<Segment> {
    let mut segments = Vec::new();

    if contig.len() < k {
        return vec![Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER)];
    }

    let mut kmer = SimpleKmer::new(k);
    let mut pos: usize = 0;
    let mut split_pos: usize = 0;
    let mut front_kmer: u64 = MISSING_KMER;

    // Main loop (C++ AGC lines 2010-2039)
    for &base in contig.iter() {
        if base > 3 {
            kmer.reset();
        } else {
            kmer.insert(base);

            if kmer.is_full() {
                let kmer_value = kmer.value();

                // Check if splitter
                if splitters.contains(&kmer_value) {
                    // Create segment
                    if pos + 1 > split_pos {
                        let segment_data = contig[split_pos..pos + 1].to_vec();
                        if !segment_data.is_empty() {
                            segments.push(Segment::new(segment_data, front_kmer, kmer_value));
                        }
                    }

                    // Update for next segment (with k-base overlap)
                    split_pos = (pos + 1).saturating_sub(k);
                    front_kmer = kmer_value;
                    kmer.reset();
                }
            }
        }

        pos += 1;
    }

    // Add final segment (C++ AGC lines 2049-2051)
    if split_pos < contig.len() {
        let segment_data = contig[split_pos..].to_vec();
        if !segment_data.is_empty() {
            segments.push(Segment::new(segment_data, front_kmer, MISSING_KMER));
        }
    }

    if segments.is_empty() {
        segments.push(Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER));
    }

    segments
}

/// RAGC's current algorithm (simplified version)
fn split_ragc(contig: &Contig, splitters: &HashSet<u64>, k: usize) -> Vec<Segment> {
    let mut segments = Vec::new();

    if contig.len() < k {
        return vec![Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER)];
    }

    let mut kmer = SimpleKmer::new(k);
    let mut segment_start = 0;
    let mut front_kmer = MISSING_KMER;
    let mut recent_kmers: Vec<(usize, u64)> = Vec::new();

    // Main loop
    for (pos, &base) in contig.iter().enumerate() {
        if base > 3 {
            kmer.reset();
            recent_kmers.clear();
        } else {
            kmer.insert(base);

            if kmer.is_full() {
                let kmer_value = kmer.value();
                recent_kmers.push((pos, kmer_value));

                if splitters.contains(&kmer_value) {
                    let segment_end = pos + 1;
                    let segment_data = contig[segment_start..segment_end].to_vec();

                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, kmer_value));
                    }

                    segment_start = (pos + 1).saturating_sub(k);
                    front_kmer = kmer_value;
                    recent_kmers.clear();
                    kmer.reset();
                }
            }
        }
    }

    // RAGC's extra backward scan (NOT in C++ AGC)
    for (pos, kmer_value) in recent_kmers.iter().rev() {
        if splitters.contains(kmer_value) {
            let segment_end = pos + 1;
            let remaining_after = contig.len() - segment_end;

            if remaining_after > k {
                if segment_end > segment_start {
                    let segment_data = contig[segment_start..segment_end].to_vec();
                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, *kmer_value));
                        segment_start = (pos + 1).saturating_sub(k);
                        front_kmer = *kmer_value;
                    }
                }
                break;
            }
        }
    }

    // Add final segment
    if segment_start < contig.len() {
        let segment_data = contig[segment_start..].to_vec();
        if !segment_data.is_empty() {
            segments.push(Segment::new(segment_data, front_kmer, MISSING_KMER));
        }
    }

    // RAGC's merging logic (NOT in C++ AGC)
    if segments.len() >= 2 {
        let last_idx = segments.len() - 1;
        if segments[last_idx].data.len() < k {
            let last_seg = segments.pop().unwrap();
            let second_last = segments.last_mut().unwrap();
            second_last.data.extend_from_slice(&last_seg.data);
            second_last.back_kmer = last_seg.back_kmer;
        }
    }

    if segments.is_empty() {
        segments.push(Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER));
    }

    segments
}

fn main() {
    println!("=== Testing C++ AGC vs RAGC Segmentation ===\n");

    // Test 1: Simple contig with one splitter
    let k = 3;
    let contig1 = vec![0u8, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]; // AAACCCGGGTTT

    // Create splitter for CCC (1,1,1)
    let mut kmer = SimpleKmer::new(k);
    kmer.insert(1);
    kmer.insert(1);
    kmer.insert(1);
    let ccc_kmer = kmer.value();

    let mut splitters1 = HashSet::new();
    splitters1.insert(ccc_kmer);

    let cpp_seg1 = split_cpp_agc(&contig1, &splitters1, k);
    let ragc_seg1 = split_ragc(&contig1, &splitters1, k);

    println!("Test 1: Contig with one splitter");
    println!("  Contig: AAACCCGGGTTT (length {})", contig1.len());
    println!("  C++ AGC segments: {}", cpp_seg1.len());
    for (i, seg) in cpp_seg1.iter().enumerate() {
        println!("    {}: len={}, front={}, back={}", i, seg.len(), seg.front_kmer, seg.back_kmer);
    }
    println!("  RAGC segments: {}", ragc_seg1.len());
    for (i, seg) in ragc_seg1.iter().enumerate() {
        println!("    {}: len={}, front={}, back={}", i, seg.len(), seg.front_kmer, seg.back_kmer);
    }

    if cpp_seg1.len() != ragc_seg1.len() {
        println!("  ❌ MISMATCH: Different segment counts!");
    } else {
        println!("  ✓ Segment counts match");
    }
    println!();

    // Test 2: Longer contig with multiple splitters
    let contig2 = vec![0u8; 100]; // 100 A's
    let mut splitters2 = HashSet::new();

    // Create splitter AAA (0,0,0)
    let mut kmer2 = SimpleKmer::new(k);
    kmer2.insert(0);
    kmer2.insert(0);
    kmer2.insert(0);
    splitters2.insert(kmer2.value());

    let cpp_seg2 = split_cpp_agc(&contig2, &splitters2, k);
    let ragc_seg2 = split_ragc(&contig2, &splitters2, k);

    println!("Test 2: Homopolymer contig (100 A's)");
    println!("  C++ AGC segments: {}", cpp_seg2.len());
    println!("  RAGC segments: {}", ragc_seg2.len());

    if cpp_seg2.len() != ragc_seg2.len() {
        println!("  ❌ MISMATCH: Different segment counts!");
        println!("  Difference: {} segments", (cpp_seg2.len() as i32 - ragc_seg2.len() as i32).abs());
    } else {
        println!("  ✓ Segment counts match");
    }
    println!();

    // Summary
    if cpp_seg1.len() == ragc_seg1.len() && cpp_seg2.len() == ragc_seg2.len() {
        println!("✅ CONCLUSION: Algorithms produce same segment counts");
    } else {
        println!("❌ CONCLUSION: Algorithms produce DIFFERENT segment counts");
        println!("   This proves RAGC has a bug in segmentation logic!");
    }
}
