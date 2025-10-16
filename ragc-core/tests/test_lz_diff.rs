// Integration test for LZ diff encoding/decoding
// Should produce output identical to C++ test_lz_diff

use ragc_core::LZDiff;
use ragc_common::types::Contig;

fn print_contig(label: &str, contig: &Contig) {
    print!("{} len:{} data:", label, contig.len());
    for (i, &byte) in contig.iter().enumerate() {
        if i >= 50 {
            break;
        }
        if i > 0 {
            print!(" ");
        }
        print!("{}", byte);
    }
    if contig.len() > 50 {
        print!(" ...");
    }
    println!();
}

fn print_encoded(label: &str, encoded: &[u8]) {
    print!("{} len:{} hex:", label, encoded.len());
    for (i, &byte) in encoded.iter().enumerate() {
        if i >= 100 {
            break;
        }
        print!("{:02x}", byte);
    }
    if encoded.len() > 100 {
        print!("...");
    }
    println!();
}

fn test_encode_decode(reference: &Contig, target: &Contig, lz_diff: &mut LZDiff) {
    let encoded = lz_diff.encode(target);

    print_contig("reference", reference);
    print_contig("target", target);
    print_encoded("encoded", &encoded);

    let mut decoded = lz_diff.decode(&encoded);

    // If encoded is empty and they're equal length, target equals reference (V2 optimization)
    if encoded.is_empty() && target.len() == reference.len() {
        decoded = reference.clone();
    }

    print_contig("decoded", &decoded);
    println!("match:{}", if *target == decoded { "true" } else { "false" });
}

fn test_simple_match() {
    println!("# Test 1: Simple match");

    // Reference: AAACCCGGGTTT (0,0,0,1,1,1,2,2,2,3,3,3)
    let reference = vec![0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3];

    // Target: identical to reference
    let target = vec![0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3];

    let mut lz_diff = LZDiff::new(18);
    lz_diff.prepare(&reference);

    test_encode_decode(&reference, &target, &mut lz_diff);
    println!();
}

fn test_with_literals() {
    println!("# Test 2: With literals");

    // Reference: 20 bases
    let reference = vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3];

    // Target: Some matches, some literals
    let target = vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0];

    let mut lz_diff = LZDiff::new(18);
    lz_diff.prepare(&reference);

    test_encode_decode(&reference, &target, &mut lz_diff);
    println!();
}

fn test_with_nruns() {
    println!("# Test 3: With N-runs");

    // Reference: 20 bases
    let reference = vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3];

    // Target: With N-run (4=N code)
    let target = vec![0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3];

    let mut lz_diff = LZDiff::new(18);
    lz_diff.prepare(&reference);

    test_encode_decode(&reference, &target, &mut lz_diff);
    println!();
}

fn test_short_sequence() {
    println!("# Test 4: Short sequence (< min_match_len)");

    // Reference: 10 bases
    let reference = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];

    // Target: 10 bases (too short for matching with min_match_len=18)
    let target = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];

    let mut lz_diff = LZDiff::new(18);
    lz_diff.prepare(&reference);

    test_encode_decode(&reference, &target, &mut lz_diff);
    println!();
}

fn test_long_match() {
    println!("# Test 5: Long match");

    // Reference: 50 bases of pattern
    let mut reference = Vec::new();
    for i in 0..50 {
        reference.push((i % 4) as u8);
    }

    // Target: identical
    let target = reference.clone();

    let mut lz_diff = LZDiff::new(18);
    lz_diff.prepare(&reference);

    test_encode_decode(&reference, &target, &mut lz_diff);
    println!();
}

fn test_partial_match() {
    println!("# Test 6: Partial match");

    // Reference: 30 bases
    let mut reference = Vec::new();
    for i in 0..30 {
        reference.push((i % 4) as u8);
    }

    // Target: First 20 bases match, then 10 different
    let mut target = Vec::new();
    for i in 0..20 {
        target.push((i % 4) as u8);
    }
    for i in 0..10 {
        target.push(((i + 1) % 4) as u8); // Different pattern
    }

    let mut lz_diff = LZDiff::new(18);
    lz_diff.prepare(&reference);

    test_encode_decode(&reference, &target, &mut lz_diff);
    println!();
}

fn main() {
    test_simple_match();
    test_with_literals();
    test_with_nruns();
    test_short_sequence();
    test_long_match();
    test_partial_match();
}
