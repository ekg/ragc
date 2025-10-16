// Integration test for hash functions
// Should produce output identical to C++ test_hash

use ragc_common::{MurMur64Hash, MurMurPair64Hash, MurMurStringsHash, MurMur32Hash};

fn main() {
    // Test 1: MurMur64Hash for various values
    println!("# Test 1: MurMur64Hash");
    let test_values: Vec<u64> = vec![0, 1, 2, 42, 100, 1000, 0xDEADBEEF, 0xFFFFFFFFFFFFFFFF];

    for val in test_values {
        println!("{:x}\t{:x}", val, MurMur64Hash::hash(val));
    }

    // Test 2: MurMurPair64Hash for pairs
    println!("\n# Test 2: MurMurPair64Hash");
    let test_pairs = vec![
        (0u64, 0u64),
        (1, 2),
        (42, 43),
        (100, 200),
        (0xDEADBEEF, 0xCAFEBABE),
        (0xFFFFFFFFFFFFFFFF, 0),
    ];

    for (first, second) in test_pairs {
        println!("{:x}\t{:x}\t{:x}", first, second, MurMurPair64Hash::hash(first, second));
    }

    // Test 3: MurMurStringsHash
    println!("\n# Test 3: MurMurStringsHash");
    let test_strings = vec![
        "",
        "a",
        "hello",
        "world",
        "The quick brown fox jumps over the lazy dog",
        "ACGT",
        "sample_name_123",
    ];

    for s in test_strings {
        println!("\"{}\"\t{:x}", s, MurMurStringsHash::hash(s));
    }

    // Test 4: MurMur32Hash
    println!("\n# Test 4: MurMur32Hash");
    let test_values32: Vec<u32> = vec![0, 1, 2, 42, 100, 1000, 0xDEADBEEF, 0xFFFFFFFF];

    for val in test_values32 {
        println!("{:x}\t{:x}", val, MurMur32Hash::hash(val));
    }
}
