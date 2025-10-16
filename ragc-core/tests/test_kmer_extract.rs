// Integration test for k-mer extraction
// Should produce output identical to C++ test_kmer_extract

use ragc_core::{enumerate_kmers, remove_non_singletons};

fn print_kmers(vec: &[u64], max_count: usize) {
    println!("Count: {}", vec.len());
    print!("K-mers: ");
    for (i, &kmer) in vec.iter().enumerate().take(max_count) {
        if i > 0 {
            print!(" ");
        }
        print!("{:x}", kmer);
    }
    println!();
}

fn main() {
    // Test 1: Simple sequence
    println!("# Test 1: Simple sequence (ACGTACGT)");
    {
        let ctg = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let kmers = enumerate_kmers(&ctg, 3);
        print_kmers(&kmers, 50);
    }

    // Test 2: Sequence with N (value > 3)
    println!("\n# Test 2: Sequence with N (AC[N]GTAC)");
    {
        let ctg = vec![0, 1, 4, 2, 3, 0, 1]; // AC[N]GTAC
        let kmers = enumerate_kmers(&ctg, 3);
        print_kmers(&kmers, 50);
    }

    // Test 3: Short sequence (too short for k-mer)
    println!("\n# Test 3: Short sequence (AC, k=3)");
    {
        let ctg = vec![0, 1]; // AC
        let kmers = enumerate_kmers(&ctg, 3);
        print_kmers(&kmers, 50);
    }

    // Test 4: Longer sequence
    println!("\n# Test 4: Longer sequence (ACGTACGTACGT, k=5)");
    {
        let ctg = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGTACGT
        let kmers = enumerate_kmers(&ctg, 5);
        print_kmers(&kmers, 50);
    }

    // Test 5: Remove non-singletons
    println!("\n# Test 5: Remove non-singletons");
    {
        let mut vec = vec![1, 2, 2, 3, 3, 3, 4, 5, 5, 6];
        print!("Before: ");
        print_kmers(&vec, 50);

        remove_non_singletons(&mut vec, 0);

        print!("After: ");
        print_kmers(&vec, 50);
    }

    // Test 6: Remove non-singletons with virtual_begin
    println!("\n# Test 6: Remove non-singletons with virtual_begin=2");
    {
        let mut vec = vec![1, 1, 2, 3, 3, 4, 5, 5];
        print!("Before: ");
        print_kmers(&vec, 50);

        remove_non_singletons(&mut vec, 2);

        print!("After: ");
        print_kmers(&vec, 50);
    }

    // Test 7: Realistic example with duplicates
    println!("\n# Test 7: Find candidate k-mers (with sorting)");
    {
        let ctg = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT (repeating)
        let mut kmers = enumerate_kmers(&ctg, 3);

        print!("All k-mers: ");
        print_kmers(&kmers, 50);

        kmers.sort_unstable();
        print!("Sorted: ");
        print_kmers(&kmers, 50);

        remove_non_singletons(&mut kmers, 0);
        print!("Singletons: ");
        print_kmers(&kmers, 50);
    }

    // Test 8: Unique sequence (all singletons)
    println!("\n# Test 8: All unique k-mers (ACGTCATG)");
    {
        let ctg = vec![0, 1, 2, 3, 1, 0, 3, 2]; // ACGTCATG
        let mut kmers = enumerate_kmers(&ctg, 3);

        print!("All k-mers: ");
        print_kmers(&kmers, 50);

        kmers.sort_unstable();
        print!("Sorted: ");
        print_kmers(&kmers, 50);

        remove_non_singletons(&mut kmers, 0);
        print!("Singletons: ");
        print_kmers(&kmers, 50);
    }
}
