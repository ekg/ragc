// Integration test for splitter identification
// Should produce output identical to C++ test_splitters

use ragc_core::determine_splitters;

fn print_splitters(splitters: &std::collections::HashSet<u64>, max_count: usize) {
    println!("Count: {}", splitters.len());
    print!("Splitters: ");

    // Convert to sorted vector for consistent output
    let mut sorted: Vec<u64> = splitters.iter().copied().collect();
    sorted.sort_unstable();

    for (i, &s) in sorted.iter().enumerate().take(max_count) {
        if i > 0 {
            print!(" ");
        }
        print!("{s:x}");
    }
    println!();
}

fn main() {
    // Test 1: Two contigs with truly different k-mers
    println!("# Test 1: Simple case (AAAC, GGGT)");
    {
        let contigs = vec![
            vec![0, 0, 0, 1], // AAAC
            vec![2, 2, 2, 3], // GGGT
        ];

        let splitters = determine_splitters(&contigs, 3);
        print_splitters(&splitters, 50);
    }

    // Test 2: Identical contigs (no singletons)
    println!("\n# Test 2: No singletons (identical contigs)");
    {
        let contigs = vec![vec![0, 1, 2, 3], vec![0, 1, 2, 3]];

        let splitters = determine_splitters(&contigs, 3);
        print_splitters(&splitters, 50);
    }

    // Test 3: All unique k-mers
    println!("\n# Test 3: All unique (AAAA, CCCC)");
    {
        let contigs = vec![
            vec![0, 0, 0, 0], // AAAA
            vec![1, 1, 1, 1], // CCCC
        ];

        let splitters = determine_splitters(&contigs, 3);
        print_splitters(&splitters, 50);
    }

    // Test 4: Mixed scenario
    println!("\n# Test 4: Mixed scenario");
    {
        let contigs = vec![
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3], // ACGTACGTACGT (repeating)
            vec![0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2], // AAAACCCCGGGG (unique patterns)
        ];

        let splitters = determine_splitters(&contigs, 3);
        print_splitters(&splitters, 50);
    }

    // Test 5: Longer realistic example
    println!("\n# Test 5: Longer sequences");
    {
        let contigs = vec![
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3], // ACGTACGTACGTACGTACGT
            vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0], // TGCATGCATGCATGCATGCA
            vec![0, 0, 1, 1, 2, 2, 3, 3, 0, 0, 1, 1, 2, 2, 3, 3, 0, 0, 1, 1], // AACCGGTTAACCGGTTAACC
        ];

        let splitters = determine_splitters(&contigs, 5);
        print_splitters(&splitters, 100);
    }

    // Test 6: Single contig (all singletons)
    println!("\n# Test 6: Single contig");
    {
        let contigs = vec![
            vec![0, 1, 2, 3, 1, 0, 3, 2, 0, 2, 1, 3], // ACGTCATGACGT (mostly unique)
        ];

        let splitters = determine_splitters(&contigs, 3);
        print_splitters(&splitters, 100);
    }
}
