// Integration test for k-mer operations
// Should produce output identical to C++ test_kmer

use ragc_core::{reverse_complement, Kmer, KmerMode};

fn main() {
    // Test 1: reverse_complement for individual bases
    println!("# Test 1: reverse_complement");
    for i in 0..4 {
        println!("{}\t{}", i, reverse_complement(i));
    }

    // Test 2: Create k-mer and insert bases (canonical mode)
    println!("\n# Test 2: Insert ACGT (canonical)");
    let mut kmer1 = Kmer::new(4, KmerMode::Canonical);

    // Insert A, C, G, T (0, 1, 2, 3)
    kmer1.insert(0); // A
    kmer1.insert(1); // C
    kmer1.insert(2); // G
    kmer1.insert(3); // T

    println!("is_full\t{}", if kmer1.is_full() { 1 } else { 0 });
    println!("cur_size\t{}", kmer1.get_cur_size());
    println!("max_size\t{}", kmer1.get_max_size());
    println!("data_dir\t{:x}", kmer1.data_dir());
    println!("data_rc\t{:x}", kmer1.data_rc());
    println!("data_canonical\t{:x}", kmer1.data_canonical());
    println!(
        "is_dir_oriented\t{}",
        if kmer1.is_dir_oriented() { 1 } else { 0 }
    );

    // Test 3: Create k-mer and insert bases (direct mode)
    println!("\n# Test 3: Insert ACGT (direct)");
    let mut kmer2 = Kmer::new(4, KmerMode::Direct);
    kmer2.insert(0); // A
    kmer2.insert(1); // C
    kmer2.insert(2); // G
    kmer2.insert(3); // T

    println!("is_full\t{}", if kmer2.is_full() { 1 } else { 0 });
    println!("data_dir\t{:x}", kmer2.data_dir());

    // Test 4: Create k-mer and insert bases (rev_comp mode)
    println!("\n# Test 4: Insert ACGT (rev_comp)");
    let mut kmer3 = Kmer::new(4, KmerMode::RevComp);
    kmer3.insert(0); // A
    kmer3.insert(1); // C
    kmer3.insert(2); // G
    kmer3.insert(3); // T

    println!("is_full\t{}", if kmer3.is_full() { 1 } else { 0 });
    println!("data_rc\t{:x}", kmer3.data_rc());

    // Test 5: Insert AAAA (all zeros)
    println!("\n# Test 5: Insert AAAA (direct)");
    let mut kmer4 = Kmer::new(4, KmerMode::Direct);
    for _ in 0..4 {
        kmer4.insert(0); // A
    }
    println!("data_dir\t{:x}", kmer4.data_dir());
    println!("data_dir_top8\t{:x}", kmer4.data_dir() >> 56);

    // Test 6: Insert TTTT (all threes)
    println!("\n# Test 6: Insert TTTT (direct)");
    let mut kmer5 = Kmer::new(4, KmerMode::Direct);
    for _ in 0..4 {
        kmer5.insert(3); // T
    }
    println!("data_dir\t{:x}", kmer5.data_dir());
    println!("data_dir_top8\t{:x}", kmer5.data_dir() >> 56);

    // Test 7: Get symbol at position
    println!("\n# Test 7: Get symbols from ACGT");
    let mut kmer6 = Kmer::new(4, KmerMode::Direct);
    kmer6.insert(0); // A
    kmer6.insert(1); // C
    kmer6.insert(2); // G
    kmer6.insert(3); // T

    for i in 0..4 {
        println!("pos_{}\t{}", i, kmer6.get_symbol(i));
    }

    // Test 8: Test swap_dir_rc
    println!("\n# Test 8: Swap dir/rc");
    let mut kmer7 = Kmer::new(4, KmerMode::Canonical);
    kmer7.insert(0); // A
    kmer7.insert(1); // C
    kmer7.insert(2); // G
    kmer7.insert(3); // T

    let dir_before = kmer7.data_dir();
    let rc_before = kmer7.data_rc();
    println!("before_dir\t{dir_before:x}");
    println!("before_rc\t{rc_before:x}");

    kmer7.swap_dir_rc();

    let dir_after = kmer7.data_dir();
    let rc_after = kmer7.data_rc();
    println!("after_dir\t{dir_after:x}");
    println!("after_rc\t{rc_after:x}");

    // Test 9: Larger k-mer (k=8)
    println!("\n# Test 9: Larger k-mer (k=8)");
    let mut kmer8 = Kmer::new(8, KmerMode::Canonical);
    // Insert ACGTACGT
    kmer8.insert(0); // A
    kmer8.insert(1); // C
    kmer8.insert(2); // G
    kmer8.insert(3); // T
    kmer8.insert(0); // A
    kmer8.insert(1); // C
    kmer8.insert(2); // G
    kmer8.insert(3); // T

    println!("is_full\t{}", if kmer8.is_full() { 1 } else { 0 });
    println!("data_dir\t{:x}", kmer8.data_dir());
    println!("data_rc\t{:x}", kmer8.data_rc());
    println!("data_canonical\t{:x}", kmer8.data_canonical());

    // Test 10: Sliding window (insert with full k-mer)
    println!("\n# Test 10: Sliding window");
    let mut kmer9 = Kmer::new(4, KmerMode::Direct);
    // Insert ACGT
    kmer9.insert(0);
    kmer9.insert(1);
    kmer9.insert(2);
    kmer9.insert(3);
    println!("after_ACGT\t{:x}", kmer9.data_dir());

    // Insert another A (should shift out the first A)
    kmer9.insert(0);
    println!("after_CGTA\t{:x}", kmer9.data_dir());
}
