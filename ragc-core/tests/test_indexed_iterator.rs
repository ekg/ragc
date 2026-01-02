#![allow(clippy::all)]
#![allow(unexpected_cfgs)]
#![allow(unused_imports)]
#[cfg(feature = "indexed-fasta")]
#[cfg(test)]
mod tests {
    use ragc_core::contig_iterator::{ContigIterator, IndexedPansnFileIterator};
    use std::path::Path;

    #[test]
    #[ignore] // Only run explicitly, requires test data
    fn test_indexed_iterator_basic() {
        let test_file = Path::new("/home/erik/scrapy/yeast235_bgzip.fa.gz");

        if !test_file.exists() {
            eprintln!("Test file not found: {}", test_file.display());
            eprintln!("Skipping test");
            return;
        }

        let mut iterator = IndexedPansnFileIterator::new(test_file)
            .expect("Failed to create IndexedPansnFileIterator");

        // Fetch first 10 contigs
        let mut count = 0;
        for i in 0..10 {
            match iterator.next_contig().expect("Error fetching contig") {
                Some((sample, header, contig)) => {
                    println!(
                        "Contig {}: Sample={}, Header={}, Length={}",
                        i + 1,
                        sample,
                        header,
                        contig.len()
                    );
                    assert!(!sample.is_empty(), "Sample name should not be empty");
                    assert!(!header.is_empty(), "Header should not be empty");
                    assert!(!contig.is_empty(), "Contig should have non-zero length");
                    count += 1;
                }
                None => {
                    panic!("Unexpected end of contigs at iteration {i}");
                }
            }
        }

        assert_eq!(count, 10, "Should have fetched 10 contigs");
        println!("\nIndexed iterator test successful!");
    }

    #[test]
    #[ignore] // Only run explicitly, requires test data
    fn test_indexed_iterator_sample_ordering() {
        let test_file = Path::new("/home/erik/scrapy/yeast235_bgzip.fa.gz");

        if !test_file.exists() {
            eprintln!("Test file not found: {}", test_file.display());
            eprintln!("Skipping test");
            return;
        }

        let mut iterator = IndexedPansnFileIterator::new(test_file)
            .expect("Failed to create IndexedPansnFileIterator");

        // Track samples to verify they appear in contiguous blocks
        let mut last_sample = String::new();
        let mut sample_transitions = Vec::new();
        let mut total_contigs = 0;

        // Read first 50 contigs to check sample grouping
        for _i in 0..50 {
            match iterator.next_contig().expect("Error fetching contig") {
                Some((sample, _header, _contig)) => {
                    if sample != last_sample {
                        sample_transitions.push(sample.clone());
                        last_sample = sample;
                    }
                    total_contigs += 1;
                }
                None => break,
            }
        }

        println!(
            "Read {} contigs with {} sample transitions",
            total_contigs,
            sample_transitions.len()
        );
        println!("Sample order: {sample_transitions:?}");

        // With proper sample grouping, we should have fewer transitions than total contigs
        assert!(
            sample_transitions.len() < total_contigs,
            "Samples should be grouped contiguously"
        );
    }
}
