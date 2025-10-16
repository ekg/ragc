// K-mer extraction and enumeration
// Rust equivalent of enumerate_kmers and related functions in agc_compressor.cpp

use crate::kmer::{Kmer, KmerMode};
use ragc_common::Contig;

/// Enumerate all canonical k-mers from a contig
///
/// This function extracts all k-mers from a contig sequence, returning their
/// canonical form (minimum of k-mer and its reverse complement).
///
/// # Arguments
/// * `contig` - The contig sequence (bases encoded as 0=A, 1=C, 2=G, 3=T)
/// * `k` - The k-mer length
///
/// # Returns
/// Vector of canonical k-mer values (as u64)
///
/// # Behavior
/// - If a non-ACGT base (value > 3) is encountered, the k-mer is reset
/// - Only full k-mers are returned (sequences shorter than k produce no k-mers)
/// - K-mers are returned in the order they appear in the contig
pub fn enumerate_kmers(contig: &Contig, k: usize) -> Vec<u64> {
    let mut vec = Vec::new();

    if contig.len() < k {
        return vec;
    }

    vec.reserve(contig.len() + 1 - k);

    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);

    for &base in contig {
        if base > 3 {
            // Non-ACGT base encountered, reset k-mer
            kmer.reset();
        } else {
            kmer.insert(base as u64);

            if kmer.is_full() {
                vec.push(kmer.data());
            }
        }
    }

    vec
}

/// Remove non-singleton k-mers from a sorted vector
///
/// This function removes all k-mers that appear more than once in the vector,
/// keeping only singletons (k-mers with exactly one occurrence).
///
/// # Arguments
/// * `vec` - Sorted vector of k-mer values
/// * `virtual_begin` - Starting index for processing (k-mers before this are preserved)
///
/// # Behavior
/// - Assumes the input vector is sorted
/// - K-mers before `virtual_begin` are not checked and remain in the output
/// - Modifies the vector in place, resizing it to contain only singletons
pub fn remove_non_singletons(vec: &mut Vec<u64>, virtual_begin: usize) {
    let mut curr_end = virtual_begin;
    let mut i = virtual_begin;

    while i < vec.len() {
        // Find the end of the run of equal values
        let mut j = i + 1;
        while j < vec.len() && vec[i] == vec[j] {
            j += 1;
        }

        // If this k-mer appears exactly once (singleton), keep it
        if i + 1 == j {
            vec[curr_end] = vec[i];
            curr_end += 1;
        }

        i = j;
    }

    vec.truncate(curr_end);
}

/// Find candidate k-mers (singletons) from a contig
///
/// This is a convenience function that combines k-mer enumeration,
/// sorting, and singleton filtering.
///
/// # Arguments
/// * `contig` - The contig sequence
/// * `k` - The k-mer length
///
/// # Returns
/// Sorted vector of singleton k-mers
pub fn find_candidate_kmers(contig: &Contig, k: usize) -> Vec<u64> {
    let mut kmers = enumerate_kmers(contig, k);
    kmers.sort_unstable();
    remove_non_singletons(&mut kmers, 0);
    kmers
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enumerate_kmers_simple() {
        // ACGT = [0, 1, 2, 3]
        let contig = vec![0, 1, 2, 3]; // ACGT
        let kmers = enumerate_kmers(&contig, 3);

        // Should get 2 k-mers: ACG and CGT
        assert_eq!(kmers.len(), 2);
    }

    #[test]
    fn test_enumerate_kmers_with_n() {
        // ACN (N=4) GT - should reset at N
        let contig = vec![0, 1, 4, 2, 3]; // AC[N]GT
        let kmers = enumerate_kmers(&contig, 3);

        // Should only get GT* (but we need 3 bases after N, so none)
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_enumerate_kmers_short_sequence() {
        let contig = vec![0, 1]; // AC (too short for k=3)
        let kmers = enumerate_kmers(&contig, 3);
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_remove_non_singletons() {
        let mut vec = vec![1, 2, 2, 3, 3, 3, 4, 5, 5, 6];
        remove_non_singletons(&mut vec, 0);

        // Only 1, 4, 6 are singletons
        assert_eq!(vec, vec![1, 4, 6]);
    }

    #[test]
    fn test_remove_non_singletons_with_virtual_begin() {
        let mut vec = vec![1, 1, 2, 3, 3, 4, 5, 5];
        remove_non_singletons(&mut vec, 2);

        // First 2 elements preserved, then only 2, 4 are singletons
        assert_eq!(vec, vec![1, 1, 2, 4]);
    }

    #[test]
    fn test_find_candidate_kmers() {
        // Create a sequence with some repeated k-mers
        // ACGTACGT = [0,1,2,3,0,1,2,3]
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let candidates = find_candidate_kmers(&contig, 3);

        // All 3-mers appear at least twice, so no singletons
        assert_eq!(candidates.len(), 0);
    }
}
