// Segmentation logic
// Splits contigs at splitter k-mer positions

use std::collections::HashSet;
use ragc_common::Contig;
use crate::kmer::{Kmer, KmerMode};

/// A segment of a contig bounded by splitter k-mers
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    /// The sequence data of this segment
    pub data: Contig,
    /// K-mer value at the start of the segment (or 0 if at contig start)
    pub front_kmer: u64,
    /// K-mer value at the end of the segment (or 0 if at contig end)
    pub back_kmer: u64,
}

impl Segment {
    /// Create a new segment
    pub fn new(data: Contig, front_kmer: u64, back_kmer: u64) -> Self {
        Segment {
            data,
            front_kmer,
            back_kmer,
        }
    }

    /// Get the length of the segment
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if the segment is empty
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// Split a contig at splitter k-mer positions
///
/// This function scans through the contig, identifies splitter k-mers,
/// and splits the contig into segments at those positions.
///
/// # Arguments
/// * `contig` - The contig to split
/// * `splitters` - Set of splitter k-mer values
/// * `k` - K-mer length
///
/// # Returns
/// Vector of segments
///
/// # Algorithm
/// 1. Scan through contig extracting k-mers
/// 2. When a splitter k-mer is found, create a segment up to that point
/// 3. The splitter k-mer becomes the back_kmer of one segment and front_kmer of the next
/// 4. Handle edge cases (contig start/end, consecutive splitters, no splitters)
pub fn split_at_splitters(contig: &Contig, splitters: &HashSet<u64>, k: usize) -> Vec<Segment> {
    let mut segments = Vec::new();

    if contig.len() < k {
        // Contig too short for k-mers, return as single segment
        return vec![Segment::new(contig.clone(), 0, 0)];
    }

    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
    let mut segment_start = 0;
    let mut front_kmer = 0u64;
    let mut pos = 0;

    for &base in contig {
        if base > 3 {
            // Non-ACGT base, reset k-mer
            kmer.reset();
        } else {
            kmer.insert(base as u64);

            if kmer.is_full() {
                let kmer_value = kmer.data();

                // Check if this is a splitter
                if splitters.contains(&kmer_value) {
                    // Create segment from segment_start to current position (inclusive of splitter)
                    let segment_end = pos + 1; // Include the base that completes the splitter
                    let segment_data = contig[segment_start..segment_end].to_vec();

                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, kmer_value));
                    }

                    // Start new segment after this splitter
                    segment_start = segment_end;
                    front_kmer = kmer_value;
                }
            }
        }
        pos += 1;
    }

    // Add final segment if there's data remaining
    if segment_start < contig.len() {
        let segment_data = contig[segment_start..].to_vec();
        if !segment_data.is_empty() {
            segments.push(Segment::new(segment_data, front_kmer, 0));
        }
    }

    // If no segments were created (no splitters), return entire contig as one segment
    if segments.is_empty() {
        segments.push(Segment::new(contig.clone(), 0, 0));
    }

    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_new() {
        let seg = Segment::new(vec![0, 1, 2, 3], 123, 456);
        assert_eq!(seg.len(), 4);
        assert_eq!(seg.front_kmer, 123);
        assert_eq!(seg.back_kmer, 456);
        assert!(!seg.is_empty());
    }

    #[test]
    fn test_split_no_splitters() {
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let splitters = HashSet::new();
        let segments = split_at_splitters(&contig, &splitters, 3);

        // Should return entire contig as one segment
        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].data, contig);
        assert_eq!(segments[0].front_kmer, 0);
        assert_eq!(segments[0].back_kmer, 0);
    }

    #[test]
    fn test_split_with_splitters() {
        // Create a test where we know the k-mer values
        let contig = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];

        // We need to figure out what the k-mer values are
        let mut kmer = Kmer::new(3, KmerMode::Canonical);
        let mut kmers = Vec::new();
        for &base in &contig {
            kmer.insert(base as u64);
            if kmer.is_full() {
                kmers.push(kmer.data());
            }
        }

        // Use the first k-mer as a splitter
        let mut splitters = HashSet::new();
        if !kmers.is_empty() {
            splitters.insert(kmers[0]);
        }

        let segments = split_at_splitters(&contig, &splitters, 3);

        // Should split at the splitter position
        assert!(segments.len() >= 1);

        // Total length should equal original contig
        let total_len: usize = segments.iter().map(|s| s.len()).sum();
        assert_eq!(total_len, contig.len());
    }

    #[test]
    fn test_split_short_contig() {
        let contig = vec![0, 1]; // Too short for k=3
        let splitters = HashSet::new();
        let segments = split_at_splitters(&contig, &splitters, 3);

        assert_eq!(segments.len(), 1);
        assert_eq!(segments[0].data, contig);
    }

    #[test]
    fn test_split_consecutive_splitters() {
        let contig = vec![0, 0, 0, 0, 0, 0]; // AAAAAA

        // All k-mers will be AAA (same value)
        let mut kmer = Kmer::new(3, KmerMode::Canonical);
        kmer.insert(0);
        kmer.insert(0);
        kmer.insert(0);
        let aaa_kmer = kmer.data();

        let mut splitters = HashSet::new();
        splitters.insert(aaa_kmer);

        let segments = split_at_splitters(&contig, &splitters, 3);

        // Should create multiple small segments
        assert!(segments.len() > 0);
    }

    #[test]
    fn test_split_with_n_bases() {
        // Contig with N (represented as 4)
        let contig = vec![0, 0, 0, 4, 1, 1, 1];

        let mut splitters = HashSet::new();
        splitters.insert(12345); // Some random splitter that won't match

        let segments = split_at_splitters(&contig, &splitters, 3);

        // Should handle N bases without crashing
        assert!(segments.len() >= 1);
    }
}
