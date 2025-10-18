// Segmentation logic
// Splits contigs at splitter k-mer positions

use crate::kmer::{Kmer, KmerMode};
use ragc_common::Contig;
use std::collections::HashSet;

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
/// This implements the C++ AGC segment splitting algorithm.
/// The splitters passed should be the ACTUALLY USED splitters from the reference
/// (not all candidates), ensuring all genomes split at the same positions.
///
/// # Arguments
/// * `contig` - The contig to split
/// * `splitters` - Set of splitter k-mer values (actually used on reference)
/// * `k` - K-mer length
/// * `_min_segment_size` - Unused (kept for API compatibility, see note below)
///
/// # Returns
/// Vector of segments
///
/// # Algorithm (matching C++ AGC)
/// 1. Scan through contig
/// 2. Split at EVERY occurrence of a splitter k-mer (no distance check!)
/// 3. Keep track of recent k-mers for end-of-contig handling
/// 4. At contig end, look backward through recent k-mers to find rightmost splitter
///
/// # Note on min_segment_size
/// The distance check only happens during SPLITTER FINDING (in splitters.rs) to select
/// which k-mers become splitters. During SEGMENTATION (this function), C++ AGC splits
/// at every occurrence of those splitters WITHOUT checking distance. This is the key
/// difference that was causing 2.3x worse compression!
pub fn split_at_splitters_with_size(contig: &Contig, splitters: &HashSet<u64>, k: usize, _min_segment_size: usize) -> Vec<Segment> {
    let mut segments = Vec::new();

    if contig.len() < k {
        // Contig too short for k-mers, return as single segment
        return vec![Segment::new(contig.clone(), 0, 0)];
    }

    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
    let mut segment_start = 0;
    let mut front_kmer = 0u64;

    // Track recent k-mers for end-of-contig handling
    // C++ AGC doesn't limit this - it accumulates all k-mers since last split
    let mut recent_kmers: Vec<(usize, u64)> = Vec::new(); // (position, kmer_value)

    for (pos, &base) in contig.iter().enumerate() {
        if base > 3 {
            // Non-ACGT base, reset k-mer
            kmer.reset();
            recent_kmers.clear();
        } else {
            kmer.insert(base as u64);

            if kmer.is_full() {
                let kmer_value = kmer.data();
                recent_kmers.push((pos, kmer_value));

                // CRITICAL FIX: C++ AGC splits at EVERY splitter occurrence!
                // The distance check (current_len >= min_segment_size) only happens
                // during SPLITTER FINDING to select which k-mers become splitters.
                // During SEGMENTATION, we split at every occurrence without distance check.
                if splitters.contains(&kmer_value) {
                    // Use this as a splitter
                    let segment_end = pos + 1;
                    let segment_data = contig[segment_start..segment_end].to_vec();

                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, kmer_value));
                    }

                    // Reset for next segment
                    // CRITICAL: Create k-base overlap so decompressor can skip first k bases
                    // The k-mer ends at position pos, occupies [pos-k+1, pos]
                    // Next segment should start at pos-k+1 to create k-base overlap
                    segment_start = (pos + 1).saturating_sub(k);
                    front_kmer = kmer_value;
                    recent_kmers.clear();
                    kmer.reset();
                }
            }
        }
    }

    // At end of contig, look backward through recent k-mers to find rightmost splitter
    // CRITICAL: Only split if the remaining data after the splitter will be >= k bytes
    // This ensures C++ AGC can skip k overlap bytes without hitting corruption
    for (pos, kmer_value) in recent_kmers.iter().rev() {
        if splitters.contains(kmer_value) {
            let segment_end = pos + 1;
            let remaining_after = contig.len() - segment_end;

            // Only split here if remaining data is > k bytes
            // (We need > k, not >= k, because after creating k-base overlap,
            // the final segment must still have > k bytes for C++ AGC decompressor)
            if remaining_after > k {
                if segment_end > segment_start {
                    let segment_data = contig[segment_start..segment_end].to_vec();
                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, *kmer_value));
                        // Create k-base overlap for next segment
                        segment_start = (pos + 1).saturating_sub(k);
                        front_kmer = *kmer_value;
                    }
                }
                break;
            } else if remaining_after == 0 {
                // Splitter is exactly at contig end - include it in current segment
                // and don't update segment_start (no next segment to create)
                if segment_end > segment_start {
                    let segment_data = contig[segment_start..segment_end].to_vec();
                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, *kmer_value));
                        // Mark that we've consumed the entire contig
                        segment_start = contig.len();
                    }
                }
                break;
            }
            // Otherwise, continue looking for an earlier splitter that leaves enough room
        }
    }

    // Add any remaining data as final segment
    // This will either be:
    // - The entire remainder if no suitable splitter was found, OR
    // - A segment >= k bytes if we split at a splitter that left enough room
    if segment_start < contig.len() {
        let segment_data = contig[segment_start..].to_vec();
        if !segment_data.is_empty() {
            segments.push(Segment::new(segment_data, front_kmer, 0));
        }
    }

    // If no segments were created, return entire contig as one segment
    if segments.is_empty() {
        segments.push(Segment::new(contig.clone(), 0, 0));
    }

    // CRITICAL FIX: Merge final segment with previous if it's too short for overlap
    // Segments after the first must be >= k bytes to handle the k-base overlap
    if segments.len() >= 2 {
        let last_idx = segments.len() - 1;
        if segments[last_idx].data.len() < k {
            // Merge last two segments
            let last_seg = segments.pop().unwrap();
            let second_last = segments.last_mut().unwrap();

            // Append last segment data to second-last
            second_last.data.extend_from_slice(&last_seg.data);
            // Keep the back_kmer from the merged segment (should be 0 for final segment)
            second_last.back_kmer = last_seg.back_kmer;
        }
    }

    segments
}

/// Split a contig at splitter k-mer positions (no minimum size constraint)
///
/// This is the old behavior - splits at every splitter regardless of segment size.
/// Kept for backwards compatibility but not recommended for use.
pub fn split_at_splitters(contig: &Contig, splitters: &HashSet<u64>, k: usize) -> Vec<Segment> {
    let mut segments = Vec::new();

    if contig.len() < k {
        // Contig too short for k-mers, return as single segment
        return vec![Segment::new(contig.clone(), 0, 0)];
    }

    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
    let mut segment_start = 0;
    let mut front_kmer = 0u64;

    for (pos, &base) in contig.iter().enumerate() {
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

                    // Start new segment with k-base overlap
                    segment_start = (pos + 1).saturating_sub(k);
                    front_kmer = kmer_value;
                }
            }
        }
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
        assert!(!segments.is_empty());

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
        assert!(!segments.is_empty());
    }

    #[test]
    fn test_split_with_n_bases() {
        // Contig with N (represented as 4)
        let contig = vec![0, 0, 0, 4, 1, 1, 1];

        let mut splitters = HashSet::new();
        splitters.insert(12345); // Some random splitter that won't match

        let segments = split_at_splitters(&contig, &splitters, 3);

        // Should handle N bases without crashing
        assert!(!segments.is_empty());
    }
}
