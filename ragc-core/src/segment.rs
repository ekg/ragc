// Segmentation logic
// Splits contigs at splitter k-mer positions

use crate::kmer::{Kmer, KmerMode};
use ragc_common::Contig;
use std::collections::HashSet;

/// Missing k-mer sentinel value (matches C++ AGC's kmer_t(-1) = u64::MAX)
/// Used when a segment doesn't have a front or back k-mer (e.g., at contig boundaries)
pub const MISSING_KMER: u64 = u64::MAX;

/// A segment of a contig bounded by splitter k-mers
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    /// The sequence data of this segment
    pub data: Contig,
    /// K-mer value at the start of the segment (or MISSING_KMER if at contig start)
    pub front_kmer: u64,
    /// K-mer value at the end of the segment (or MISSING_KMER if at contig end)
    pub back_kmer: u64,
    /// Whether front k-mer is in direct orientation (for C++ AGC's is_dir_oriented())
    pub front_kmer_is_dir: bool,
    /// Whether back k-mer is in direct orientation (for C++ AGC's is_dir_oriented())
    pub back_kmer_is_dir: bool,
}

impl Segment {
    /// Create a new segment
    pub fn new(data: Contig, front_kmer: u64, back_kmer: u64, front_kmer_is_dir: bool, back_kmer_is_dir: bool) -> Self {
        Segment {
            data,
            front_kmer,
            back_kmer,
            front_kmer_is_dir,
            back_kmer_is_dir,
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
pub fn split_at_splitters_with_size(
    contig: &Contig,
    splitters: &HashSet<u64>,
    k: usize,
    _min_segment_size: usize,
) -> Vec<Segment> {
    let debug = crate::env_cache::debug_segment_coverage();
    if debug {
        eprintln!("\n=== SEGMENTATION START: contig_len={} k={} ===", contig.len(), k);
    }

    let mut segments = Vec::new();

    if contig.len() < k {
        // Contig too short for k-mers, return as single segment
        if debug {
            eprintln!("Contig too short for k-mers, returning as single segment");
        }
        return vec![Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER, false, false)];
    }

    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
    let mut segment_start = 0;
    let mut front_kmer = MISSING_KMER;
    let mut front_kmer_is_dir = false;

    // Track recent k-mers for end-of-contig handling
    // C++ AGC doesn't limit this - it accumulates all k-mers since last split
    // Store (position, kmer_value, is_dir_oriented) to match C++ AGC's orientation logic
    let mut recent_kmers: Vec<(usize, u64, bool)> = Vec::new();

    for (pos, &base) in contig.iter().enumerate() {
        if base > 3 {
            // Non-ACGT base, reset k-mer
            kmer.reset();
            recent_kmers.clear();
        } else {
            kmer.insert(base as u64);

            if kmer.is_full() {
                let kmer_value = kmer.data();
                let is_dir = kmer.is_dir_oriented();
                recent_kmers.push((pos, kmer_value, is_dir));

                // CRITICAL FIX: C++ AGC splits at EVERY splitter occurrence!
                // The distance check (current_len >= min_segment_size) only happens
                // during SPLITTER FINDING to select which k-mers become splitters.
                // During SEGMENTATION, we split at every occurrence without distance check.

                // DEBUG: Trace positions near end of contig
                if crate::env_cache::debug_endpos() && pos + 50 >= contig.len() {
                    let is_splitter = splitters.contains(&kmer_value);
                    let bytes_left = contig.len() - (pos + 1);
                    eprintln!("ENDPOS_TRACE: pos={} kmer={:#x} is_splitter={} bytes_left={} k={} contig_len={}",
                              pos, kmer_value, is_splitter, bytes_left, k, contig.len());
                }

                if splitters.contains(&kmer_value) {
                    // Comprehensive split logging for debugging
                    if crate::env_cache::trace_all_splits() {
                        let segment_len = (pos + 1) - segment_start;
                        eprintln!("RAGC_SPLIT: pos={} kmer={} segment_start={} segment_len={} contig_len={}",
                                  pos, kmer_value, segment_start, segment_len, contig.len());
                    }
                    // Use this as a splitter
                    let segment_end = pos + 1;
                    let segment_data = contig[segment_start..segment_end].to_vec();

                    if !segment_data.is_empty() {
                        // Determine front and back k-mers and their orientations
                        let (seg_front, seg_back, seg_front_is_dir, seg_back_is_dir) = if front_kmer == MISSING_KMER {
                            // First segment of contig - ALWAYS back k-mer only
                            // The splitter is at the END of the segment, so it's the back k-mer
                            (MISSING_KMER, kmer_value, false, is_dir)
                        } else {
                            // Normal case: both k-mers present
                            (front_kmer, kmer_value, front_kmer_is_dir, is_dir)
                        };
                        // Note: Contig name is logged at call site, this is just position info
                        if debug {
                            eprintln!("  MAIN_LOOP_SPLIT: segment=[{}..{}) len={}", segment_start, segment_end, segment_data.len());
                        }
                        #[cfg(feature = "verbose_debug")]
                        if crate::env_cache::debug_overlap() {
                            let first_5: Vec<u8> = segment_data.iter().take(5).copied().collect();
                            let last_5: Vec<u8> = segment_data.iter().rev().take(5).rev().copied().collect();
                            eprintln!("RAGC_SEG_SPLIT: pos={} splitter={} segment=[{}..{}) len={} front={} back={} first_5={:?} last_5={:?}",
                                pos, kmer_value, segment_start, segment_end, segment_data.len(),
                                if seg_front == MISSING_KMER { "MISSING".to_string() } else { seg_front.to_string() },
                                if seg_back == MISSING_KMER { "MISSING".to_string() } else { seg_back.to_string() },
                                first_5, last_5);
                        } else {
                            #[cfg(feature = "verbose_debug")]
                            eprintln!("RAGC_SEG_SPLIT: pos={} splitter={} segment=[{}..{}) len={} front={} back={}",
                                pos, kmer_value, segment_start, segment_end, segment_data.len(),
                                if seg_front == MISSING_KMER { "MISSING".to_string() } else { seg_front.to_string() },
                                if seg_back == MISSING_KMER { "MISSING".to_string() } else { seg_back.to_string() });
                        }
                        segments.push(Segment::new(segment_data, seg_front, seg_back, seg_front_is_dir, seg_back_is_dir));
                    }

                    // Reset for next segment
                    // CRITICAL: Create k-byte overlap to match C++ AGC behavior
                    // Each segment must include the FULL k-mer at the start
                    let new_start = (pos + 1).saturating_sub(k);
                    if debug {
                        eprintln!("  Setting segment_start: {} -> {} (overlap of {} bytes)", segment_start, new_start, (pos + 1) - new_start);
                    }
                    segment_start = new_start;
                    front_kmer = kmer_value;
                    front_kmer_is_dir = is_dir;
                    recent_kmers.clear();
                    kmer.reset();
                }
            }
        }
    }

    // End-of-contig handling: C++ AGC segmentation does not perform any
    // backward search for a last-minute splitter. It simply splits at every
    // occurrence encountered in the main loop and then emits the final segment.
    // We therefore intentionally do nothing here to match C++ behavior exactly.

    // Add any remaining data as final segment
    // This will either be:
    // - The entire remainder if no suitable splitter was found, OR
    // - A segment >= k bytes if we split at a splitter that left enough room
    if debug {
        eprintln!("\n=== FINAL SEGMENT ===");
        eprintln!("  segment_start={}, contig.len()={}", segment_start, contig.len());
    }
    if segment_start < contig.len() {
        let segment_data = contig[segment_start..].to_vec();
        if !segment_data.is_empty() {
            // Final segment k-mer extraction
            let (final_front, final_back, final_front_is_dir, final_back_is_dir) = if front_kmer == MISSING_KMER {
                // No splitters found in entire contig - this is an orphan segment
                // C++ AGC uses (MISSING_KMER, MISSING_KMER) for orphan segments,
                // which routes them to raw group 0 (see distribute_segments in agc_compressor.h)
                // FIX 15: Match C++ AGC behavior - use MISSING for both k-mers
                (MISSING_KMER, MISSING_KMER, false, false)
            } else {
                // Front k-mer present (from last splitter)
                // C++ AGC uses MISSING_KMER for back k-mer of final segments
                // (see agc_compressor.cpp line 2188)
                (front_kmer, MISSING_KMER, front_kmer_is_dir, false)
            };
            if debug {
                eprintln!("  FINAL: segment=[{}..{}) len={}", segment_start, contig.len(), segment_data.len());
            }
            #[cfg(feature = "verbose_debug")]
            eprintln!("RAGC_SEG_FINAL: segment=[{}..{}) len={} front={} back={}",
                segment_start, contig.len(), segment_data.len(),
                if final_front == MISSING_KMER { "MISSING".to_string() } else { final_front.to_string() },
                if final_back == MISSING_KMER { "MISSING".to_string() } else { final_back.to_string() });
            // Debug logging for Case 3 is_dir investigation
            if crate::env_cache::debug_is_dir() && final_back == MISSING_KMER && final_front != MISSING_KMER {
                eprintln!("RAGC_FINAL_SEG_IS_DIR: front_kmer={} front_kmer_is_dir={}", final_front, final_front_is_dir);
            }
            segments.push(Segment::new(segment_data, final_front, final_back, final_front_is_dir, final_back_is_dir));
        }
    }

    // If no segments were created, return entire contig as one segment
    if segments.is_empty() {
        if debug {
            eprintln!("  NO_SPLIT: Returning entire contig as single segment");
        }
        #[cfg(feature = "verbose_debug")]
        eprintln!("RAGC_SEG_NOSPLIT: len={} front=MISSING back=MISSING", contig.len());
        segments.push(Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER, false, false));
    }

    // Summary
    if debug {
        eprintln!("\n=== SEGMENTATION SUMMARY ===");
        eprintln!("  Original contig length: {}", contig.len());
        eprintln!("  Number of segments: {}", segments.len());
        let total_segment_bytes: usize = segments.iter().map(|s| s.len()).sum();
        eprintln!("  Total segment bytes: {}", total_segment_bytes);

        // Calculate expected reconstructed size
        let expected_size = if segments.is_empty() {
            0
        } else if segments.len() == 1 {
            segments[0].len()
        } else {
            // First segment contributes all bytes, subsequent segments skip (k-1) overlap
            segments[0].len() + segments[1..].iter().map(|s| s.len().saturating_sub(k - 1)).sum::<usize>()
        };
        eprintln!("  Expected reconstructed size: {} (with {} overlaps of {} bytes)",
                  expected_size, segments.len().saturating_sub(1), k - 1);

        if expected_size != contig.len() {
            eprintln!("  ⚠️  SIZE MISMATCH: Expected {} but contig is {} (diff: {})",
                      expected_size, contig.len(), contig.len() as i64 - expected_size as i64);
        } else {
            eprintln!("  ✓ Size matches!");
        }

        // Show segment ranges to identify gaps
        eprintln!("\n  Segment coverage:");
        for (i, seg) in segments.iter().enumerate() {
            eprintln!("    Segment {}: len={}", i, seg.len());
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
        return vec![Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER, false, false)];
    }

    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
    let mut segment_start = 0;
    let mut front_kmer = MISSING_KMER;
    let mut front_kmer_is_dir = false;

    for (pos, &base) in contig.iter().enumerate() {
        if base > 3 {
            // Non-ACGT base, reset k-mer
            kmer.reset();
        } else {
            kmer.insert(base as u64);

            if kmer.is_full() {
                let kmer_value = kmer.data();
                let is_dir = kmer.is_dir_oriented();

                // Check if this is a splitter
                if splitters.contains(&kmer_value) {
                    // Create segment from segment_start to current position (inclusive of splitter)
                    let segment_end = pos + 1; // Include the base that completes the splitter
                    let segment_data = contig[segment_start..segment_end].to_vec();

                    if !segment_data.is_empty() {
                        #[cfg(feature = "verbose_debug")]
                        eprintln!("RAGC_SEG_SPLIT: pos={} splitter={} segment=[{}..{}) len={} front={} back={}",
                            pos, kmer_value, segment_start, segment_end, segment_data.len(),
                            if front_kmer == MISSING_KMER { "MISSING".to_string() } else { front_kmer.to_string() },
                            kmer_value);
                        segments.push(Segment::new(segment_data, front_kmer, kmer_value, front_kmer_is_dir, is_dir));
                    }

                    // Start new segment with k-base overlap
                    segment_start = (pos + 1).saturating_sub(k);
                    front_kmer = kmer_value;
                    front_kmer_is_dir = is_dir;
                }
            }
        }
    }

    // Add final segment if there's data remaining
    if segment_start < contig.len() {
        let segment_data = contig[segment_start..].to_vec();
        if !segment_data.is_empty() {
            #[cfg(feature = "verbose_debug")]
            eprintln!("RAGC_SEG_FINAL: segment=[{}..{}) len={} front={} back=MISSING",
                segment_start, contig.len(), segment_data.len(),
                if front_kmer == MISSING_KMER { "MISSING".to_string() } else { front_kmer.to_string() });
            segments.push(Segment::new(segment_data, front_kmer, MISSING_KMER, front_kmer_is_dir, false));
        }
    }

    // If no segments were created (no splitters), return entire contig as one segment
    if segments.is_empty() {
        #[cfg(feature = "verbose_debug")]
        eprintln!("RAGC_SEG_NOSPLIT: len={} front=MISSING back=MISSING", contig.len());
        segments.push(Segment::new(contig.clone(), MISSING_KMER, MISSING_KMER, false, false));
    }

    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_new() {
        let seg = Segment::new(vec![0, 1, 2, 3], 123, 456, true, false);
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
        assert_eq!(segments[0].front_kmer, MISSING_KMER);
        assert_eq!(segments[0].back_kmer, MISSING_KMER);
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

        // Verify reconstructed length (accounting for k-base overlaps)
        // First segment contributes all its bytes, subsequent segments skip first k bytes
        let k = 3;
        let reconstructed_len: usize = if segments.is_empty() {
            0
        } else {
            segments[0].len()
                + segments[1..]
                    .iter()
                    .map(|s| s.len().saturating_sub(k))
                    .sum::<usize>()
        };
        assert_eq!(reconstructed_len, contig.len());
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
