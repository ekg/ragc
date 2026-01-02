// FFI implementation for finding splitters in a contig
// Replaces C++ AGC's find_splitters_in_contig() function

use crate::kmer::{Kmer, KmerMode};
use ragc_common::hash::MurMur64Hash;
use std::slice;

/// Result structure for find_splitters_in_contig
/// Contains both splitters and fallback minimizer mappings
#[repr(C)]
pub struct FindSplittersResult {
    /// Array of splitter k-mer values
    pub splitters_data: *mut u64,
    pub splitters_len: usize,
    /// Flattened array of fallback mappings (4 u64s per entry)
    /// Each entry: [prev_splitter, next_splitter, fallback_kmer, is_dir_oriented]
    pub fallbacks_data: *mut u64,
    pub fallbacks_len: usize,
}

/// Find splitters in a contig
///
/// Matches C++ AGC's find_splitters_in_contig() algorithm:
/// 1. Iterate through contig building k-mers
/// 2. Reset on non-ACGT bases (> 3)
/// 3. Track fallback minimizers (filtered by hash threshold)
/// 4. At segment boundaries, find the best splitter k-mer
/// 5. Return both splitters and fallback mappings
///
/// # Parameters
/// - contig_data: Pointer to contig sequence (numeric encoding: A=0, C=1, G=2, T=3)
/// - contig_len: Length of contig
/// - k: K-mer length
/// - segment_size: Target segment size
/// - candidate_kmers_ptr: Pointer to sorted candidate splitter k-mers
/// - candidate_kmers_len: Length of candidate k-mers array
/// - fallback_threshold: Hash threshold for fallback filtering (0 = disabled)
///
/// # Returns
/// FindSplittersResult containing splitters and fallback mappings
/// Must be freed with ragc_free_find_splitters_result()
///
/// # Safety
/// - All pointers must be valid for their specified lengths
/// - candidate_kmers must be sorted (for binary search)
#[no_mangle]
pub extern "C" fn ragc_find_splitters_in_contig(
    contig_data: *const u8,
    contig_len: usize,
    k: u32,
    segment_size: u64,
    candidate_kmers_ptr: *const u64,
    candidate_kmers_len: usize,
    fallback_threshold: u64,
) -> FindSplittersResult {
    unsafe {
        let contig = slice::from_raw_parts(contig_data, contig_len);
        let candidate_kmers = slice::from_raw_parts(candidate_kmers_ptr, candidate_kmers_len);

        eprintln!("[RUST FFI] find_splitters_in_contig called:");
        eprintln!("[RUST FFI]   contig_len: {}", contig_len);
        eprintln!("[RUST FFI]   k: {}", k);
        eprintln!("[RUST FFI]   segment_size: {}", segment_size);
        eprintln!("[RUST FFI]   candidate_kmers_len: {}", candidate_kmers_len);
        if candidate_kmers_len > 0 {
            eprintln!(
                "[RUST FFI]   first 5 candidates: {:?}",
                &candidate_kmers[..candidate_kmers_len.min(5)]
            );
        }

        // Initialize output vectors
        let mut splitters = Vec::new();
        let mut fallbacks: Vec<[u64; 4]> = Vec::new();

        // Initialization to large value to add 1st candidate k-mer
        let mut current_len = segment_size;
        let mut v_recent_kmers = Vec::new();
        let mut kmer = Kmer::new(k, KmerMode::Canonical);

        let mut prev_splitter = u64::MAX; // ~0ull in C++
        let mut fallback_kmers_in_segment: Vec<(u64, bool)> = Vec::new();

        const RND: u64 = 0xD73F8BF11046C40E;

        kmer.reset();

        let mut checked_count = 0;
        let mut found_count = 0;

        for &base in contig {
            // Reset k-mer on non-ACGT bases
            if base > 3 {
                kmer.reset();
            } else {
                kmer.insert(base as u64);

                if kmer.is_full() {
                    let kmer_value = kmer.data();
                    v_recent_kmers.push(kmer_value);

                    // Check if k-mer should be tracked as fallback
                    // C++: fallback_filter(kmer.data()) && kmer.data_dir() != kmer.data_rc()
                    let passes_fallback_filter = if fallback_threshold == 0 {
                        false
                    } else {
                        (MurMur64Hash::hash(kmer_value) ^ RND) < fallback_threshold
                    };

                    // Only track symmetric k-mers (where orientation is clear)
                    if passes_fallback_filter && kmer.data_dir() != kmer.data_rc() {
                        fallback_kmers_in_segment.push((kmer_value, kmer.is_dir_oriented()));
                    }

                    // Check if we've reached a segment boundary
                    if current_len >= segment_size {
                        // Check if this k-mer is a splitter
                        checked_count += 1;
                        if is_splitter(kmer_value, candidate_kmers) {
                            found_count += 1;
                            if found_count <= 3 {
                                eprintln!("[RUST FFI] Found splitter #{}: kmer_value=0x{:016x} at current_len={}",
                                    found_count, kmer_value, current_len);
                                eprintln!(
                                    "[RUST FFI]   kmer_dir=0x{:016x}, kmer_rc=0x{:016x}",
                                    kmer.data_dir(),
                                    kmer.data_rc()
                                );
                            }
                            splitters.push(kmer_value);

                            // Add fallback mappings for this segment
                            for &(fallback_kmer, is_dir) in &fallback_kmers_in_segment {
                                fallbacks.push([
                                    prev_splitter,
                                    kmer_value,
                                    fallback_kmer,
                                    is_dir as u64,
                                ]);
                            }

                            fallback_kmers_in_segment.clear();
                            prev_splitter = kmer_value;
                            current_len = 0;
                            kmer.reset();
                            v_recent_kmers.clear();
                        }
                    }
                }
            }

            current_len += 1;
        }

        // Try to add the rightmost candidate k-mer
        // Iterate backwards through recent k-mers to find a splitter
        for &kmer_value in v_recent_kmers.iter().rev() {
            if is_splitter(kmer_value, candidate_kmers) {
                found_count += 1;
                eprintln!(
                    "[RUST FFI] Found rightmost splitter: kmer_value=0x{:016x}",
                    kmer_value
                );
                splitters.push(kmer_value);
                for &(fallback_kmer, is_dir) in &fallback_kmers_in_segment {
                    fallbacks.push([prev_splitter, kmer_value, fallback_kmer, is_dir as u64]);
                }
                break;
            }
        }

        eprintln!(
            "[RUST FFI] Total: checked {} boundary positions, found {} splitters",
            checked_count, found_count
        );

        // Flatten fallbacks array for FFI
        let mut fallbacks_flat: Vec<u64> = Vec::with_capacity(fallbacks.len() * 4);
        for entry in fallbacks {
            fallbacks_flat.extend_from_slice(&entry);
        }

        // Prepare result
        let splitters_ptr = splitters.as_mut_ptr();
        let splitters_len = splitters.len();
        std::mem::forget(splitters);

        let fallbacks_ptr = fallbacks_flat.as_mut_ptr();
        let fallbacks_len = fallbacks_flat.len() / 4; // Number of entries, not total elements
        std::mem::forget(fallbacks_flat);

        FindSplittersResult {
            splitters_data: splitters_ptr,
            splitters_len,
            fallbacks_data: fallbacks_ptr,
            fallbacks_len,
        }
    }
}

/// Free memory allocated by ragc_find_splitters_in_contig()
///
/// # Safety
/// - Must only be called once per FindSplittersResult
/// - Result must have been created by ragc_find_splitters_in_contig()
#[no_mangle]
pub extern "C" fn ragc_free_find_splitters_result(result: FindSplittersResult) {
    unsafe {
        if !result.splitters_data.is_null() && result.splitters_len > 0 {
            let _ = Vec::from_raw_parts(
                result.splitters_data,
                result.splitters_len,
                result.splitters_len,
            );
        }

        if !result.fallbacks_data.is_null() && result.fallbacks_len > 0 {
            // Reconstruct flattened array (4 u64s per entry)
            let flat_len = result.fallbacks_len * 4;
            let _ = Vec::from_raw_parts(result.fallbacks_data, flat_len, flat_len);
        }
    }
}

/// Check if a k-mer is a splitter using binary search
#[inline]
fn is_splitter(kmer_value: u64, candidate_kmers: &[u64]) -> bool {
    candidate_kmers.binary_search(&kmer_value).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_splitters_simple() {
        // Simple test case: ACGTACGT sequence
        let contig = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let k = 3;
        let segment_size = 5;

        // No candidate k-mers - should find no splitters
        let candidate_kmers: Vec<u64> = vec![];

        let result = ragc_find_splitters_in_contig(
            contig.as_ptr(),
            contig.len(),
            k,
            segment_size,
            candidate_kmers.as_ptr(),
            candidate_kmers.len(),
            0, // No fallback filtering
        );

        unsafe {
            let splitters = slice::from_raw_parts(result.splitters_data, result.splitters_len);
            assert_eq!(splitters.len(), 0);
        }

        ragc_free_find_splitters_result(result);
    }

    #[test]
    fn test_find_splitters_with_candidates() {
        // ACGTACGT sequence
        let contig = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let k = 3;
        let segment_size = 3;

        // Extract actual k-mers from the sequence to use as candidates
        let mut kmer = Kmer::new(k, KmerMode::Canonical);
        let mut candidate_kmers = Vec::new();
        for &base in &contig {
            if base <= 3 {
                kmer.insert(base as u64);
                if kmer.is_full() {
                    candidate_kmers.push(kmer.data());
                }
            }
        }
        candidate_kmers.sort_unstable();
        candidate_kmers.dedup();

        let result = ragc_find_splitters_in_contig(
            contig.as_ptr(),
            contig.len(),
            k,
            segment_size,
            candidate_kmers.as_ptr(),
            candidate_kmers.len(),
            0, // No fallback filtering
        );

        unsafe {
            let splitters = slice::from_raw_parts(result.splitters_data, result.splitters_len);
            // Should find at least one splitter
            assert!(splitters.len() > 0);
        }

        ragc_free_find_splitters_result(result);
    }

    #[test]
    fn test_fallback_filtering() {
        // Test with fallback threshold enabled
        let contig = vec![0u8, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let k = 3;
        let segment_size = 3;

        let mut kmer = Kmer::new(k, KmerMode::Canonical);
        let mut candidate_kmers = Vec::new();
        for &base in &contig {
            if base <= 3 {
                kmer.insert(base as u64);
                if kmer.is_full() {
                    candidate_kmers.push(kmer.data());
                }
            }
        }
        candidate_kmers.sort_unstable();
        candidate_kmers.dedup();

        // Use a high threshold to capture many k-mers
        let fallback_threshold = u64::MAX / 2;

        let result = ragc_find_splitters_in_contig(
            contig.as_ptr(),
            contig.len(),
            k,
            segment_size,
            candidate_kmers.as_ptr(),
            candidate_kmers.len(),
            fallback_threshold,
        );

        // Check that we have fallback mappings
        // May or may not have fallbacks depending on hash function
        // Just verify the structure is valid (no crash when freeing)
        let _fallbacks_len = result.fallbacks_len;

        ragc_free_find_splitters_result(result);
    }
}
