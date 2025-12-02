// FFI helpers for k-mer operations - micro-functions callable from C++
// Tests that Rust k-mer canonicalization matches C++ CKmer exactly

use crate::kmer::{Kmer, KmerMode};
use crate::kmer_extract::{remove_non_singletons, remove_non_singletons_with_duplicates, find_new_splitters_kmers};
use std::slice;

/// Extract all canonical k-mer values from a contig
///
/// Matches C++ CKmer behavior:
/// - Scans through contig building rolling k-mers
/// - Resets on non-ACGT bases (> 3)
/// - Returns canonical representation of each k-mer
///
/// # Safety
/// - Caller must ensure contig_data points to valid memory of contig_len bytes
/// - Returned array must be freed with ragc_free_kmer_array()
#[repr(C)]
pub struct KmerArray {
    /// Array of k-mer values (canonical)
    pub data: *mut u64,
    /// Number of k-mers
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn ragc_extract_canonical_kmers(
    contig_data: *const u8,
    contig_len: usize,
    k: u32,
) -> KmerArray {
    unsafe {
        let contig = slice::from_raw_parts(contig_data, contig_len);
        let mut kmers = Vec::new();
        let mut kmer = Kmer::new(k, KmerMode::Canonical);

        for &base in contig {
            if base > 3 {
                // Non-ACGT base, reset k-mer
                kmer.reset();
            } else {
                kmer.insert(base as u64);

                if kmer.is_full() {
                    kmers.push(kmer.data());
                }
            }
        }

        let mut result = kmers;
        let result_ptr = result.as_mut_ptr();
        let result_len = result.len();

        // Prevent Rust from freeing the allocation
        std::mem::forget(result);

        KmerArray {
            data: result_ptr,
            len: result_len,
        }
    }
}

/// Free a k-mer array allocated by ragc_extract_canonical_kmers()
///
/// # Safety
/// - Must only be called once per KmerArray
/// - array.data must be a valid pointer from ragc_extract_canonical_kmers()
#[no_mangle]
pub extern "C" fn ragc_free_kmer_array(array: KmerArray) {
    unsafe {
        if !array.data.is_null() && array.len > 0 {
            // Reconstruct the Vec and let it drop
            let _ = Vec::from_raw_parts(array.data, array.len, array.len);
        }
    }
}

/// Extract a single k-mer at a specific position
///
/// Returns the canonical k-mer value at position `pos` in the contig,
/// or u64::MAX if the k-mer cannot be extracted (position out of bounds,
/// non-ACGT base in k-mer window).
///
/// # Safety
/// - Caller must ensure contig_data points to valid memory of contig_len bytes
#[no_mangle]
pub extern "C" fn ragc_extract_kmer_at_position(
    contig_data: *const u8,
    contig_len: usize,
    k: u32,
    pos: usize,
) -> u64 {
    unsafe {
        let contig = slice::from_raw_parts(contig_data, contig_len);
        let k = k as usize;

        // Check bounds
        if pos + k > contig_len {
            return u64::MAX;
        }

        // Check for non-ACGT bases in the k-mer window
        for i in 0..k {
            if contig[pos + i] > 3 {
                return u64::MAX;
            }
        }

        // Build k-mer
        let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
        for i in 0..k {
            kmer.insert(contig[pos + i] as u64);
        }

        kmer.data()
    }
}

/// Remove non-singleton k-mers from a sorted vector
///
/// Modifies the vector in place to keep only k-mers that appear exactly once.
/// K-mers before `virtual_begin` are not checked and remain in the output.
///
/// Returns the new length of the vector.
///
/// # Safety
/// - vec_ptr must point to a valid vector allocation of at least vec_capacity elements
/// - The vector must be sorted (for correct singleton detection)
/// - The caller must resize the vector to the returned length
/// - Caller maintains ownership of the allocation
#[no_mangle]
pub extern "C" fn ragc_remove_non_singletons(
    vec_ptr: *mut u64,
    vec_len: usize,
    vec_capacity: usize,
    virtual_begin: usize,
) -> usize {
    unsafe {
        // Reconstruct the Vec from raw parts (temporarily borrow ownership)
        let mut vec = Vec::from_raw_parts(vec_ptr, vec_len, vec_capacity);

        // Call the Rust implementation
        remove_non_singletons(&mut vec, virtual_begin);

        // Get the new length after modification
        let new_len = vec.len();

        // Prevent Rust from freeing the allocation (C++ owns it)
        std::mem::forget(vec);

        new_len
    }
}

/// Result of remove_non_singletons_with_duplicates containing both new lengths
#[repr(C)]
pub struct RemoveSingletonsResult {
    /// New length of the main vector (singletons only)
    pub vec_new_len: usize,
    /// New length of the duplicates vector
    pub dup_new_len: usize,
}

/// Remove non-singleton k-mers from a sorted vector, collecting duplicates
///
/// Modifies the vector in place to keep only k-mers that appear exactly once,
/// and collects k-mers appearing more than once into a separate vector.
/// K-mers before `virtual_begin` are not checked and remain in the output.
///
/// Returns a struct containing the new lengths of both vectors.
///
/// # Safety
/// - vec_ptr must point to a valid vector allocation of at least vec_capacity elements
/// - dup_ptr must point to a valid vector allocation
/// - The main vector must be sorted (for correct singleton detection)
/// - The caller must resize both vectors to the returned lengths
/// - Caller maintains ownership of both allocations
/// - The duplicates vector will be cleared and filled with new values
#[no_mangle]
pub extern "C" fn ragc_remove_non_singletons_with_duplicates(
    vec_ptr: *mut u64,
    vec_len: usize,
    vec_capacity: usize,
    dup_ptr: *mut u64,
    dup_len: usize,
    dup_capacity: usize,
    virtual_begin: usize,
) -> RemoveSingletonsResult {
    unsafe {
        // Reconstruct both Vecs from raw parts (temporarily borrow ownership)
        let mut vec = Vec::from_raw_parts(vec_ptr, vec_len, vec_capacity);
        let mut duplicated = Vec::from_raw_parts(dup_ptr, dup_len, dup_capacity);

        // Call the Rust implementation
        remove_non_singletons_with_duplicates(&mut vec, &mut duplicated, virtual_begin);

        // Get the new lengths after modification
        let new_vec_len = vec.len();
        let new_dup_len = duplicated.len();

        // Prevent Rust from freeing the allocations (C++ owns them)
        std::mem::forget(vec);
        std::mem::forget(duplicated);

        RemoveSingletonsResult {
            vec_new_len: new_vec_len,
            dup_new_len: new_dup_len,
        }
    }
}

/// Find new splitter k-mers from a contig by excluding reference k-mers
///
/// Implements the k-mer filtering workflow from C++ AGC's find_new_splitters():
/// 1. Extract canonical k-mers from contig
/// 2. Filter to singletons only
/// 3. Exclude k-mers that appear in reference singletons
/// 4. Exclude k-mers that appear in reference duplicates
///
/// # Parameters
/// - contig_data: Pointer to contig sequence data (numeric encoding: A=0, C=1, G=2, T=3)
/// - contig_len: Length of contig
/// - k: K-mer length
/// - candidate_kmers_ptr: Pointer to sorted reference singleton k-mers
/// - candidate_kmers_len: Length of reference singleton array
/// - candidate_kmers_offset: Offset to start reading from candidate k-mers
/// - duplicated_kmers_ptr: Pointer to sorted reference duplicate k-mers
/// - duplicated_kmers_len: Length of reference duplicate array
///
/// # Returns
/// KmerArray containing novel k-mer values (must be freed with ragc_free_kmer_array)
///
/// # Safety
/// - Caller must ensure all pointers point to valid memory
/// - All k-mer arrays must be sorted
/// - Returned array must be freed with ragc_free_kmer_array()
#[no_mangle]
pub extern "C" fn ragc_find_new_splitters_kmers(
    contig_data: *const u8,
    contig_len: usize,
    k: u32,
    candidate_kmers_ptr: *const u64,
    candidate_kmers_len: usize,
    candidate_kmers_offset: usize,
    duplicated_kmers_ptr: *const u64,
    duplicated_kmers_len: usize,
) -> KmerArray {
    unsafe {
        let contig = slice::from_raw_parts(contig_data, contig_len);
        let candidate_kmers = slice::from_raw_parts(candidate_kmers_ptr, candidate_kmers_len);
        let duplicated_kmers = slice::from_raw_parts(duplicated_kmers_ptr, duplicated_kmers_len);

        let mut result = find_new_splitters_kmers(
            contig,
            k,
            candidate_kmers,
            candidate_kmers_offset,
            duplicated_kmers,
        );

        let len = result.len();
        let ptr = result.as_mut_ptr();
        std::mem::forget(result);

        KmerArray { data: ptr, len }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_canonical_kmers() {
        // ACGT sequence
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let k = 3;
        let array = ragc_extract_canonical_kmers(contig.as_ptr(), contig.len(), k);

        unsafe {
            let kmers = slice::from_raw_parts(array.data, array.len);
            // Should have 6 k-mers (length - k + 1)
            assert_eq!(kmers.len(), 6);
        }

        ragc_free_kmer_array(array);
    }

    #[test]
    fn test_extract_kmers_with_reset() {
        // Sequence with non-ACGT base (N = 4)
        let contig = vec![0, 1, 2, 4, 0, 1, 2, 3]; // ACGTNACGT
        let k = 3;
        let array = ragc_extract_canonical_kmers(contig.as_ptr(), contig.len(), k);

        unsafe {
            let kmers = slice::from_raw_parts(array.data, array.len);
            // Should have 3 k-mers: ACG (0-2), then reset at N, then ACG (4-6), CGT (5-7)
            assert_eq!(kmers.len(), 3);
        }

        ragc_free_kmer_array(array);
    }

    #[test]
    fn test_extract_kmer_at_position() {
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let k = 3;

        // Extract at position 0: ACG
        let kmer0 = ragc_extract_kmer_at_position(contig.as_ptr(), contig.len(), k, 0);
        assert_ne!(kmer0, u64::MAX);

        // Extract at position 5: CGT
        let kmer5 = ragc_extract_kmer_at_position(contig.as_ptr(), contig.len(), k, 5);
        assert_ne!(kmer5, u64::MAX);

        // Out of bounds
        let kmer_oob = ragc_extract_kmer_at_position(contig.as_ptr(), contig.len(), k, 10);
        assert_eq!(kmer_oob, u64::MAX);
    }

    #[test]
    fn test_extract_kmer_with_n() {
        let contig = vec![0, 4, 2, 3]; // ANGT
        let k = 3;

        // Position 0 contains N, should return MAX
        let kmer = ragc_extract_kmer_at_position(contig.as_ptr(), contig.len(), k, 0);
        assert_eq!(kmer, u64::MAX);
    }
}
