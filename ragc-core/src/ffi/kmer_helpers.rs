// FFI helpers for k-mer operations - micro-functions callable from C++
// Tests that Rust k-mer canonicalization matches C++ CKmer exactly

use crate::kmer::{Kmer, KmerMode};
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
            // Should have 4 k-mers: ACG, then reset at N, then ACG, CGT
            assert_eq!(kmers.len(), 4);
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
