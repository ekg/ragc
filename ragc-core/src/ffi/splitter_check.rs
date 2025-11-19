// FFI helper for splitter checking - micro-function callable from C++
// Tests that our splitter data structures work correctly

use std::collections::HashSet;

/// Check if a k-mer is a splitter
///
/// Matches C++ AGC's splitter checking logic (agc_compressor.cpp:2034):
/// ```cpp
/// if (bloom_splitters.check(d) && hs_splitters.check(d))
/// ```
///
/// This is the critical decision point in compress_contig() that determines
/// where to split segments.
///
/// # Safety
/// - splitters_ptr must point to a valid array of splitters_len elements
/// - This function does not take ownership of the splitters array
#[no_mangle]
pub extern "C" fn ragc_is_splitter(
    kmer_value: u64,
    splitters_ptr: *const u64,
    splitters_len: usize,
) -> bool {
    unsafe {
        let splitters = std::slice::from_raw_parts(splitters_ptr, splitters_len);

        // Simple linear search for now (could optimize with HashSet if needed)
        // But C++ AGC uses hash set + bloom filter, so we should match that
        splitters.binary_search(&kmer_value).is_ok()
    }
}

/// Create a splitter checker with proper data structures
///
/// This matches C++ AGC's use of both bloom filter and hash set.
/// For now, we'll use a simpler approach and just return a sorted vector
/// that can be binary searched.
///
/// Returns a pointer to a sorted array of splitters that must be freed
/// with ragc_free_splitter_checker().
#[repr(C)]
pub struct SplitterChecker {
    pub splitters: *mut u64,
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn ragc_create_splitter_checker(
    splitters_ptr: *const u64,
    splitters_len: usize,
) -> SplitterChecker {
    unsafe {
        let splitters_slice = std::slice::from_raw_parts(splitters_ptr, splitters_len);
        let mut splitters_vec: Vec<u64> = splitters_slice.to_vec();

        // Sort for binary search
        splitters_vec.sort_unstable();

        let ptr = splitters_vec.as_mut_ptr();
        let len = splitters_vec.len();

        std::mem::forget(splitters_vec);

        SplitterChecker {
            splitters: ptr,
            len,
        }
    }
}

#[no_mangle]
pub extern "C" fn ragc_free_splitter_checker(checker: SplitterChecker) {
    unsafe {
        if !checker.splitters.is_null() && checker.len > 0 {
            let _ = Vec::from_raw_parts(checker.splitters, checker.len, checker.len);
        }
    }
}

/// Batch check if multiple k-mers are splitters
///
/// More efficient than calling ragc_is_splitter() repeatedly.
/// Returns array of bools indicating which k-mers are splitters.
///
/// # Safety
/// - kmers_ptr must point to valid array of kmers_len elements
/// - splitters_ptr must point to valid array of splitters_len elements
/// - Returned array must be freed with ragc_free_bool_array()
#[repr(C)]
pub struct BoolArray {
    pub data: *mut bool,
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn ragc_check_splitters_batch(
    kmers_ptr: *const u64,
    kmers_len: usize,
    splitters_ptr: *const u64,
    splitters_len: usize,
) -> BoolArray {
    unsafe {
        let kmers = std::slice::from_raw_parts(kmers_ptr, kmers_len);
        let splitters_slice = std::slice::from_raw_parts(splitters_ptr, splitters_len);

        // Build HashSet for O(1) lookups
        let splitter_set: HashSet<u64> = splitters_slice.iter().copied().collect();

        // Check each k-mer
        let mut results: Vec<bool> = kmers
            .iter()
            .map(|kmer| splitter_set.contains(kmer))
            .collect();

        let ptr = results.as_mut_ptr();
        let len = results.len();

        std::mem::forget(results);

        BoolArray { data: ptr, len }
    }
}

#[no_mangle]
pub extern "C" fn ragc_free_bool_array(array: BoolArray) {
    unsafe {
        if !array.data.is_null() && array.len > 0 {
            let _ = Vec::from_raw_parts(array.data, array.len, array.len);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_splitter() {
        let splitters = vec![100u64, 200, 300, 400, 500];

        // K-mer is a splitter
        assert!(ragc_is_splitter(300, splitters.as_ptr(), splitters.len()));

        // K-mer is not a splitter
        assert!(!ragc_is_splitter(350, splitters.as_ptr(), splitters.len()));

        // First and last
        assert!(ragc_is_splitter(100, splitters.as_ptr(), splitters.len()));
        assert!(ragc_is_splitter(500, splitters.as_ptr(), splitters.len()));
    }

    #[test]
    fn test_splitter_checker() {
        let splitters = vec![300u64, 100, 500, 200, 400]; // Unsorted
        let checker = ragc_create_splitter_checker(splitters.as_ptr(), splitters.len());

        unsafe {
            let sorted = std::slice::from_raw_parts(checker.splitters, checker.len);

            // Should be sorted
            assert_eq!(sorted, &[100, 200, 300, 400, 500]);
        }

        ragc_free_splitter_checker(checker);
    }

    #[test]
    fn test_batch_check() {
        let splitters = vec![100u64, 200, 300];
        let kmers = vec![50u64, 100, 150, 200, 250, 300, 350];

        let results = ragc_check_splitters_batch(
            kmers.as_ptr(),
            kmers.len(),
            splitters.as_ptr(),
            splitters.len(),
        );

        unsafe {
            let checks = std::slice::from_raw_parts(results.data, results.len);
            assert_eq!(checks, &[false, true, false, true, false, true, false]);
        }

        ragc_free_bool_array(results);
    }
}
