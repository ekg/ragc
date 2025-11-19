// FFI wrapper for preprocessing operations

use crate::preprocessing::preprocess_raw_contig;

/// Preprocess raw contig (in-place) by converting ASCII to numeric codes
///
/// Matches C++ AGC's preprocess_raw_contig() behavior:
/// - Filters bytes >= 64 (letter characters)
/// - Converts using cnv_num lookup table
/// - Removes non-letter characters
/// - Returns new length after preprocessing
///
/// # Safety
/// - contig_ptr must point to a valid vector allocation of at least contig_capacity elements
/// - The caller must resize the vector to the returned length
/// - Caller maintains ownership of the allocation
#[no_mangle]
pub extern "C" fn ragc_preprocess_raw_contig(
    contig_ptr: *mut u8,
    contig_len: usize,
    contig_capacity: usize,
) -> usize {
    unsafe {
        // Reconstruct the Vec from raw parts (temporarily borrow ownership)
        let mut contig = Vec::from_raw_parts(contig_ptr, contig_len, contig_capacity);

        // Call the Rust implementation
        preprocess_raw_contig(&mut contig);

        // Get the new length after preprocessing
        let new_len = contig.len();

        // Prevent Rust from freeing the allocation (C++ owns it)
        std::mem::forget(contig);

        new_len
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ffi_preprocess_raw_contig() {
        let mut data = b"ACGT".to_vec();
        let ptr = data.as_mut_ptr();
        let len = data.len();
        let cap = data.capacity();

        // Prevent Rust from freeing - we're simulating C++ ownership
        std::mem::forget(data);

        // Call FFI function
        let new_len = ragc_preprocess_raw_contig(ptr, len, cap);

        // Reconstruct to check result
        let result = unsafe { Vec::from_raw_parts(ptr, new_len, cap) };
        assert_eq!(result, vec![0, 1, 2, 3]);
    }
}
