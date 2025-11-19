// FFI helpers for segmentation - micro-functions callable from C++
// Ship of Theseus approach: replace tiny pieces one at a time

use std::os::raw::c_char;
use std::slice;

/// C representation of a contig slice result
#[repr(C)]
pub struct ContigSlice {
    /// Pointer to the slice data
    pub data: *mut u8,
    /// Length of the slice
    pub len: usize,
}

/// Extract a slice of a contig
///
/// Matches C++ AGC's get_part() function (agc_compressor.cpp:2101-2107):
/// ```cpp
/// contig_t CAGCCompressor::get_part(const contig_t& contig, uint64_t pos, uint64_t len)
/// {
///     if (pos + len < contig.size())
///         return contig_t(contig.begin() + pos, contig.begin() + pos + len);
///     else
///         return contig_t(contig.begin() + pos, contig.end());
/// }
/// ```
///
/// # Safety
/// - Caller must ensure contig_data points to valid memory of contig_len bytes
/// - Returned slice is owned by caller and must be freed with ragc_free_contig_slice()
#[no_mangle]
pub extern "C" fn ragc_get_contig_part(
    contig_data: *const u8,
    contig_len: usize,
    pos: u64,
    len: u64,
) -> ContigSlice {
    unsafe {
        let contig = slice::from_raw_parts(contig_data, contig_len);

        let pos = pos as usize;
        let len = len as usize;

        // Match C++ logic exactly
        let result_slice = if pos + len < contig.len() {
            // Can extract full requested length
            &contig[pos..pos + len]
        } else {
            // Would exceed contig end, return until end
            &contig[pos..]
        };

        // Allocate owned copy for C++
        let mut result = result_slice.to_vec();
        let result_ptr = result.as_mut_ptr();
        let result_len = result.len();

        // Prevent Rust from freeing the allocation
        std::mem::forget(result);

        ContigSlice {
            data: result_ptr,
            len: result_len,
        }
    }
}

/// Free a contig slice allocated by ragc_get_contig_part()
///
/// # Safety
/// - Must only be called once per ContigSlice returned from ragc_get_contig_part()
/// - slice.data must be a valid pointer from ragc_get_contig_part()
#[no_mangle]
pub extern "C" fn ragc_free_contig_slice(slice: ContigSlice) {
    unsafe {
        if !slice.data.is_null() && slice.len > 0 {
            // Reconstruct the Vec and let it drop
            let _ = Vec::from_raw_parts(slice.data, slice.len, slice.len);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_part_full_length() {
        let contig = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let slice = ragc_get_contig_part(contig.as_ptr(), contig.len(), 2, 5);

        unsafe {
            let result = slice::from_raw_parts(slice.data, slice.len);
            assert_eq!(result, &[2, 3, 4, 5, 6]);
        }

        ragc_free_contig_slice(slice);
    }

    #[test]
    fn test_get_part_exceeds_length() {
        let contig = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let slice = ragc_get_contig_part(contig.as_ptr(), contig.len(), 7, 10);

        unsafe {
            let result = slice::from_raw_parts(slice.data, slice.len);
            assert_eq!(result, &[7, 8, 9]);
        }

        ragc_free_contig_slice(slice);
    }

    #[test]
    fn test_get_part_exact_end() {
        let contig = vec![0, 1, 2, 3, 4];
        let slice = ragc_get_contig_part(contig.as_ptr(), contig.len(), 2, 3);

        unsafe {
            let result = slice::from_raw_parts(slice.data, slice.len);
            assert_eq!(result, &[2, 3, 4]);
        }

        ragc_free_contig_slice(slice);
    }
}
