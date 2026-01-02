// FFI helper for DNA base validation
// Simple but fundamental check used throughout segmentation

/// Check if a base is a valid ACGT nucleotide
///
/// Matches C++ AGC's base validation (agc_compressor.cpp:2025):
/// ```cpp
/// if (x >> 2)         // x > 3
///     kmer.Reset();
/// ```
///
/// Valid bases are encoded as:
/// - A = 0
/// - C = 1
/// - G = 2
/// - T = 3
/// - N or other = 4+ (invalid, triggers k-mer reset)
///
/// # Arguments
/// * `base` - Encoded base value
///
/// # Returns
/// true if base is valid ACGT (0-3), false otherwise
#[no_mangle]
pub extern "C" fn ragc_is_valid_base(base: u8) -> bool {
    // Match C++ AGC logic exactly
    // C++ uses: if (x >> 2) to check if x > 3
    // We can use: base <= 3 (equivalent and clearer)
    base <= 3
}

/// Check if a base requires k-mer reset
///
/// This is the inverse of is_valid_base() - returns true if the base
/// should trigger a k-mer reset (non-ACGT).
///
/// Matches C++ AGC: if (x >> 2) kmer.Reset();
#[no_mangle]
pub extern "C" fn ragc_should_reset_kmer(base: u8) -> bool {
    base > 3
}

/// Validate an entire sequence, counting valid and invalid bases
///
/// Returns (n_valid, n_invalid) counts
#[repr(C)]
pub struct BaseCounts {
    pub n_valid: usize,
    pub n_invalid: usize,
}

#[no_mangle]
pub extern "C" fn ragc_count_base_validity(sequence: *const u8, length: usize) -> BaseCounts {
    unsafe {
        let seq = std::slice::from_raw_parts(sequence, length);

        let mut n_valid = 0;
        let mut n_invalid = 0;

        for &base in seq {
            if ragc_is_valid_base(base) {
                n_valid += 1;
            } else {
                n_invalid += 1;
            }
        }

        BaseCounts { n_valid, n_invalid }
    }
}

/// Find positions of invalid bases (N or other) in a sequence
///
/// Returns array of positions where invalid bases occur.
/// Useful for debugging or understanding where k-mer resets happen.
///
/// # Safety
/// - sequence must point to valid memory of length bytes
/// - Returned array must be freed with ragc_free_position_array()
#[repr(C)]
pub struct PositionArray {
    pub data: *mut usize,
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn ragc_find_invalid_base_positions(
    sequence: *const u8,
    length: usize,
) -> PositionArray {
    unsafe {
        let seq = std::slice::from_raw_parts(sequence, length);

        let mut positions: Vec<usize> = seq
            .iter()
            .enumerate()
            .filter(|(_, &base)| !ragc_is_valid_base(base))
            .map(|(pos, _)| pos)
            .collect();

        let ptr = positions.as_mut_ptr();
        let len = positions.len();

        std::mem::forget(positions);

        PositionArray { data: ptr, len }
    }
}

#[no_mangle]
pub extern "C" fn ragc_free_position_array(array: PositionArray) {
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
    fn test_valid_bases() {
        // Valid ACGT bases
        assert!(ragc_is_valid_base(0)); // A
        assert!(ragc_is_valid_base(1)); // C
        assert!(ragc_is_valid_base(2)); // G
        assert!(ragc_is_valid_base(3)); // T

        // Invalid bases
        assert!(!ragc_is_valid_base(4)); // N
        assert!(!ragc_is_valid_base(5));
        assert!(!ragc_is_valid_base(255));
    }

    #[test]
    fn test_should_reset() {
        // Should NOT reset for valid bases
        assert!(!ragc_should_reset_kmer(0));
        assert!(!ragc_should_reset_kmer(1));
        assert!(!ragc_should_reset_kmer(2));
        assert!(!ragc_should_reset_kmer(3));

        // SHOULD reset for invalid bases
        assert!(ragc_should_reset_kmer(4));
        assert!(ragc_should_reset_kmer(5));
        assert!(ragc_should_reset_kmer(255));
    }

    #[test]
    fn test_count_validity() {
        let sequence = vec![0, 1, 2, 3, 4, 0, 1, 5]; // ACGTNACX
        let counts = ragc_count_base_validity(sequence.as_ptr(), sequence.len());

        assert_eq!(counts.n_valid, 6); // A,C,G,T,A,C
        assert_eq!(counts.n_invalid, 2); // N,X
    }

    #[test]
    fn test_find_invalid_positions() {
        let sequence = vec![0, 1, 4, 2, 3, 5, 0]; // ACNGTXA
        let positions = ragc_find_invalid_base_positions(sequence.as_ptr(), sequence.len());

        unsafe {
            let pos_slice = std::slice::from_raw_parts(positions.data, positions.len);
            assert_eq!(pos_slice, &[2, 5]); // Positions of N and X
        }

        ragc_free_position_array(positions);
    }

    #[test]
    fn test_all_valid() {
        let sequence = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let positions = ragc_find_invalid_base_positions(sequence.as_ptr(), sequence.len());

        assert_eq!(positions.len, 0); // No invalid bases

        ragc_free_position_array(positions);
    }
}
