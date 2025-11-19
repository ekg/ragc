// FFI helper for reverse complement operations
// Core DNA sequence transformation used in segmentation

/// Complement a single base
///
/// A (0) ↔ T (3)
/// C (1) ↔ G (2)
/// N (4+) → N (unchanged)
#[inline]
fn complement_base(base: u8) -> u8 {
    if base < 4 {
        3 - base
    } else {
        base  // N or invalid stays unchanged
    }
}

/// Reverse complement a sequence in-place
///
/// Matches C++ AGC's CAGCBasic::reverse_complement() (agc_basic.cpp:257)
///
/// # Safety
/// - sequence must point to valid mutable memory of length bytes
#[no_mangle]
pub extern "C" fn ragc_reverse_complement_inplace(
    sequence: *mut u8,
    length: usize,
) {
    if length == 0 {
        return;
    }

    unsafe {
        let seq = std::slice::from_raw_parts_mut(sequence, length);

        let mut i = 0;
        let mut j = length - 1;

        while i < j {
            let x = complement_base(seq[j]);
            let y = complement_base(seq[i]);

            seq[i] = x;
            seq[j] = y;

            i += 1;
            j -= 1;
        }

        // Handle middle element if odd length
        if i == j {
            seq[i] = complement_base(seq[i]);
        }
    }
}

/// Reverse complement copy - create reverse complement in new buffer
///
/// Matches C++ AGC's CAGCBasic::reverse_complement_copy() (agc_basic.cpp:282)
///
/// # Safety
/// - src must point to valid memory of src_len bytes
/// - Returned buffer must be freed with ragc_free_sequence()
#[repr(C)]
pub struct Sequence {
    pub data: *mut u8,
    pub len: usize,
}

#[no_mangle]
pub extern "C" fn ragc_reverse_complement_copy(
    src: *const u8,
    src_len: usize,
) -> Sequence {
    if src_len == 0 {
        return Sequence {
            data: std::ptr::null_mut(),
            len: 0,
        };
    }

    unsafe {
        let src_slice = std::slice::from_raw_parts(src, src_len);
        let mut dest: Vec<u8> = src_slice
            .iter()
            .rev()
            .map(|&base| complement_base(base))
            .collect();

        let ptr = dest.as_mut_ptr();
        let len = dest.len();

        std::mem::forget(dest);

        Sequence { data: ptr, len }
    }
}

#[no_mangle]
pub extern "C" fn ragc_free_sequence(seq: Sequence) {
    unsafe {
        if !seq.data.is_null() && seq.len > 0 {
            let _ = Vec::from_raw_parts(seq.data, seq.len, seq.len);
        }
    }
}

/// Complement a single base (public wrapper)
#[no_mangle]
pub extern "C" fn ragc_complement_base(base: u8) -> u8 {
    complement_base(base)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_base() {
        assert_eq!(complement_base(0), 3); // A -> T
        assert_eq!(complement_base(1), 2); // C -> G
        assert_eq!(complement_base(2), 1); // G -> C
        assert_eq!(complement_base(3), 0); // T -> A
        assert_eq!(complement_base(4), 4); // N -> N
        assert_eq!(complement_base(5), 5); // Invalid -> Invalid
    }

    #[test]
    fn test_reverse_complement_inplace() {
        // Test with non-palindrome: ACG -> CGT
        let mut seq = vec![0, 1, 2]; // ACG (0,1,2)
        ragc_reverse_complement_inplace(seq.as_mut_ptr(), seq.len());
        assert_eq!(seq, vec![1, 2, 3]); // CGT (1,2,3)

        let mut seq2 = vec![0, 0, 0]; // AAA
        ragc_reverse_complement_inplace(seq2.as_mut_ptr(), seq2.len());
        assert_eq!(seq2, vec![3, 3, 3]); // AAA -> TTT
    }

    #[test]
    fn test_reverse_complement_copy() {
        // ACG = 0,1,2
        // reverse = 2,1,0 = GCA
        // complement of GCA = 1,2,3 = CGT
        let src = vec![0, 1, 2]; // ACG
        let result = ragc_reverse_complement_copy(src.as_ptr(), src.len());

        unsafe {
            let dest = std::slice::from_raw_parts(result.data, result.len);
            assert_eq!(dest, &[1, 2, 3]); // CGT
        }

        ragc_free_sequence(result);
    }

    #[test]
    fn test_reverse_complement_with_n() {
        // ACNG = 0,1,4,2
        // reverse = 2,4,1,0 = GNCA
        // complement of GNCA = 1,4,2,3 = CNGT
        let src = vec![0, 1, 4, 2]; // ACNG
        let result = ragc_reverse_complement_copy(src.as_ptr(), src.len());

        unsafe {
            let dest = std::slice::from_raw_parts(result.data, result.len);
            assert_eq!(dest, &[1, 4, 2, 3]); // CNGT
        }

        ragc_free_sequence(result);
    }

    #[test]
    fn test_empty_sequence() {
        let empty: Vec<u8> = vec![];
        let result = ragc_reverse_complement_copy(empty.as_ptr(), 0);

        assert_eq!(result.len, 0);
        assert!(result.data.is_null());
    }
}
