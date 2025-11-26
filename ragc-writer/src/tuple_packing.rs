// Tuple Packing for AGC Segment Data
// Matches C++ AGC implementation in segment.h

/// Pack bytes into tuples based on maximum value
/// This matches C++ AGC's bytes2tuples implementation
pub fn bytes_to_tuples(bytes: &[u8]) -> Vec<u8> {
    if bytes.is_empty() {
        return vec![0x10]; // Empty with marker
    }

    let max_elem = *bytes.iter().max().unwrap();

    if max_elem < 4 {
        pack_tuples::<4, 4>(bytes)
    } else if max_elem < 6 {
        pack_tuples::<3, 6>(bytes)
    } else if max_elem < 16 {
        pack_tuples::<2, 16>(bytes)
    } else {
        // No packing needed
        let mut result = bytes.to_vec();
        result.push(0x10); // Marker: no packing
        result
    }
}

/// Unpack tuples back to bytes
/// This matches C++ AGC's tuples2bytes implementation
pub fn tuples_to_bytes(tuples: &[u8]) -> Vec<u8> {
    if tuples.is_empty() {
        return Vec::new();
    }

    let marker = tuples[tuples.len() - 1];
    let no_bytes = marker >> 4;
    let trailing_bytes = marker & 0xf;

    if no_bytes == 1 {
        // No packing was used, just return data without marker
        return tuples[..tuples.len() - 1].to_vec();
    }

    // Output size calculation matches C++ AGC line 99
    let output_size = (tuples.len() - 2) * (no_bytes as usize) + (trailing_bytes as usize);
    let mut result = vec![0u8; output_size];

    match no_bytes {
        2 => unpack_tuples::<2, 16>(&tuples[..tuples.len() - 1], &mut result, output_size),
        3 => unpack_tuples::<3, 6>(&tuples[..tuples.len() - 1], &mut result, output_size),
        4 => unpack_tuples::<4, 4>(&tuples[..tuples.len() - 1], &mut result, output_size),
        _ => panic!("Invalid no_bytes: {no_bytes}"),
    }

    result
}

/// Pack N values per byte, where each value is in range [0, MAX)
/// Matches C++ bytes2tuples_impl
fn pack_tuples<const N: usize, const MAX: u8>(bytes: &[u8]) -> Vec<u8> {
    let mut result = Vec::new();
    let mut i = 0;

    // Pack full tuples (C++ lines 122-130)
    while i + N <= bytes.len() {
        let mut c: u32 = 0;
        for j in 0..N {
            c = c * (MAX as u32) + (bytes[i + j] as u32);
        }
        result.push(c as u8);
        i += N;
    }

    // Pack trailing bytes (C++ lines 132-135)
    // ALWAYS add trailing tuple (even if c=0)
    let mut c: u32 = 0;
    while i < bytes.len() {
        c = c * (MAX as u32) + (bytes[i] as u32);
        i += 1;
    }
    result.push(c as u8);

    // Add marker byte (C++ line 137): (NO_BYTES << 4) + (v_bytes.size() % NO_BYTES)
    let marker = ((N as u8) << 4) | ((bytes.len() % N) as u8);
    result.push(marker);

    result
}

/// Unpack N values per byte
/// Matches C++ tuples2bytes_impl
fn unpack_tuples<const N: usize, const MAX: u8>(
    tuples: &[u8],
    output: &mut [u8],
    output_size: usize,
) {
    let mut i = 0; // tuple index
    let mut j = 0; // output index

    // Unpack full tuples (C++ lines 148-157)
    while j + N <= output_size {
        let mut c = tuples[i] as u32;

        // Extract N values in reverse order (C++ lines 152-156)
        for k in (0..N).rev() {
            output[j + k] = (c % (MAX as u32)) as u8;
            c /= MAX as u32;
        }

        i += 1;
        j += N;
    }

    // Handle trailing bytes (C++ lines 159-168)
    let n = output_size % N;
    if n > 0 {
        let mut c = tuples[i] as u32;

        for k in (0..n).rev() {
            output[j + k] = (c % (MAX as u32)) as u8;
            c /= MAX as u32;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_dna() {
        // DNA: values 0-3 (A=0, C=1, G=2, T=3)
        let dna = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1];
        let packed = bytes_to_tuples(&dna);
        let unpacked = tuples_to_bytes(&packed);
        assert_eq!(dna, unpacked);
    }

    #[test]
    fn test_pack_unpack_with_n() {
        // DNA with N (N=4)
        let data = vec![0, 1, 2, 3, 4, 1, 2, 3];
        let packed = bytes_to_tuples(&data);
        let unpacked = tuples_to_bytes(&packed);
        assert_eq!(data, unpacked);
    }

    #[test]
    fn test_pack_unpack_threshold_6() {
        // Values 0-5 (should use 3-per-byte packing)
        let data = vec![0, 1, 2, 3, 4, 5, 0, 1, 2];
        let packed = bytes_to_tuples(&data);
        let unpacked = tuples_to_bytes(&packed);
        assert_eq!(data, unpacked);
    }

    #[test]
    fn test_pack_unpack_threshold_7() {
        // Values 0-7 (should use 2-per-byte packing, NOT 3-per-byte!)
        let data = vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 1];
        let packed = bytes_to_tuples(&data);
        let unpacked = tuples_to_bytes(&packed);
        assert_eq!(data, unpacked);
    }

    #[test]
    fn test_empty() {
        let data = vec![];
        let packed = bytes_to_tuples(&data);
        let unpacked = tuples_to_bytes(&packed);
        assert_eq!(data, unpacked);
    }

    #[test]
    fn test_no_packing_needed() {
        // Values >= 16 should not be packed
        let data = vec![16, 20, 100, 200];
        let packed = bytes_to_tuples(&data);
        // Should be: data + 0x10 marker
        assert_eq!(packed.len(), data.len() + 1);
        assert_eq!(packed[packed.len() - 1], 0x10);

        let unpacked = tuples_to_bytes(&packed);
        assert_eq!(data, unpacked);
    }
}
