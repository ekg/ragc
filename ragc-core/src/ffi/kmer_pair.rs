// FFI helper for k-mer pair operations
// Used for grouping segments by their terminal k-mers

use std::cmp::{min, max};

/// Create an ordered k-mer pair (min, max)
///
/// Matches C++ AGC's minmax() function usage (agc_compressor.cpp:1419, 1427):
/// ```cpp
/// pk = minmax(split_match.first, kmer2.data());
/// pk = minmax(kmer1.data(), split_match.first);
/// ```
///
/// Returns pair with smaller value first, larger value second.
/// This ensures consistent grouping regardless of k-mer order.
///
/// # Arguments
/// * `kmer1` - First k-mer value
/// * `kmer2` - Second k-mer value
///
/// # Returns
/// Ordered pair (min, max)
#[repr(C)]
pub struct KmerPair {
    pub first: u64,
    pub second: u64,
}

#[no_mangle]
pub extern "C" fn ragc_create_kmer_pair(kmer1: u64, kmer2: u64) -> KmerPair {
    KmerPair {
        first: min(kmer1, kmer2),
        second: max(kmer1, kmer2),
    }
}

/// Check if two k-mer pairs are equal
///
/// Used for map lookups and comparisons in C++ AGC.
#[no_mangle]
pub extern "C" fn ragc_kmer_pair_equals(
    pair1_first: u64,
    pair1_second: u64,
    pair2_first: u64,
    pair2_second: u64,
) -> bool {
    pair1_first == pair2_first && pair1_second == pair2_second
}

/// Check if a k-mer value is the empty/invalid marker (~0ull)
///
/// Matches C++ AGC: if (pk.first == ~0ull || pk.second == ~0ull)
#[no_mangle]
pub extern "C" fn ragc_is_empty_kmer(kmer: u64) -> bool {
    kmer == u64::MAX  // ~0ull in C++ = u64::MAX in Rust
}

/// Check if a k-mer pair is valid (neither value is empty marker)
///
/// Matches C++ AGC: pk.first != ~0ull && pk.second != ~0ull
#[no_mangle]
pub extern "C" fn ragc_is_valid_kmer_pair(first: u64, second: u64) -> bool {
    first != u64::MAX && second != u64::MAX
}

/// Create empty k-mer pair (~0ull, ~0ull)
///
/// Matches C++ AGC's pk_empty constant
#[no_mangle]
pub extern "C" fn ragc_create_empty_kmer_pair() -> KmerPair {
    KmerPair {
        first: u64::MAX,
        second: u64::MAX,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_kmer_pair() {
        // Forward order
        let pair = ragc_create_kmer_pair(100, 200);
        assert_eq!(pair.first, 100);
        assert_eq!(pair.second, 200);

        // Reverse order - should swap
        let pair = ragc_create_kmer_pair(200, 100);
        assert_eq!(pair.first, 100);
        assert_eq!(pair.second, 200);

        // Equal values
        let pair = ragc_create_kmer_pair(150, 150);
        assert_eq!(pair.first, 150);
        assert_eq!(pair.second, 150);
    }

    #[test]
    fn test_kmer_pair_equals() {
        assert!(ragc_kmer_pair_equals(100, 200, 100, 200));
        assert!(!ragc_kmer_pair_equals(100, 200, 100, 201));
        assert!(!ragc_kmer_pair_equals(100, 200, 101, 200));
        assert!(!ragc_kmer_pair_equals(100, 200, 200, 100));
    }

    #[test]
    fn test_is_empty_kmer() {
        assert!(ragc_is_empty_kmer(u64::MAX));
        assert!(!ragc_is_empty_kmer(0));
        assert!(!ragc_is_empty_kmer(100));
        assert!(!ragc_is_empty_kmer(u64::MAX - 1));
    }

    #[test]
    fn test_is_valid_kmer_pair() {
        assert!(ragc_is_valid_kmer_pair(100, 200));
        assert!(ragc_is_valid_kmer_pair(0, 0));
        assert!(!ragc_is_valid_kmer_pair(u64::MAX, 200));
        assert!(!ragc_is_valid_kmer_pair(100, u64::MAX));
        assert!(!ragc_is_valid_kmer_pair(u64::MAX, u64::MAX));
    }

    #[test]
    fn test_create_empty_kmer_pair() {
        let pair = ragc_create_empty_kmer_pair();
        assert_eq!(pair.first, u64::MAX);
        assert_eq!(pair.second, u64::MAX);
        assert!(ragc_is_empty_kmer(pair.first));
        assert!(ragc_is_empty_kmer(pair.second));
        assert!(!ragc_is_valid_kmer_pair(pair.first, pair.second));
    }

    #[test]
    fn test_realistic_kmer_values() {
        // Realistic k-mer values (21-mers have max value 2^42-1)
        let kmer1: u64 = 0x3FFFFFFFFFF; // Max 21-mer
        let kmer2: u64 = 0x123456789AB;

        let pair = ragc_create_kmer_pair(kmer1, kmer2);
        assert_eq!(pair.first, min(kmer1, kmer2));
        assert_eq!(pair.second, max(kmer1, kmer2));
        assert!(ragc_is_valid_kmer_pair(pair.first, pair.second));
    }
}
