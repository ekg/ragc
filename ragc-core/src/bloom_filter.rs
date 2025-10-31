// Bloom filter for fast splitter k-mer membership testing
// Matches C++ AGC's bloom filter (agc_compressor.h)
//
// Purpose: Fast probabilistic membership test before expensive hash set lookup
// Pattern: bloom_filter.check(kmer) && hash_set.contains(kmer)

use std::sync::Mutex;

/// Simple Bloom filter for k-mer lookups
/// Matches C++ AGC interface: check(), insert(), filling_factor(), resize()
pub struct BloomFilter {
    /// Bit array for bloom filter
    bits: Mutex<Vec<u64>>,
    /// Number of hash functions (C++ AGC uses 3)
    num_hashes: usize,
    /// Number of bits in filter
    size_bits: usize,
    /// Number of items inserted
    num_items: usize,
}

impl BloomFilter {
    /// Create new Bloom filter with specified bit size
    ///
    /// C++ AGC initialization: bloom_splitters.resize(expected_splitters * 8)
    /// Uses ~8 bits per expected item for good false positive rate
    pub fn new(size_bits: usize) -> Self {
        let num_words = (size_bits + 63) / 64; // Round up to word boundary
        BloomFilter {
            bits: Mutex::new(vec![0u64; num_words]),
            num_hashes: 3,                // C++ AGC uses 3 hash functions
            size_bits: size_bits.max(64), // Minimum size
            num_items: 0,
        }
    }

    /// Check if k-mer might be in the set (probabilistic)
    /// Returns true if POSSIBLY present, false if DEFINITELY absent
    ///
    /// Matches C++ AGC: bloom_splitters.check(kmer)
    #[inline]
    pub fn check(&self, kmer: u64) -> bool {
        let bits = self.bits.lock().unwrap();

        for i in 0..self.num_hashes {
            let hash = self.hash(kmer, i);
            let bit_index = hash % self.size_bits;
            let word_index = bit_index / 64;
            let bit_offset = bit_index % 64;

            if (bits[word_index] & (1u64 << bit_offset)) == 0 {
                return false; // Definitely not in set
            }
        }

        true // Possibly in set
    }

    /// Insert k-mer into bloom filter
    ///
    /// Matches C++ AGC: bloom_splitters.insert(kmer)
    pub fn insert(&mut self, kmer: u64) {
        let mut bits = self.bits.lock().unwrap();

        for i in 0..self.num_hashes {
            let hash = self.hash(kmer, i);
            let bit_index = hash % self.size_bits;
            let word_index = bit_index / 64;
            let bit_offset = bit_index % 64;

            bits[word_index] |= 1u64 << bit_offset;
        }

        self.num_items += 1;
    }

    /// Get current filling factor (ratio of items to capacity)
    ///
    /// Matches C++ AGC: bloom_splitters.filling_factor()
    /// C++ AGC resizes when filling_factor() > 0.3
    pub fn filling_factor(&self) -> f64 {
        // Approximate: items / (bits / bits_per_item)
        // With 3 hash functions, optimal is ~8 bits per item
        let capacity = self.size_bits / 8;
        if capacity == 0 {
            return 1.0;
        }
        self.num_items as f64 / capacity as f64
    }

    /// Resize bloom filter to new bit size
    /// Used when filling factor exceeds threshold (C++ AGC: > 0.3)
    ///
    /// Matches C++ AGC: bloom_splitters.resize(new_size)
    /// Note: Must re-insert all items after resize (caller's responsibility)
    pub fn resize(&mut self, new_size_bits: usize) {
        let num_words = (new_size_bits + 63) / 64;
        let mut bits = self.bits.lock().unwrap();
        *bits = vec![0u64; num_words];
        self.size_bits = new_size_bits.max(64);
        self.num_items = 0; // Reset count (items must be re-inserted)
    }

    /// Clear all bits (used for testing or re-initialization)
    pub fn clear(&mut self) {
        let mut bits = self.bits.lock().unwrap();
        bits.fill(0);
        self.num_items = 0;
    }

    /// Get current size in bits
    pub fn size_bits(&self) -> usize {
        self.size_bits
    }

    /// Compute hash for k-mer with seed
    /// Simple hash mixing for multiple independent hash functions
    #[inline]
    fn hash(&self, kmer: u64, seed: usize) -> usize {
        // Mix kmer with seed using simple multiplicative hash
        let mut h = kmer.wrapping_mul(0x9e3779b97f4a7c15u64); // Golden ratio
        h = h.wrapping_add(seed as u64);
        h ^= h >> 32;
        h = h.wrapping_mul(0x9e3779b97f4a7c15u64);
        h ^= h >> 32;
        h as usize
    }
}

impl Default for BloomFilter {
    fn default() -> Self {
        // Default size: 64KB = 512K bits (~64K items at 8 bits/item)
        Self::new(512 * 1024)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom_filter_basic() {
        let mut bloom = BloomFilter::new(1024);

        // Insert some k-mers
        bloom.insert(12345);
        bloom.insert(67890);
        bloom.insert(11111);

        // Check inserted items (should return true)
        assert!(bloom.check(12345));
        assert!(bloom.check(67890));
        assert!(bloom.check(11111));

        // Check non-inserted items (may return false positives, but likely false)
        // We can't assert false because bloom filters have false positives
    }

    #[test]
    fn test_bloom_filter_filling_factor() {
        let mut bloom = BloomFilter::new(1024); // ~128 item capacity

        // Insert items
        for i in 0..50 {
            bloom.insert(i);
        }

        let ff = bloom.filling_factor();
        assert!(ff > 0.3); // 50 items / 128 capacity = 0.39
        assert!(ff < 0.5);
    }

    #[test]
    fn test_bloom_filter_resize() {
        let mut bloom = BloomFilter::new(1024);

        bloom.insert(12345);
        assert!(bloom.check(12345));

        // Resize (clears all bits)
        bloom.resize(2048);

        // After resize, items are gone (must re-insert)
        // Note: check might still return true due to false positives
        assert_eq!(bloom.num_items, 0);
    }

    #[test]
    fn test_bloom_filter_clear() {
        let mut bloom = BloomFilter::new(1024);

        bloom.insert(12345);
        bloom.clear();

        assert_eq!(bloom.num_items, 0);
    }
}
