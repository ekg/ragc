// Hash functions for AGC
// Rust equivalent of MurMurHash implementations in src/common/utils.h

/// MurMurHash3 for 64-bit values
/// Direct equivalent of MurMur64Hash in C++
#[derive(Debug, Clone, Copy)]
pub struct MurMur64Hash;

impl MurMur64Hash {
    #[inline]
    pub fn hash(mut h: u64) -> u64 {
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51afd7ed558ccd_u64);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ceb9fe1a85ec53_u64);
        h ^= h >> 33;
        h
    }
}

impl std::hash::Hasher for MurMur64Hash {
    fn finish(&self) -> u64 {
        0 // Not used in this context
    }

    fn write(&mut self, _bytes: &[u8]) {
        // Not implemented for this use case
    }
}

/// MurMurHash3 for pairs of 64-bit values
/// Direct equivalent of MurMurPair64Hash in C++
#[derive(Debug, Clone, Copy)]
pub struct MurMurPair64Hash;

impl MurMurPair64Hash {
    #[inline]
    pub fn hash(first: u64, second: u64) -> u64 {
        let mut h = first;

        h ^= h >> 33;
        h = h.wrapping_mul(0xff51afd7ed558ccd_u64);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ceb9fe1a85ec53_u64);
        h ^= h >> 33;

        h ^= second;

        h ^= h >> 33;
        h = h.wrapping_mul(0xff51afd7ed558ccd_u64);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ceb9fe1a85ec53_u64);
        h ^= h >> 33;

        h
    }
}

/// MurMurHash3 for strings
/// Direct equivalent of MurMurStringsHash in C++
#[derive(Debug, Clone, Copy)]
pub struct MurMurStringsHash;

impl MurMurStringsHash {
    const C1: u64 = 0x87c37b91114253d5_u64;
    const C2: u64 = 0x4cf5ad432745937f_u64;

    #[inline]
    fn rotl64(x: u64, r: i8) -> u64 {
        x.rotate_left(r as u32)
    }

    #[inline]
    fn fmix64(mut k: u64) -> u64 {
        k ^= k >> 33;
        k = k.wrapping_mul(0xff51afd7ed558ccd_u64);
        k ^= k >> 33;
        k = k.wrapping_mul(0xc4ceb9fe1a85ec53_u64);
        k ^= k >> 33;
        k
    }

    #[inline]
    fn load64(data: &[u8], offset: usize) -> u64 {
        let mut x = 0u64;
        for i in 0..8 {
            if offset + i < data.len() {
                x = (x << 8) | (data[offset + i] as u64);
            }
        }
        x
    }

    pub fn hash(s: &str) -> u64 {
        let data = s.as_bytes();
        let mut h1 = 0u64;
        let mut h2 = 0u64;

        let n_blocks = data.len() / 16;

        // Process 16-byte blocks
        for i in 0..n_blocks {
            let offset = i * 16;
            let mut k1 = Self::load64(data, offset);
            let mut k2 = Self::load64(data, offset + 8);

            k1 = k1.wrapping_mul(Self::C1);
            k1 = Self::rotl64(k1, 31);
            k1 = k1.wrapping_mul(Self::C2);
            h1 ^= k1;

            h1 = Self::rotl64(h1, 27);
            h1 = h1.wrapping_add(h2);
            h1 = h1.wrapping_mul(5).wrapping_add(0x52dce729);

            k2 = k2.wrapping_mul(Self::C2);
            k2 = Self::rotl64(k2, 33);
            k2 = k2.wrapping_mul(Self::C1);
            h2 ^= k2;

            h2 = Self::rotl64(h2, 31);
            h2 = h2.wrapping_add(h1);
            h2 = h2.wrapping_mul(5).wrapping_add(0x38495ab5);
        }

        // Process remaining bytes
        let tail = data.len() % 16;
        let tail_offset = n_blocks * 16;

        let mut k1 = 0u64;
        let mut k2 = 0u64;

        if tail >= 15 {
            k2 ^= (data[tail_offset + 14] as u64) << 48;
        }
        if tail >= 14 {
            k2 ^= (data[tail_offset + 13] as u64) << 40;
        }
        if tail >= 13 {
            k2 ^= (data[tail_offset + 12] as u64) << 32;
        }
        if tail >= 12 {
            k2 ^= (data[tail_offset + 11] as u64) << 24;
        }
        if tail >= 11 {
            k2 ^= (data[tail_offset + 10] as u64) << 16;
        }
        if tail >= 10 {
            k2 ^= (data[tail_offset + 9] as u64) << 8;
        }
        if tail >= 9 {
            k2 ^= data[tail_offset + 8] as u64;
            k2 = k2.wrapping_mul(Self::C2);
            k2 = Self::rotl64(k2, 33);
            k2 = k2.wrapping_mul(Self::C1);
            h2 ^= k2;
        }
        if tail >= 8 {
            k1 ^= (data[tail_offset + 7] as u64) << 56;
        }
        if tail >= 7 {
            k1 ^= (data[tail_offset + 6] as u64) << 48;
        }
        if tail >= 6 {
            k1 ^= (data[tail_offset + 5] as u64) << 40;
        }
        if tail >= 5 {
            k1 ^= (data[tail_offset + 4] as u64) << 32;
        }
        if tail >= 4 {
            k1 ^= (data[tail_offset + 3] as u64) << 24;
        }
        if tail >= 3 {
            k1 ^= (data[tail_offset + 2] as u64) << 16;
        }
        if tail >= 2 {
            k1 ^= (data[tail_offset + 1] as u64) << 8;
        }
        if tail >= 1 {
            k1 ^= data[tail_offset] as u64;
            k1 = k1.wrapping_mul(Self::C1);
            k1 = Self::rotl64(k1, 31);
            k1 = k1.wrapping_mul(Self::C2);
            h1 ^= k1;
        }

        h1 ^= data.len() as u64;
        h2 ^= data.len() as u64;

        h1 = h1.wrapping_add(h2);
        h2 = h2.wrapping_add(h1);

        h1 = Self::fmix64(h1);
        h2 = Self::fmix64(h2);

        h1 = h1.wrapping_add(h2);
        h2 = h2.wrapping_add(h1);

        h1 ^ h2
    }
}

/// MurMurHash32 for 32-bit values
#[derive(Debug, Clone, Copy)]
pub struct MurMur32Hash;

impl MurMur32Hash {
    #[inline]
    pub fn hash(mut h: u32) -> u32 {
        h ^= h >> 16;
        h = h.wrapping_mul(0x85ebca6b);
        h ^= h >> 13;
        h = h.wrapping_mul(0xc2b2ae35);
        h ^= h >> 16;
        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_murmurhash64() {
        // Test basic properties (actual values will be validated against C++)
        assert_eq!(MurMur64Hash::hash(0), 0);
        assert_ne!(MurMur64Hash::hash(1), 0);
        assert_ne!(MurMur64Hash::hash(42), 0);
    }

    #[test]
    fn test_murmurhash_pair64() {
        // Test basic properties
        let h1 = MurMurPair64Hash::hash(0, 0);
        assert_eq!(h1, 0);

        let h2 = MurMurPair64Hash::hash(1, 2);
        assert_ne!(h2, 0);
    }

    #[test]
    fn test_murmurhash_strings() {
        // Test basic properties
        let h1 = MurMurStringsHash::hash("");
        assert_eq!(h1, 0);

        let h2 = MurMurStringsHash::hash("hello");
        assert_ne!(h2, 0);

        let h3 = MurMurStringsHash::hash("world");
        assert_ne!(h3, 0);
        assert_ne!(h2, h3);
    }

    #[test]
    fn test_murmurhash32() {
        // Test basic properties
        assert_eq!(MurMur32Hash::hash(0), 0);
        assert_ne!(MurMur32Hash::hash(1), 0);
        assert_ne!(MurMur32Hash::hash(42), 0);
    }
}
