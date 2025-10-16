// K-mer operations
// Rust equivalent of src/core/kmer.h

use ragc_common::Base;

/// K-mer mode (orientation)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerMode {
    Direct,
    RevComp,
    Canonical,
}

/// K-mer structure that can store both direct and reverse complement representations
/// K-mers are stored in a compact 64-bit representation, with 2 bits per base.
#[derive(Debug, Clone)]
pub struct Kmer {
    kmer_dir: u64,     // Direct k-mer
    kmer_rc: u64,      // Reverse complement k-mer
    cur_size: u32,     // Current size of k-mer
    max_size: u32,     // Maximum k-mer size
    variant: KmerMode, // Which variant to use
    mask: u64,         // Mask for maximum k-mer size
    shift: u32,        // Shift amount
}

impl Kmer {
    /// Create a new empty k-mer with specified maximum size and mode
    pub fn new(max_size: u32, variant: KmerMode) -> Self {
        let shift = 64 - 2 * max_size;
        let mask = (!0u64) << shift;

        Self {
            kmer_dir: 0,
            kmer_rc: 0,
            cur_size: 0,
            max_size,
            variant,
            mask,
            shift,
        }
    }

    /// Create a k-mer from existing direct and reverse complement values
    pub fn from_values(
        kmer_dir: u64,
        kmer_rc: u64,
        max_size: u32,
        cur_size: u32,
        variant: KmerMode,
    ) -> Self {
        let shift = 64 - 2 * max_size;
        let mask = (!0u64) << shift;

        Self {
            kmer_dir,
            kmer_rc,
            cur_size,
            max_size,
            variant,
            mask,
            shift,
        }
    }

    /// Reset the k-mer to empty state
    pub fn reset(&mut self) {
        self.kmer_dir = 0;
        self.kmer_rc = 0;
        self.cur_size = 0;
    }

    /// Insert a symbol (base) at the end of the k-mer (canonical mode)
    #[inline]
    pub fn insert_canonical(&mut self, symbol: u64) {
        // Reverse complement code
        self.kmer_rc >>= 2;
        self.kmer_rc += reverse_complement(symbol) << 62;
        self.kmer_rc &= self.mask;

        // Direct code
        if self.cur_size == self.max_size {
            self.kmer_dir <<= 2;
            self.kmer_dir += symbol << self.shift;
        } else {
            self.cur_size += 1;
            self.kmer_dir += symbol << (64 - 2 * self.cur_size);
        }
    }

    /// Insert a symbol based on the current mode
    #[inline]
    pub fn insert(&mut self, symbol: u64) {
        match self.variant {
            KmerMode::Direct => self.insert_direct(symbol),
            KmerMode::RevComp => self.insert_rev_comp(symbol),
            KmerMode::Canonical => self.insert_canonical(symbol),
        }
    }

    /// Insert a symbol (direct mode)
    #[inline]
    fn insert_direct(&mut self, symbol: u64) {
        if self.cur_size == self.max_size {
            self.kmer_dir <<= 2;
            self.kmer_dir += symbol << self.shift;
        } else {
            self.cur_size += 1;
            self.kmer_dir += symbol << (64 - 2 * self.cur_size);
        }
    }

    /// Insert a symbol (reverse complement mode)
    #[inline]
    fn insert_rev_comp(&mut self, symbol: u64) {
        self.kmer_rc >>= 2;
        self.kmer_rc += reverse_complement(symbol) << 62;
        self.kmer_rc &= self.mask;

        if self.cur_size < self.max_size {
            self.cur_size += 1;
        }
    }

    /// Get the k-mer data based on current mode
    #[inline]
    pub fn data(&self) -> u64 {
        match self.variant {
            KmerMode::Direct => self.kmer_dir,
            KmerMode::RevComp => self.kmer_rc,
            KmerMode::Canonical => self.kmer_dir.min(self.kmer_rc),
        }
    }

    /// Get the canonical k-mer (minimum of direct and reverse complement)
    #[inline]
    pub fn data_canonical(&self) -> u64 {
        self.kmer_dir.min(self.kmer_rc)
    }

    /// Get the direct k-mer
    #[inline]
    pub fn data_dir(&self) -> u64 {
        match self.variant {
            KmerMode::Direct | KmerMode::Canonical => self.kmer_dir,
            KmerMode::RevComp => 0,
        }
    }

    /// Get the reverse complement k-mer
    #[inline]
    pub fn data_rc(&self) -> u64 {
        match self.variant {
            KmerMode::RevComp | KmerMode::Canonical => self.kmer_rc,
            KmerMode::Direct => 0,
        }
    }

    /// Check if k-mer is full (reached maximum size)
    #[inline]
    pub fn is_full(&self) -> bool {
        self.cur_size == self.max_size
    }

    /// Get current size
    #[inline]
    pub fn get_cur_size(&self) -> u32 {
        self.cur_size
    }

    /// Get maximum size
    #[inline]
    pub fn get_max_size(&self) -> u32 {
        self.max_size
    }

    /// Get a symbol at a specific position
    #[inline]
    pub fn get_symbol(&self, pos: u32) -> u64 {
        match self.variant {
            KmerMode::Direct | KmerMode::Canonical => {
                let sym_shift = 62 - 2 * pos;
                (self.kmer_dir >> sym_shift) & 3
            }
            KmerMode::RevComp => {
                let sym_shift = 64 - 2 * self.cur_size + 2 * pos;
                (self.kmer_rc >> sym_shift) & 3
            }
        }
    }

    /// Swap direct and reverse complement representations
    pub fn swap_dir_rc(&mut self) {
        if self.variant == KmerMode::Canonical {
            std::mem::swap(&mut self.kmer_dir, &mut self.kmer_rc);
        }
    }

    /// Check if the k-mer is in direct orientation (for canonical mode)
    pub fn is_dir_oriented(&self) -> bool {
        if self.variant == KmerMode::Canonical {
            self.kmer_dir <= self.kmer_rc
        } else {
            false
        }
    }
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        match self.variant {
            KmerMode::Direct | KmerMode::Canonical => self.kmer_dir == other.kmer_dir,
            KmerMode::RevComp => self.kmer_rc == other.kmer_rc,
        }
    }
}

/// Compute reverse complement of a 2-bit encoded base
/// A(0) <-> T(3), C(1) <-> G(2)
#[inline]
pub const fn reverse_complement(base: u64) -> u64 {
    match base {
        0 => 3, // A -> T
        1 => 2, // C -> G
        2 => 1, // G -> C
        3 => 0, // T -> A
        _ => 4, // Invalid
    }
}

/// Encode a character to a 2-bit value (A=0, C=1, G=2, T=3)
#[inline]
pub fn encode_base(c: char) -> Option<u64> {
    Base::from_char(c).map(|b| b as u64)
}

/// Decode a 2-bit value to a character (0=A, 1=C, 2=G, 3=T)
#[inline]
pub fn decode_base(val: u64) -> Option<char> {
    Base::from_u8(val as u8).map(|b| b.to_char())
}

/// Extract a k-mer from a sequence at a given position
pub fn extract_kmer(sequence: &[u8], pos: usize, k: usize, mode: KmerMode) -> Option<Kmer> {
    if pos + k > sequence.len() {
        return None;
    }

    let mut kmer = Kmer::new(k as u32, mode);

    for i in 0..k {
        let base = encode_base(sequence[pos + i] as char)?;
        kmer.insert(base);
    }

    Some(kmer)
}

/// Compute the canonical k-mer (minimum of k-mer and its reverse complement)
#[inline]
pub fn canonical_kmer(kmer: u64, k: u32) -> u64 {
    let rc = reverse_complement_kmer(kmer, k);
    kmer.min(rc)
}

/// Compute reverse complement of an entire k-mer
pub fn reverse_complement_kmer(kmer: u64, k: u32) -> u64 {
    let mut result = 0u64;
    let shift = 64 - 2 * k;

    for i in 0..k {
        let base = (kmer >> (shift + 2 * i)) & 3;
        let rc_base = reverse_complement(base);
        result |= rc_base << (shift + 2 * (k - 1 - i));
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(0), 3); // A -> T
        assert_eq!(reverse_complement(1), 2); // C -> G
        assert_eq!(reverse_complement(2), 1); // G -> C
        assert_eq!(reverse_complement(3), 0); // T -> A
    }

    #[test]
    fn test_encode_decode() {
        assert_eq!(encode_base('A'), Some(0));
        assert_eq!(encode_base('C'), Some(1));
        assert_eq!(encode_base('G'), Some(2));
        assert_eq!(encode_base('T'), Some(3));

        assert_eq!(decode_base(0), Some('A'));
        assert_eq!(decode_base(1), Some('C'));
        assert_eq!(decode_base(2), Some('G'));
        assert_eq!(decode_base(3), Some('T'));
    }

    #[test]
    fn test_kmer_insertion() {
        let mut kmer = Kmer::new(4, KmerMode::Canonical);

        // Insert ACGT
        kmer.insert(0); // A
        kmer.insert(1); // C
        kmer.insert(2); // G
        kmer.insert(3); // T

        assert!(kmer.is_full());
        assert_eq!(kmer.get_cur_size(), 4);
    }

    #[test]
    fn test_kmer_data() {
        let mut kmer = Kmer::new(4, KmerMode::Direct);

        // Insert AAAA (all zeros)
        for _ in 0..4 {
            kmer.insert(0);
        }

        // AAAA should be 0x0000... in the top bits
        assert_eq!(kmer.data_dir() >> 56, 0);
    }

    #[test]
    fn test_canonical_kmer() {
        // Test that ACGT and its reverse complement ACGT give the same canonical
        let mut kmer1 = Kmer::new(4, KmerMode::Canonical);
        kmer1.insert(0); // A
        kmer1.insert(1); // C
        kmer1.insert(2); // G
        kmer1.insert(3); // T

        let canonical = kmer1.data_canonical();
        // The canonical should be the smaller of forward and RC
        assert!(canonical == kmer1.kmer_dir || canonical == kmer1.kmer_rc);
        assert_eq!(canonical, kmer1.kmer_dir.min(kmer1.kmer_rc));
    }
}
