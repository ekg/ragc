// AGC Core Types and Constants
// Rust equivalent of src/common/defs.h

/// Version information
pub const AGC_VER_MAJOR: u32 = 3;
pub const AGC_VER_MINOR: u32 = 2;
pub const AGC_VER_BUGFIX: u32 = 1;
pub const AGC_VER_BUILD: &str = "20241125.1";

/// Archive file format version
pub const AGC_FILE_MAJOR: u32 = 3;
pub const AGC_FILE_MINOR: u32 = 0;

/// Contig separator byte in packed-contig mode (0xFF)
/// C++ uses this to pack multiple contigs into a single compressed part
pub const CONTIG_SEPARATOR: u8 = 0xFF;

/// Full version string
pub fn agc_version() -> String {
    format!(
        "AGC (Assembled Genomes Compressor) v. {}.{}.{} [build {}]",
        AGC_VER_MAJOR, AGC_VER_MINOR, AGC_VER_BUGFIX, AGC_VER_BUILD
    )
}

/// A contig (genome sequence) represented as a vector of bytes
/// Each byte typically represents a nucleotide base
pub type Contig = Vec<u8>;

/// A packed/compressed data block
pub type PackedBlock = Vec<u8>;

/// Nucleotide encoding (A=0, C=1, G=2, T=3)
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    /// Convert a character to a base (A/C/G/T -> 0/1/2/3)
    #[inline]
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' => Some(Base::T),
            _ => None,
        }
    }

    /// Convert a base to a character (0/1/2/3 -> A/C/G/T)
    #[inline]
    pub fn to_char(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
        }
    }

    /// Get the complement base
    #[inline]
    pub fn complement(self) -> Self {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
        }
    }

    /// Convert u8 value to Base
    #[inline]
    pub fn from_u8(val: u8) -> Option<Self> {
        match val {
            0 => Some(Base::A),
            1 => Some(Base::C),
            2 => Some(Base::G),
            3 => Some(Base::T),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_encoding() {
        assert_eq!(Base::from_char('A'), Some(Base::A));
        assert_eq!(Base::from_char('C'), Some(Base::C));
        assert_eq!(Base::from_char('G'), Some(Base::G));
        assert_eq!(Base::from_char('T'), Some(Base::T));
        assert_eq!(Base::from_char('N'), None);
    }

    #[test]
    fn test_base_complement() {
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::T.complement(), Base::A);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::G.complement(), Base::C);
    }

    #[test]
    fn test_version_string() {
        let version = agc_version();
        assert!(version.contains("3.2.1"));
        assert!(version.contains("20241125.1"));
    }
}
