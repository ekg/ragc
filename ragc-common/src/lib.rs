// AGC Common Library
// Shared types, utilities, and data structures

pub mod types;
pub mod hash;
pub mod varint;
pub mod archive;
pub mod collection;
pub mod stream_naming;

// Re-export commonly used types
pub use types::{Contig, PackedBlock, Base};
pub use types::{AGC_VER_MAJOR, AGC_VER_MINOR, AGC_VER_BUGFIX, AGC_FILE_MAJOR, AGC_FILE_MINOR, CONTIG_SEPARATOR};

// Re-export hash functions
pub use hash::{MurMur64Hash, MurMurPair64Hash, MurMurStringsHash, MurMur32Hash};

// Re-export varint functions
pub use varint::{write_varint, read_varint, write_fixed_u64, read_fixed_u64, encode_varint, decode_varint};

// Re-export archive
pub use archive::Archive;

// Re-export collection
pub use collection::{CollectionV3, SegmentDesc, ContigDesc, SampleDesc, CollectionVarInt};
pub use collection::{zigzag_encode, zigzag_decode, zigzag_encode_i64, zigzag_decode_i64};

// Re-export stream naming
pub use stream_naming::{int_to_base64, stream_prefix, stream_base, stream_ref_name, stream_delta_name, stream_ref_ext, stream_delta_ext};
