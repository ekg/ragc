// AGC Common Library
// Shared types, utilities, and data structures

pub mod archive;
pub mod collection;
pub mod hash;
pub mod stream_naming;
pub mod types;
pub mod varint;

// Re-export commonly used types
pub use types::{Base, Contig, PackedBlock};
pub use types::{
    AGC_FILE_MAJOR, AGC_FILE_MINOR, AGC_VER_BUGFIX, AGC_VER_MAJOR, AGC_VER_MINOR, CONTIG_SEPARATOR,
};

// Re-export hash functions
pub use hash::{MurMur32Hash, MurMur64Hash, MurMurPair64Hash, MurMurStringsHash};

// Re-export varint functions
pub use varint::{
    decode_varint, encode_varint, read_fixed_u64, read_varint, write_fixed_u64, write_varint,
};

// Re-export archive
pub use archive::Archive;

// Re-export collection
pub use collection::{zigzag_decode, zigzag_decode_i64, zigzag_encode, zigzag_encode_i64};
pub use collection::{CollectionV3, CollectionVarInt, ContigDesc, SampleDesc, SegmentDesc};

// Re-export stream naming
pub use stream_naming::{
    int_to_base64, stream_base, stream_delta_ext, stream_delta_name, stream_prefix, stream_ref_ext,
    stream_ref_name,
};
