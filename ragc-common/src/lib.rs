//! Common data structures and utilities for AGC genome compression.
//!
//! This crate provides the foundational types and utilities used across the ragc project:
//!
//! - **Archive I/O** - Reading and writing AGC archive format
//! - **Collection metadata** - Managing samples, contigs, and segment descriptors
//! - **Variable-length integers** - Space-efficient encoding/decoding
//! - **Hash functions** - MurmurHash implementations for k-mer hashing
//! - **Stream naming** - Archive version-aware stream identification
//!
//! # Examples
//!
//! ## Creating and reading an archive
//!
//! ```no_run
//! use ragc_common::Archive;
//!
//! // Create a new archive for writing
//! let mut archive = Archive::new_writer();
//! archive.open("output.agc").expect("Failed to create archive");
//!
//! // Register a stream and add data
//! let stream_id = archive.register_stream("my_stream");
//! let data = b"Hello, AGC!";
//! archive.add_part(stream_id, data, data.len() as u64).expect("Failed to add data");
//!
//! archive.close().expect("Failed to close archive");
//!
//! // Read it back
//! let mut archive = Archive::new_reader();
//! archive.open("output.agc").expect("Failed to open archive");
//!
//! let stream_id = archive.get_stream_id("my_stream").expect("Stream not found");
//! let (data, _) = archive.get_part_by_id(stream_id, 0).expect("Failed to read data");
//!
//! assert_eq!(&data, b"Hello, AGC!");
//! ```
//!
//! ## Variable-length integer encoding
//!
//! ```
//! use ragc_common::{write_varint, read_varint};
//! use std::io::Cursor;
//!
//! let mut buffer = Vec::new();
//! write_varint(&mut buffer, 12345).expect("Failed to encode");
//!
//! let mut cursor = Cursor::new(&buffer);
//! let (value, bytes_read) = read_varint(&mut cursor).expect("Failed to decode");
//!
//! assert_eq!(value, 12345);
//! ```
//!
//! ## Using hash functions
//!
//! ```
//! use ragc_common::MurMur64Hash;
//!
//! let kmer_value = 0x12345678u64;
//! let hash = MurMur64Hash::hash(kmer_value);
//! ```

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
