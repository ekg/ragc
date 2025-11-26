//! AGC archive reader/decompressor
//!
//! This crate provides efficient reading of AGC genome archives.

mod decompressor;
mod lz_decode;
mod segment_decompression;

pub use decompressor::{Decompressor, DecompressorConfig};
