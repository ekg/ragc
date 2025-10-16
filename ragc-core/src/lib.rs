// AGC Core Library
// Core compression and decompression algorithms

pub mod compressor;
pub mod decompressor;
pub mod genome_io;
pub mod kmer;
pub mod kmer_extract;
pub mod lz_diff;
pub mod segment;
pub mod segment_compression;
pub mod splitters;

// Re-export commonly used types
pub use compressor::{Compressor, CompressorConfig};
pub use decompressor::{Decompressor, DecompressorConfig};
pub use genome_io::{GenomeIO, GenomeWriter};
pub use kmer::{
    canonical_kmer, decode_base, encode_base, reverse_complement, reverse_complement_kmer,
};
pub use kmer::{Kmer, KmerMode};
pub use kmer_extract::{enumerate_kmers, find_candidate_kmers, remove_non_singletons};
pub use lz_diff::LZDiff;
pub use segment::{split_at_splitters, Segment};
pub use segment_compression::{compress_segment, compress_segment_with_level, decompress_segment};
pub use splitters::{determine_splitters, find_candidate_kmers_multi, is_splitter};
