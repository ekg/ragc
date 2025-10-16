// AGC Core Library
// Core compression and decompression algorithms

pub mod kmer;
pub mod genome_io;
pub mod kmer_extract;
pub mod splitters;
pub mod segment;
pub mod lz_diff;
pub mod segment_compression;
pub mod compressor;
pub mod decompressor;

// Re-export commonly used types
pub use kmer::{Kmer, KmerMode};
pub use kmer::{encode_base, decode_base, canonical_kmer, reverse_complement, reverse_complement_kmer};
pub use genome_io::{GenomeIO, GenomeWriter};
pub use kmer_extract::{enumerate_kmers, remove_non_singletons, find_candidate_kmers};
pub use splitters::{determine_splitters, find_candidate_kmers_multi, is_splitter};
pub use segment::{Segment, split_at_splitters};
pub use lz_diff::LZDiff;
pub use segment_compression::{compress_segment, decompress_segment, compress_segment_with_level};
pub use compressor::{Compressor, CompressorConfig};
pub use decompressor::{Decompressor, DecompressorConfig};
