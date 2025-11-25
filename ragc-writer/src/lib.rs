//! AGC archive writer/compressor
//!
//! This crate provides compression functionality for creating AGC genome archives.

#[allow(dead_code, unused_variables)]
pub mod _compressor_streaming_old;
pub mod bloom_filter;
pub mod compressor;
pub mod contig_compression;
pub mod contig_iterator;
pub mod decompressor;
pub mod genome_io;
pub mod kmer;
pub mod kmer_extract;
pub mod lz_diff;
pub mod lz_matcher;
pub mod memory_bounded_queue;
pub mod priority_queue;
pub mod segment;
pub mod segment_buffer;
pub mod segment_compression;
pub mod splitters;
pub mod streaming_compressor_queue;
pub mod task;
pub mod tuple_packing;
pub mod worker;
pub mod zstd_pool;

// Re-export commonly used types
pub use _compressor_streaming_old::{StreamingCompressor, StreamingCompressorConfig};
pub use compressor::{Compressor, CompressorConfig};
pub use contig_iterator::{MultiFileIterator, PansnFileIterator};
pub use decompressor::{Decompressor, DecompressorConfig};
pub use genome_io::{GenomeIO, GenomeWriter};
pub use kmer::{
    canonical_kmer, decode_base, encode_base, reverse_complement, reverse_complement_kmer,
};
pub use kmer::{Kmer, KmerMode};
pub use kmer_extract::{enumerate_kmers, find_candidate_kmers, remove_non_singletons};
pub use lz_diff::LZDiff;
pub use memory_bounded_queue::MemoryBoundedQueue;
pub use segment::{split_at_splitters, split_at_splitters_with_size, Segment};
pub use segment_compression::{
    compress_reference_segment, compress_segment, compress_segment_configured, decompress_segment,
    decompress_segment_with_marker,
};
pub use splitters::{
    determine_splitters, determine_splitters_streaming, find_candidate_kmers_multi, is_splitter,
};
pub use streaming_compressor_queue::{QueueStats, StreamingQueueCompressor, StreamingQueueConfig};
pub use worker::create_agc_archive;

// Re-export core types from ragc-core
pub use ragc_core::{
    Archive, CollectionV3, Contig, ContigDesc, SampleDesc, SegmentDesc,
    AGC_FILE_MAJOR, AGC_FILE_MINOR, CONTIG_SEPARATOR,
};
