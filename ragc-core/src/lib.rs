#![allow(clippy::all)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(private_interfaces)]
//! Core compression and decompression algorithms for the AGC genome compression format.
//!
//! This crate implements the complete AGC compression pipeline with full C++ AGC
//! format compatibility. Archives created by this library can be read by the C++
//! implementation and vice versa.
//!
//! # Features
//!
//! - **Compression** - Create AGC archives from FASTA files
//! - **Decompression** - Extract genomes from AGC archives
//! - **C++ Compatibility** - Bidirectional format interoperability
//! - **Multi-sample support** - Handle multiple genomes in one archive
//! - **LZ differential encoding** - Efficient encoding against reference sequences
//! - **ZSTD compression** - High-ratio compression of segments
//!
//! # Examples
//!
//! ## Compressing genomes
//!
//! ```ignore
//! use ragc_core::{Compressor, CompressorConfig};
//! use std::path::Path;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Create a compressor
//! let config = CompressorConfig::default();
//! let mut compressor = Compressor::new("output.agc", config)?;
//!
//! // Add FASTA files
//! compressor.add_fasta_file("sample1", Path::new("genome1.fasta"))?;
//! compressor.add_fasta_file("sample2", Path::new("genome2.fasta"))?;
//!
//! // Finalize the archive
//! compressor.finalize()?;
//! # Ok(())
//! # }
//! ```
//!
//! ## Decompressing genomes
//!
//! ```no_run
//! use ragc_core::{Decompressor, DecompressorConfig};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Open an archive
//! let config = DecompressorConfig::default();
//! let mut decompressor = Decompressor::open("archive.agc", config)?;
//!
//! // List available samples
//! let samples = decompressor.list_samples();
//! println!("Found {} samples", samples.len());
//!
//! // Extract a sample
//! let contigs = decompressor.get_sample("sample1")?;
//! for (name, sequence) in contigs {
//!     println!(">{}",  name);
//!     // sequence is Vec<u8> with numeric encoding (A=0, C=1, G=2, T=3)
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## Working with k-mers
//!
//! ```
//! use ragc_core::{Kmer, KmerMode};
//!
//! // Create a canonical k-mer
//! let mut kmer = Kmer::new(21, KmerMode::Canonical);
//!
//! // Insert bases (0=A, 1=C, 2=G, 3=T)
//! kmer.insert(0); // A
//! kmer.insert(1); // C
//! kmer.insert(2); // G
//!
//! if kmer.is_full() {
//!     let value = kmer.data();
//!     println!("K-mer value: {}", value);
//! }
//! ```
//!
//! ## Custom compression settings
//!
//! ```ignore
//! use ragc_core::CompressorConfig;
//!
//! let config = CompressorConfig {
//!     kmer_length: 25,        // Use 25-mers instead of default 21
//!     segment_size: 2000,     // Larger segments
//!     min_match_len: 20,      // Minimum LZ match length
//!     verbosity: 2,           // More verbose output
//! };
//! ```
//!
//! # Archive Format
//!
//! The AGC format organizes data into streams:
//!
//! - **file_type_info** - Version and producer metadata
//! - **params** - Compression parameters (k-mer length, segment size)
//! - **splitters** - Singleton k-mers used for segmentation (future)
//! - **seg-NN** or **seg_dNN** - Compressed genome segments
//! - **collection** - Sample and contig metadata
//!
//! # Compatibility
//!
//! This implementation is tested for compatibility with C++ AGC:
//!
//! - Archives created by ragc can be read by C++ AGC
//! - Archives created by C++ AGC can be read by ragc
//! - Format version 3.0 support
//! - SHA256-verified roundtrip testing

pub mod bloom_filter;
pub mod contig_compression;
pub mod contig_iterator;
pub mod decompressor;
pub mod genome_io;
pub mod kmer;
pub mod kmer_extract;
pub mod lz_diff;
pub mod lz_matcher;
pub mod memory_bounded_queue;
pub mod preprocessing;
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

// C++ AGC FFI (always available for byte-identical archives)
#[path = "ffi/agc_index.rs"]
pub mod agc_index_ffi;

#[path = "ffi/agc_compress.rs"]
pub mod agc_compress_ffi;

#[path = "ffi/splitters.rs"]
pub mod splitters_ffi;

#[path = "ffi/segment_helpers.rs"]
pub mod segment_helpers_ffi;

#[path = "ffi/kmer_helpers.rs"]
pub mod kmer_helpers_ffi;

#[path = "ffi/splitter_check.rs"]
pub mod splitter_check_ffi;

#[path = "ffi/segment_boundary.rs"]
pub mod segment_boundary_ffi;

#[path = "ffi/base_validation.rs"]
pub mod base_validation_ffi;

#[path = "ffi/reverse_complement.rs"]
pub mod reverse_complement_ffi;

#[path = "ffi/segment_split.rs"]
pub mod segment_split_ffi;

#[path = "ffi/kmer_pair.rs"]
pub mod kmer_pair_ffi;

#[path = "ffi/preprocessing.rs"]
pub mod preprocessing_ffi;

#[path = "ffi/find_splitters_in_contig.rs"]
pub mod find_splitters_in_contig_ffi;

#[cfg(feature = "ffi_cost")]
pub mod ragc_ffi {
    extern "C" {
        pub fn agc_cost_vector(
            prefix: i32,
            ref_ptr: *const u8,
            ref_len: usize,
            text_ptr: *const u8,
            text_len: usize,
            min_match_len: u32,
            out_costs: *mut u32,
        ) -> usize;

        pub fn agc_best_split(
            left_ref: *const u8,
            left_len: usize,
            right_ref: *const u8,
            right_len: usize,
            text_ptr: *const u8,
            text_len: usize,
            min_match_len: u32,
            k: u32,
            front_lt_mid: i32,
            mid_lt_back: i32,
            should_reverse: i32,
            out_best_pos: *mut u32,
            out_seg2_start: *mut u32,
            out_should_split: *mut i32,
        ) -> i32;

        pub fn agc_find_middle(
            front_list: *const u64,
            n_front: usize,
            back_list: *const u64,
            n_back: usize,
            out_middle: *mut u64,
        ) -> i32;

        pub fn agc_decide_split(
            front_list: *const u64,
            n_front: usize,
            back_list: *const u64,
            n_back: usize,
            left_ref: *const u8,
            left_len: usize,
            right_ref: *const u8,
            right_len: usize,
            text_ptr: *const u8,
            text_len: usize,
            front_kmer: u64,
            back_kmer: u64,
            min_match_len: u32,
            k: u32,
            should_reverse: i32,
            out_has_middle: *mut i32,
            out_middle: *mut u64,
            out_best_pos: *mut u32,
            out_seg2_start: *mut u32,
            out_should_split: *mut i32,
        ) -> i32;

        // Grouping Engine FFI
        pub fn agc_grouping_engine_create(k: u32, start_group_id: u32) -> *mut std::ffi::c_void;
        pub fn agc_grouping_engine_destroy(engine: *mut std::ffi::c_void);
        pub fn agc_grouping_engine_register(engine: *mut std::ffi::c_void, kmer_front: u64, kmer_back: u64, group_id: u32);
        pub fn agc_grouping_engine_find_middle(engine: *mut std::ffi::c_void, front: u64, back: u64, out_middle: *mut u64) -> i32;
        pub fn agc_grouping_engine_group_exists(engine: *mut std::ffi::c_void, kmer_front: u64, kmer_back: u64) -> i32;
        pub fn agc_grouping_engine_get_group_id(engine: *mut std::ffi::c_void, kmer_front: u64, kmer_back: u64) -> u32;
        pub fn agc_grouping_engine_alloc_id(engine: *mut std::ffi::c_void) -> u32;
    }

    pub fn cost_vector(prefix: bool, reference: &[u8], text: &[u8], min_match_len: u32) -> Vec<u32> {
        unsafe {
            let mut out = vec![0u32; text.len()];
            let _ = agc_cost_vector(
                if prefix { 1 } else { 0 },
                reference.as_ptr(),
                reference.len(),
                text.as_ptr(),
                text.len(),
                min_match_len,
                out.as_mut_ptr(),
            );
            out
        }
    }

    pub fn best_split(left_ref: &[u8], right_ref: &[u8], text: &[u8], min_match_len: u32, k: u32, front_lt_mid: bool, mid_lt_back: bool, should_reverse: bool) -> Option<(usize, usize, bool)> {
        unsafe {
            let mut best: u32 = 0;
            let mut seg2: u32 = 0;
            let mut should: i32 = 0;
            let ok = agc_best_split(
                left_ref.as_ptr(), left_ref.len(),
                right_ref.as_ptr(), right_ref.len(),
                text.as_ptr(), text.len(),
                min_match_len,
                k,
                if front_lt_mid {1} else {0},
                if mid_lt_back {1} else {0},
                if should_reverse {1} else {0},
                &mut best as *mut u32,
                &mut seg2 as *mut u32,
                &mut should as *mut i32,
            );
            if ok != 0 { Some((best as usize, seg2 as usize, should != 0)) } else { None }
        }
    }

    pub fn find_middle(front_neighbors: &[u64], back_neighbors: &[u64]) -> Option<u64> {
        unsafe {
            let mut out: u64 = 0;
            let ok = agc_find_middle(
                front_neighbors.as_ptr(), front_neighbors.len(),
                back_neighbors.as_ptr(), back_neighbors.len(),
                &mut out as *mut u64,
            );
            if ok != 0 { Some(out) } else { None }
        }
    }

    pub fn decide_split(
        front_neighbors: &[u64],
        back_neighbors: &[u64],
        left_ref: &[u8],
        right_ref: &[u8],
        text: &[u8],
        front_kmer: u64,
        back_kmer: u64,
        min_match_len: u32,
        k: u32,
        should_reverse: bool,
    ) -> Option<(bool, u64, usize, usize, bool)> {
        unsafe {
            let mut has_mid: i32 = 0;
            let mut middle: u64 = 0;
            let mut best: u32 = 0;
            let mut seg2: u32 = 0;
            let mut should: i32 = 0;
            let ok = agc_decide_split(
                front_neighbors.as_ptr(), front_neighbors.len(),
                back_neighbors.as_ptr(), back_neighbors.len(),
                left_ref.as_ptr(), left_ref.len(),
                right_ref.as_ptr(), right_ref.len(),
                text.as_ptr(), text.len(),
                front_kmer, back_kmer,
                min_match_len, k,
                if should_reverse {1} else {0},
                &mut has_mid as *mut i32,
                &mut middle as *mut u64,
                &mut best as *mut u32,
                &mut seg2 as *mut u32,
                &mut should as *mut i32,
            );
            if ok != 0 {
                Some((has_mid != 0, middle, best as usize, seg2 as usize, should != 0))
            } else { None }
        }
    }

    /// Safe Rust wrapper for C++ GroupingEngine
    pub struct GroupingEngine {
        ptr: *mut std::ffi::c_void,
    }

    impl GroupingEngine {
        pub fn new(k: u32, start_group_id: u32) -> Self {
            unsafe {
                Self {
                    ptr: agc_grouping_engine_create(k, start_group_id),
                }
            }
        }

        pub fn register_group(&mut self, kmer_front: u64, kmer_back: u64, group_id: u32) {
            unsafe {
                agc_grouping_engine_register(self.ptr, kmer_front, kmer_back, group_id);
            }
        }

        pub fn find_middle(&self, front: u64, back: u64) -> Option<u64> {
            unsafe {
                let mut out: u64 = 0;
                let ok = agc_grouping_engine_find_middle(self.ptr, front, back, &mut out as *mut u64);
                if ok != 0 { Some(out) } else { None }
            }
        }

        pub fn group_exists(&self, kmer_front: u64, kmer_back: u64) -> bool {
            unsafe {
                agc_grouping_engine_group_exists(self.ptr, kmer_front, kmer_back) != 0
            }
        }

        pub fn get_group_id(&self, kmer_front: u64, kmer_back: u64) -> Option<u32> {
            unsafe {
                let gid = agc_grouping_engine_get_group_id(self.ptr, kmer_front, kmer_back);
                if gid == u32::MAX { None } else { Some(gid) }
            }
        }

        pub fn alloc_group_id(&mut self) -> u32 {
            unsafe {
                agc_grouping_engine_alloc_id(self.ptr)
            }
        }
    }

    impl Drop for GroupingEngine {
        fn drop(&mut self) {
            unsafe {
                agc_grouping_engine_destroy(self.ptr);
            }
        }
    }

    // GroupingEngine is Send + Sync because C++ implementation is thread-safe
    unsafe impl Send for GroupingEngine {}
    unsafe impl Sync for GroupingEngine {}
}

// Re-export commonly used types
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
