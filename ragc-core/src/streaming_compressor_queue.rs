// Queue-based streaming compressor API
// Provides simple push() interface with automatic backpressure and constant memory usage

use crate::lz_diff::LZDiff;
use crate::memory_bounded_queue::MemoryBoundedQueue;
use crate::segment::{split_at_splitters_with_size, MISSING_KMER};
use crate::splitters::determine_splitters;
use anyhow::{Context, Result};
use ragc_common::{Archive, CollectionV3, Contig, CONTIG_SEPARATOR};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

/// MurmurHash64A implementation matching C++ AGC's MurMur64Hash
/// This is the same hash function used by C++ AGC for fallback filtering
fn murmur_hash_64a(key: u64) -> u64 {
    const M: u64 = 0xc6a4a7935bd1e995;
    const R: u32 = 47;

    let mut h: u64 = 0xc70f6907u64.wrapping_mul(M);
    let mut k = key;

    k = k.wrapping_mul(M);
    k ^= k >> R;
    k = k.wrapping_mul(M);
    h ^= k;
    h = h.wrapping_mul(M);

    h ^= h >> R;
    h = h.wrapping_mul(M);
    h ^= h >> R;

    h
}

/// Fallback k-mer filter matching C++ AGC's kmer_filter_t
/// Used to select a fraction of k-mers for fallback grouping
#[derive(Debug, Clone)]
struct FallbackFilter {
    /// Threshold for hash comparison (0 = disabled, u64::MAX = all pass)
    threshold: u64,
    /// Random seed for hash mixing (matches C++ AGC's rnd constant)
    rnd: u64,
}

impl FallbackFilter {
    /// Create a new fallback filter with the given fraction
    /// Matches C++ AGC's kmer_filter_t constructor
    fn new(fraction: f64) -> Self {
        let threshold = if fraction == 0.0 {
            0
        } else {
            (u64::MAX as f64 * fraction) as u64
        };
        Self {
            threshold,
            rnd: 0xD73F8BF11046C40E, // Matches C++ AGC constant
        }
    }

    /// Check if the filter is enabled (fraction > 0)
    fn is_enabled(&self) -> bool {
        self.threshold != 0
    }

    /// Check if a k-mer passes the filter
    /// Matches C++ AGC's kmer_filter_t::operator()
    fn passes(&self, kmer: u64) -> bool {
        (murmur_hash_64a(kmer) ^ self.rnd) < self.threshold
    }
}

/// Configuration for the streaming queue-based compressor
#[derive(Debug, Clone)]
pub struct StreamingQueueConfig {
    /// K-mer length for splitters
    pub k: usize,

    /// Segment size for splitting contigs
    pub segment_size: usize,

    /// Minimum match length for LZ encoding
    pub min_match_len: usize,

    /// ZSTD compression level (1-22)
    pub compression_level: i32,

    /// Number of worker threads
    pub num_threads: usize,

    /// Queue capacity in bytes (default: 2 GB, like C++ AGC)
    pub queue_capacity: usize,

    /// Verbosity level
    pub verbosity: usize,

    /// Adaptive mode: find new splitters for samples that can't be segmented well
    /// (matches C++ AGC -a flag)
    pub adaptive_mode: bool,

    /// Fallback fraction: fraction of minimizers to use for fallback grouping
    /// (matches C++ AGC --fallback-frac parameter, default 0.0)
    pub fallback_frac: f64,
}

impl Default for StreamingQueueConfig {
    fn default() -> Self {
        Self {
            k: 31,
            segment_size: 60_000,
            min_match_len: 20,
            compression_level: 17,
            num_threads: rayon::current_num_threads().max(4),
            queue_capacity: 2 * 1024 * 1024 * 1024, // 2 GB like C++ AGC
            verbosity: 1,
            adaptive_mode: false, // Default matches C++ AGC (adaptive mode off)
            fallback_frac: 0.0,   // Default matches C++ AGC (fallback disabled)
        }
    }
}

/// Task to be processed by workers
/// Note: Contig is type alias for Vec<u8>, so we store the name separately
///
/// Priority ordering matches C++ AGC:
/// - Higher sample_priority first (sample1 > sample2 > sample3...)
/// - Within same sample, lexicographic order on contig_name (ascending)
///
/// NOTE: C++ AGC uses a multimap<pair<priority, cost>, T> where cost=contig.size().
/// Since multimap iterates in ascending key order, smaller names come first.
/// This results in lexicographic ordering: chrI, chrII, chrIII, chrIV, chrIX, chrMT, chrV...
/// RAGC must match this ordering for byte-identical archives.
#[derive(Clone)]
struct ContigTask {
    sample_name: String,
    contig_name: String,
    data: Contig, // Vec<u8>
    sample_priority: i32, // Higher = process first (decreases for each sample)
    sequence: u64, // Insertion order within sample - lower = processed first (FASTA order)
}

// Implement priority ordering for BinaryHeap (max-heap)
// BinaryHeap pops the "greatest" element, so we want:
// - Higher sample_priority = greater (first sample processed first)
// - Lexicographically SMALLER contig_name = greater (to be popped first)
//
// C++ AGC uses multimap which iterates in ascending order, so "chrI" < "chrIX" < "chrMT" < "chrV"
// To match this with a max-heap, we reverse the contig_name comparison.
impl PartialEq for ContigTask {
    fn eq(&self, other: &Self) -> bool {
        self.sample_priority == other.sample_priority && self.contig_name == other.contig_name
    }
}

impl Eq for ContigTask {}

impl PartialOrd for ContigTask {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ContigTask {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // First compare by sample_priority (higher priority first)
        match self.sample_priority.cmp(&other.sample_priority) {
            std::cmp::Ordering::Equal => {
                // Then by contig_name (lexicographic ascending order)
                // REVERSE comparison: other.contig_name vs self.contig_name
                // This makes lexicographically smaller names "greater" so they're popped first
                // E.g., "chrI" < "chrIX" < "chrMT" < "chrV" in ascending order
                // So we want "chrI" to be popped first from the max-heap
                other.contig_name.cmp(&self.contig_name)
            }
            other_ord => other_ord,
        }
    }
}

/// Segment group identified by flanking k-mers (matching batch mode)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Ord, PartialOrd)]
struct SegmentGroupKey {
    kmer_front: u64,
    kmer_back: u64,
}

/// Pending segment for batch-local processing (before group assignment)
#[derive(Debug, Clone)]
struct PendingSegment {
    key: SegmentGroupKey,
    segment_data: Vec<u8>,
    should_reverse: bool,
    sample_name: String,
    contig_name: String,
    place: usize,
}

/// Buffered segment waiting to be packed
#[derive(Debug, Clone, PartialEq, Eq)]
struct BufferedSegment {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    data: Contig,
    is_rev_comp: bool,
}

// Match C++ AGC sorting order (agc_compressor.h lines 112-119)
// Sort by: sample_name, then contig_name, then seg_part_no
impl PartialOrd for BufferedSegment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BufferedSegment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // First compare by sample_name
        match self.sample_name.cmp(&other.sample_name) {
            std::cmp::Ordering::Equal => {
                // Then by contig_name
                match self.contig_name.cmp(&other.contig_name) {
                    std::cmp::Ordering::Equal => {
                        // Finally by seg_part_no
                        self.seg_part_no.cmp(&other.seg_part_no)
                    }
                    other => other,
                }
            }
            other => other,
        }
    }
}

/// Buffer for a segment group (packs 50 segments together)
struct SegmentGroupBuffer {
    group_id: u32,
    stream_id: usize,                           // Delta stream for packed segments
    ref_stream_id: usize,                       // Reference stream for first segment
    reference_segment: Option<BufferedSegment>, // First segment (reference for LZ encoding)
    segments: Vec<BufferedSegment>, // Up to PACK_CARDINALITY segments (EXCLUDING reference)
    ref_written: bool,              // Whether reference has been written
    segments_written: u32,          // Counter for delta segments written (NOT including reference)
    lz_diff: Option<LZDiff>,        // LZ encoder prepared once with reference, reused for all segments (matches C++ AGC CSegment::lz_diff)
    // CRITICAL: Partial pack persistence to ensure pack alignment with decompression expectations
    // Pack N must contain entries for in_group_ids (N*50)+1 to (N+1)*50
    // These fields persist unique deltas until we have exactly 50 for a complete pack
    pending_deltas: Vec<Vec<u8>>,   // Unique deltas waiting to be written (< 50)
    pending_delta_ids: Vec<u32>,    // in_group_ids for pending deltas (for deduplication)
    raw_placeholder_written: bool,  // Whether raw group placeholder has been written
}

impl SegmentGroupBuffer {
    fn new(group_id: u32, stream_id: usize, ref_stream_id: usize) -> Self {
        Self {
            group_id,
            stream_id,
            ref_stream_id,
            reference_segment: None,
            segments: Vec::new(),
            ref_written: false,
            segments_written: 0,
            lz_diff: None,  // Prepared when reference is written (matches C++ AGC segment.cpp line 43)
            pending_deltas: Vec::new(),
            pending_delta_ids: Vec::new(),
            raw_placeholder_written: false,
        }
    }
}

/// Pack size (C++ AGC default)
const PACK_CARDINALITY: usize = 50;
/// First 16 groups are raw-only (no LZ encoding)
const NO_RAW_GROUPS: u32 = 16;

/// Streaming compressor with queue-based API
///
/// # Example
/// ```no_run
/// use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
/// use std::collections::HashSet;
///
/// # fn main() -> anyhow::Result<()> {
/// let config = StreamingQueueConfig::default();
/// let splitters = HashSet::new(); // Normally from reference
/// let mut compressor = StreamingQueueCompressor::with_splitters(
///     "output.agc",
///     config,
///     splitters
/// )?;
///
/// // Push sequences (blocks when queue is full - automatic backpressure!)
/// # let sequences = vec![("sample1".to_string(), "chr1".to_string(), vec![0u8; 1000])];
/// for (sample, contig_name, data) in sequences {
///     compressor.push(sample, contig_name, data)?;
/// }
///
/// // Finalize - waits for all compression to complete
/// compressor.finalize()?;
/// # Ok(())
/// # }
/// ```
pub struct StreamingQueueCompressor {
    queue: Arc<MemoryBoundedQueue<ContigTask>>,
    workers: Vec<JoinHandle<Result<()>>>,
    collection: Arc<Mutex<CollectionV3>>,
    splitters: Arc<HashSet<u64>>,
    config: StreamingQueueConfig,
    archive: Arc<Mutex<Archive>>,
    segment_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    group_counter: Arc<AtomicU32>, // Starts at 16 for LZ groups
    raw_group_counter: Arc<AtomicU32>, // Round-robin counter for raw groups (0-15)
    reference_sample_name: Arc<Mutex<Option<String>>>, // First sample becomes reference
    // Segment splitting support (Phase 1)
    map_segments: Arc<Mutex<std::collections::HashMap<SegmentGroupKey, u32>>>, // (front, back) -> group_id
    map_segments_terminators: Arc<Mutex<std::collections::HashMap<u64, Vec<u64>>>>, // kmer -> [connected kmers]

    // FFI Grouping Engine - C++ AGC-compatible group assignment
    #[cfg(feature = "cpp_agc")]
    grouping_engine: Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,

    // Persistent reference segment storage (matches C++ AGC v_segments)
    // Stores reference segment data even after groups are flushed, enabling LZ cost estimation
    // for subsequent samples (fixes multi-sample group fragmentation bug)
    reference_segments: Arc<Mutex<std::collections::HashMap<u32, Vec<u8>>>>, // group_id -> reference segment data

    // Reference orientation tracking - stores is_rev_comp for each group's reference segment
    // When a delta segment joins an existing group, it MUST use the same orientation as the reference
    // to ensure LZ encoding works correctly (fixes ZERO_MATCH bug in Case 3 terminator segments)
    reference_orientations: Arc<Mutex<std::collections::HashMap<u32, bool>>>, // group_id -> reference is_rev_comp

    // Track segment splits for renumbering subsequent segments
    // Maps (sample_name, contig_name, original_place) -> number of splits inserted before this position
    split_offsets: Arc<Mutex<std::collections::HashMap<(String, String, usize), usize>>>,

    // Priority assignment for interleaved processing (matches C++ AGC)
    // Higher priority = processed first (sample1 > sample2 > sample3...)
    sample_priorities: Arc<Mutex<std::collections::HashMap<String, i32>>>, // sample_name -> priority

    // Batch-local group assignment (matches C++ AGC m_kmers per-batch behavior)
    // When sample changes, we flush pending segments and clear batch-local state
    current_batch_sample: Arc<Mutex<Option<String>>>, // Current sample being processed
    batch_local_groups: Arc<Mutex<HashMap<SegmentGroupKey, u32>>>, // Batch-local m_kmers equivalent
    batch_local_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>>, // Batch-local terminators - only merged to global at batch end
    pending_batch_segments: Arc<Mutex<Vec<PendingSegment>>>, // Buffer segments until batch boundary
    // Fallback minimizers map for segments with no terminator match (matches C++ AGC map_fallback_minimizers)
    map_fallback_minimizers: Arc<Mutex<HashMap<u64, Vec<(u64, u64)>>>>, // kmer -> [(front, back)] candidate group keys
    next_priority: Arc<Mutex<i32>>, // Decreases for each new sample (starts at i32::MAX)
    next_sequence: Arc<std::sync::atomic::AtomicU64>, // Increases for each contig (FASTA order)

    // Deferred metadata streams - written AFTER segment data (C++ AGC compatibility)
    // C++ AGC writes segment data first, then metadata streams at the end
    deferred_file_type_info: (usize, Vec<u8>),    // (stream_id, data)
    deferred_params: (usize, Vec<u8>),            // (stream_id, data)
    deferred_splitters: (usize, Vec<u8>),         // (stream_id, data)
    deferred_segment_splitters: (usize, Vec<u8>), // (stream_id, data)
}

impl StreamingQueueCompressor {
    /// Create a new streaming compressor with pre-computed splitters
    ///
    /// Use this when you already have splitters (e.g., from a reference genome)
    ///
    /// # Arguments
    /// * `output_path` - Path to output AGC archive
    /// * `config` - Compression configuration
    /// * `splitters` - Pre-computed splitter k-mers
    pub fn with_splitters(
        output_path: impl AsRef<Path>,
        config: StreamingQueueConfig,
        splitters: HashSet<u64>,
    ) -> Result<Self> {
        let output_path = output_path.as_ref();
        let archive_path = output_path.to_string_lossy().to_string();

        if config.verbosity > 0 {
            eprintln!("Initializing streaming compressor...");
            eprintln!(
                "  Queue capacity: {} GB",
                config.queue_capacity / (1024 * 1024 * 1024)
            );
            eprintln!("  Worker threads: {}", config.num_threads);
            eprintln!("  Splitters: {}", splitters.len());
        }

        // Create archive
        let mut archive = Archive::new_writer();
        archive.open(output_path)?;

        // Create collection
        let mut collection = CollectionV3::new();
        collection.set_config(config.segment_size as u32, config.k as u32, None);

        // CRITICAL: Register collection streams FIRST (C++ AGC compatibility)
        // C++ AGC expects collection-samples at stream 0, collection-contigs at 1, collection-details at 2
        collection.prepare_for_compression(&mut archive)?;

        // DEFERRED METADATA STREAMS (C++ AGC compatibility)
        // C++ AGC writes segment data FIRST, then metadata streams at the END.
        // We register streams now but defer writing data until finalize().

        // Prepare file_type_info data (defer write)
        let deferred_file_type_info = {
            let mut data = Vec::new();
            let append_str = |data: &mut Vec<u8>, s: &str| {
                data.extend_from_slice(s.as_bytes());
                data.push(0);
            };

            append_str(&mut data, "producer");
            append_str(&mut data, "ragc");
            append_str(&mut data, "producer_version_major");
            append_str(&mut data, &ragc_common::AGC_FILE_MAJOR.to_string());
            append_str(&mut data, "producer_version_minor");
            append_str(&mut data, &ragc_common::AGC_FILE_MINOR.to_string());
            append_str(&mut data, "producer_version_build");
            append_str(&mut data, "0");
            append_str(&mut data, "file_version_major");
            append_str(&mut data, &ragc_common::AGC_FILE_MAJOR.to_string());
            append_str(&mut data, "file_version_minor");
            append_str(&mut data, &ragc_common::AGC_FILE_MINOR.to_string());
            append_str(&mut data, "comment");
            append_str(
                &mut data,
                &format!(
                    "RAGC v.{}.{}",
                    ragc_common::AGC_FILE_MAJOR,
                    ragc_common::AGC_FILE_MINOR
                ),
            );

            let stream_id = archive.register_stream("file_type_info");
            // DEFERRED: archive.add_part(stream_id, &data, 7) will be called in finalize()
            (stream_id, data)
        };

        // Prepare params data (defer write)
        let deferred_params = {
            let stream_id = archive.register_stream("params");
            let mut data = Vec::new();
            data.extend_from_slice(&(config.k as u32).to_le_bytes());
            data.extend_from_slice(&(config.min_match_len as u32).to_le_bytes());
            data.extend_from_slice(&50u32.to_le_bytes()); // pack_cardinality (default)
            data.extend_from_slice(&(config.segment_size as u32).to_le_bytes());
            // DEFERRED: archive.add_part(stream_id, &data, 0) will be called in finalize()
            (stream_id, data)
        };

        // Prepare empty splitters stream (defer write)
        let deferred_splitters = {
            let stream_id = archive.register_stream("splitters");
            let data = Vec::new();
            // DEFERRED: archive.add_part(stream_id, &data, 0) will be called in finalize()
            (stream_id, data)
        };

        // Prepare empty segment-splitters stream (defer write)
        let deferred_segment_splitters = {
            let stream_id = archive.register_stream("segment-splitters");
            let data = Vec::new();
            // DEFERRED: archive.add_part(stream_id, &data, 0) will be called in finalize()
            (stream_id, data)
        };

        let collection = Arc::new(Mutex::new(collection));
        let archive = Arc::new(Mutex::new(archive));

        // Create memory-bounded queue
        let queue = Arc::new(MemoryBoundedQueue::new(config.queue_capacity));

        let splitters = Arc::new(splitters);

        // Segment grouping for LZ packing (using BTreeMap for better memory efficiency)
        let segment_groups = Arc::new(Mutex::new(BTreeMap::new()));
        let group_counter = Arc::new(AtomicU32::new(NO_RAW_GROUPS)); // Start at 16 (LZ groups), raw groups 0-15 handled separately
        let raw_group_counter = Arc::new(AtomicU32::new(0)); // Round-robin counter for raw groups (0-15)
        let reference_sample_name = Arc::new(Mutex::new(None)); // Shared across all workers

        // Segment splitting support (Phase 1)
        let map_segments: Arc<Mutex<HashMap<SegmentGroupKey, u32>>> = Arc::new(Mutex::new(HashMap::new()));
        let map_segments_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>> = Arc::new(Mutex::new(HashMap::new()));
        let split_offsets: Arc<Mutex<HashMap<(String, String, usize), usize>>> = Arc::new(Mutex::new(HashMap::new()));

        // Persistent reference segment storage (matches C++ AGC v_segments)
        let reference_segments: Arc<Mutex<HashMap<u32, Vec<u8>>>> = Arc::new(Mutex::new(HashMap::new()));

        // Reference orientation tracking (fixes ZERO_MATCH bug in Case 3 terminator segments)
        let reference_orientations: Arc<Mutex<HashMap<u32, bool>>> = Arc::new(Mutex::new(HashMap::new()));

        // FFI Grouping Engine - C++ AGC-compatible group assignment
        #[cfg(feature = "cpp_agc")]
        let grouping_engine = Arc::new(Mutex::new(crate::ragc_ffi::GroupingEngine::new(
            config.k as u32,
            NO_RAW_GROUPS,  // Start group IDs at 16 (raw groups 0-15 handled separately)
        )));

        // Priority tracking for interleaved processing (matches C++ AGC)
        let sample_priorities: Arc<Mutex<HashMap<String, i32>>> = Arc::new(Mutex::new(HashMap::new()));
        let next_priority = Arc::new(Mutex::new(i32::MAX)); // Start high, decrease for each sample
        let next_sequence = Arc::new(std::sync::atomic::AtomicU64::new(0)); // Increases for each contig (FASTA order)

        // Batch-local group assignment (matches C++ AGC m_kmers per-batch behavior)
        let current_batch_sample: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));
        let batch_local_groups: Arc<Mutex<HashMap<SegmentGroupKey, u32>>> = Arc::new(Mutex::new(HashMap::new()));
        let batch_local_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>> = Arc::new(Mutex::new(HashMap::new()));
        let pending_batch_segments: Arc<Mutex<Vec<PendingSegment>>> = Arc::new(Mutex::new(Vec::new()));
        let map_fallback_minimizers: Arc<Mutex<HashMap<u64, Vec<(u64, u64)>>>> = Arc::new(Mutex::new(HashMap::new()));

        // Spawn worker threads
        let mut workers = Vec::new();
        for worker_id in 0..config.num_threads {
            let queue = Arc::clone(&queue);
            let collection = Arc::clone(&collection);
            let splitters = Arc::clone(&splitters);
            let archive = Arc::clone(&archive);
            let segment_groups = Arc::clone(&segment_groups);
            let group_counter = Arc::clone(&group_counter);
            let raw_group_counter = Arc::clone(&raw_group_counter);
            let reference_sample_name = Arc::clone(&reference_sample_name);
            let map_segments = Arc::clone(&map_segments);
            let map_segments_terminators = Arc::clone(&map_segments_terminators);
            let reference_segments = Arc::clone(&reference_segments);
            let reference_orientations = Arc::clone(&reference_orientations);
            let split_offsets = Arc::clone(&split_offsets);
            #[cfg(feature = "cpp_agc")]
            let grouping_engine = Arc::clone(&grouping_engine);
            let current_batch_sample = Arc::clone(&current_batch_sample);
            let batch_local_groups = Arc::clone(&batch_local_groups);
            let batch_local_terminators = Arc::clone(&batch_local_terminators);
            let pending_batch_segments = Arc::clone(&pending_batch_segments);
            let map_fallback_minimizers = Arc::clone(&map_fallback_minimizers);
            let config = config.clone();

            let handle = thread::spawn(move || {
                worker_thread(
                    worker_id,
                    queue,
                    collection,
                    splitters,
                    archive,
                    segment_groups,
                    group_counter,
                    raw_group_counter,
                    reference_sample_name,
                    map_segments,
                    map_segments_terminators,
                    reference_segments,
                    reference_orientations,
                    split_offsets,
                    #[cfg(feature = "cpp_agc")]
                    grouping_engine,
                    current_batch_sample,
                    batch_local_groups,
                    batch_local_terminators,
                    pending_batch_segments,
                    map_fallback_minimizers,
                    config,
                )
            });

            workers.push(handle);
        }

        if config.verbosity > 0 {
            eprintln!("Ready to receive sequences!");
        }

        Ok(Self {
            queue,
            workers,
            collection,
            splitters,
            config,
            archive,
            segment_groups,
            group_counter,
            raw_group_counter,
            reference_sample_name,
            map_segments,
            map_segments_terminators,
            #[cfg(feature = "cpp_agc")]
            grouping_engine,
            reference_segments,
            reference_orientations,
            split_offsets,
            sample_priorities,
            next_priority,
            current_batch_sample,
            batch_local_groups,
            batch_local_terminators,
            pending_batch_segments,
            map_fallback_minimizers,
            next_sequence,
            // Deferred metadata streams (written at end for C++ AGC compatibility)
            deferred_file_type_info,
            deferred_params,
            deferred_splitters,
            deferred_segment_splitters,
        })
    }

    /// Create compressor and determine splitters from first contig
    ///
    /// **Note**: This requires at least one contig to be pushed before workers start.
    /// Consider using `with_splitters()` instead if you have a reference genome.
    pub fn new(output_path: impl AsRef<Path>, config: StreamingQueueConfig) -> Result<Self> {
        // Start with empty splitters - will be determined from first push
        Self::with_splitters(output_path, config, HashSet::new())
    }

    /// Push a contig to the compression queue
    ///
    /// **BLOCKS** if the queue is full (automatic backpressure!)
    ///
    /// # Arguments
    /// * `sample_name` - Name of the sample
    /// * `contig_name` - Name of the contig
    /// * `data` - Contig sequence data (Vec<u8>)
    ///
    /// # Example
    /// ```no_run
    /// # use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
    /// # use std::collections::HashSet;
    /// # let mut compressor = StreamingQueueCompressor::with_splitters("out.agc", StreamingQueueConfig::default(), HashSet::new())?;
    /// compressor.push("sample1".to_string(), "chr1".to_string(), vec![b'A', b'T', b'G', b'C'])?;
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn push(&mut self, sample_name: String, contig_name: String, data: Contig) -> Result<()> {
        // If no splitters yet, determine from this contig
        if self.splitters.is_empty() && self.workers.is_empty() {
            if self.config.verbosity > 0 {
                eprintln!("Determining splitters from first contig...");
            }

            let (splitters, _, _) =
                determine_splitters(&[data.clone()], self.config.k, self.config.segment_size);

            if self.config.verbosity > 0 {
                eprintln!("Found {} splitters", splitters.len());
            }

            // Update splitters and spawn workers
            self.splitters = Arc::new(splitters);

            // Spawn workers now that we have splitters
            for worker_id in 0..self.config.num_threads {
                let queue = Arc::clone(&self.queue);
                let collection = Arc::clone(&self.collection);
                let splitters = Arc::clone(&self.splitters);
                let archive = Arc::clone(&self.archive);
                let segment_groups = Arc::clone(&self.segment_groups);
                let group_counter = Arc::clone(&self.group_counter);
                let raw_group_counter = Arc::clone(&self.raw_group_counter);
                let reference_sample_name = Arc::clone(&self.reference_sample_name);
                let map_segments = Arc::clone(&self.map_segments);
                let map_segments_terminators = Arc::clone(&self.map_segments_terminators);
                let reference_segments = Arc::clone(&self.reference_segments);
                let reference_orientations = Arc::clone(&self.reference_orientations);
                let split_offsets = Arc::clone(&self.split_offsets);
                #[cfg(feature = "cpp_agc")]
                let grouping_engine = Arc::clone(&self.grouping_engine);
                let current_batch_sample = Arc::clone(&self.current_batch_sample);
                let batch_local_groups = Arc::clone(&self.batch_local_groups);
                let batch_local_terminators = Arc::clone(&self.batch_local_terminators);
                let pending_batch_segments = Arc::clone(&self.pending_batch_segments);
                let map_fallback_minimizers = Arc::clone(&self.map_fallback_minimizers);
                let config = self.config.clone();

                let handle = thread::spawn(move || {
                    worker_thread(
                        worker_id,
                        queue,
                        collection,
                        splitters,
                        archive,
                        segment_groups,
                        group_counter,
                        raw_group_counter,
                        reference_sample_name,
                        map_segments,
                        map_segments_terminators,
                        reference_segments,
                        reference_orientations,
                        split_offsets,
                        #[cfg(feature = "cpp_agc")]
                        grouping_engine,
                        current_batch_sample,
                        batch_local_groups,
                        batch_local_terminators,
                        pending_batch_segments,
                        map_fallback_minimizers,
                        config,
                    )
                });

                self.workers.push(handle);
            }

            if self.config.verbosity > 0 {
                eprintln!("Workers spawned and ready!");
            }
        }

        // Register contig in collection
        {
            let mut collection = self.collection.lock().unwrap();
            collection
                .register_sample_contig(&sample_name, &contig_name)
                .context("Failed to register contig")?;
        }

        // Set first sample as reference (multi-file mode)
        {
            let mut ref_sample = self.reference_sample_name.lock().unwrap();
            if ref_sample.is_none() {
                if self.config.verbosity > 0 {
                    eprintln!("Using first sample ({}) as reference", sample_name);
                }
                *ref_sample = Some(sample_name.clone());
            }
        }

        // Calculate task size
        let task_size = data.len();

        // Get sequence number for FASTA ordering (lower = earlier = higher priority)
        let sequence = self.next_sequence.fetch_add(1, std::sync::atomic::Ordering::SeqCst);

        // Get or assign priority for this sample (matches C++ AGC priority queue)
        // Higher priority = processed first (decreases for each new sample)
        let sample_priority = {
            let mut priorities = self.sample_priorities.lock().unwrap();
            *priorities.entry(sample_name.clone()).or_insert_with(|| {
                // First time seeing this sample - assign new priority
                let mut next_p = self.next_priority.lock().unwrap();
                let priority = *next_p;
                *next_p -= 1; // Decrement for next sample (C++ AGC uses --sample_priority)
                priority
            })
        };

        // Create task with priority information
        // NOTE: sequence is used for FASTA ordering (lower = processed first)
        let task = ContigTask {
            sample_name,
            contig_name,
            data,
            sample_priority,
            sequence,
        };

        // Push to queue (BLOCKS if queue is full!)
        // Queue is now a priority queue - highest priority processed first
        self.queue
            .push(task, task_size)
            .context("Failed to push to queue")?;

        Ok(())
    }

    /// Finalize compression
    ///
    /// This will:
    /// 1. Close the queue (no more pushes allowed)
    /// 2. Wait for all worker threads to finish processing
    /// 3. Write metadata to the archive
    /// 4. Close the archive file
    ///
    /// # Example
    /// ```no_run
    /// # use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
    /// # use std::collections::HashSet;
    /// # let mut compressor = StreamingQueueCompressor::with_splitters("out.agc", StreamingQueueConfig::default(), HashSet::new())?;
    /// // ... push sequences ...
    /// compressor.finalize()?;
    /// # Ok::<(), antml:Error>(())
    /// ```
    pub fn drain(&self) -> Result<()> {
        if self.config.verbosity > 0 {
            eprintln!("Draining queue (waiting for {} items to be processed)...", self.queue.len());
        }

        // Wait for queue to empty
        // Poll every 100ms until queue is empty
        while self.queue.len() > 0 {
            std::thread::sleep(std::time::Duration::from_millis(100));
        }

        if self.config.verbosity > 0 {
            eprintln!("Queue drained - all queued contigs processed");
        }

        Ok(())
    }

    pub fn finalize(self) -> Result<()> {
        if self.config.verbosity > 0 {
            eprintln!("Finalizing compression...");
            eprintln!("  Closing queue...");
        }

        // Close queue - no more pushes allowed
        self.queue.close();

        if self.config.verbosity > 0 {
            eprintln!("  Waiting for {} workers to finish...", self.workers.len());
        }

        // Wait for all workers to finish
        for (i, handle) in self.workers.into_iter().enumerate() {
            handle
                .join()
                .expect("Worker thread panicked")
                .with_context(|| format!("Worker {} failed", i))?;
        }

        if self.config.verbosity > 0 {
            eprintln!("All workers finished!");
            eprintln!("Flushing remaining segment packs...");
        }

        // Flush all remaining partial packs
        {
            let mut groups = self.segment_groups.lock().unwrap();
            let num_groups = groups.len();

            for (key, buffer) in groups.iter_mut() {
                // Flush if there are delta segments OR if reference hasn't been written
                if !buffer.segments.is_empty() || !buffer.ref_written {
                    if self.config.verbosity > 1 {
                        eprintln!(
                            "Flushing group {} with {} segments (k-mers: {:#x}, {:#x})",
                            buffer.group_id,
                            buffer.segments.len(),
                            key.kmer_front,
                            key.kmer_back
                        );
                    }
                    flush_pack(buffer, &self.collection, &self.archive, &self.config, &self.reference_segments)
                        .context("Failed to flush remaining pack")?;
                }

                // CRITICAL: Write any remaining pending_deltas as final partial pack
                // This is where partial packs (< 50 entries) are written, ensuring
                // pack boundaries align with decompression expectations
                if !buffer.pending_deltas.is_empty() {
                    use crate::segment_compression::compress_segment_configured;

                    let use_lz_encoding = buffer.group_id >= NO_RAW_GROUPS;
                    let mut packed_data = Vec::new();

                    // Raw groups need placeholder if this is the first pack
                    if !use_lz_encoding && !buffer.raw_placeholder_written {
                        packed_data.push(0x7f);
                        packed_data.push(CONTIG_SEPARATOR);
                    }

                    for delta in buffer.pending_deltas.iter() {
                        packed_data.extend_from_slice(delta);
                        packed_data.push(CONTIG_SEPARATOR);
                    }

                    if !packed_data.is_empty() {
                        let total_raw_size = packed_data.len();
                        let mut compressed = compress_segment_configured(&packed_data, self.config.compression_level)
                            .context("Failed to compress final partial pack")?;
                        compressed.push(0); // Marker 0 = plain ZSTD

                        let mut arch = self.archive.lock().unwrap();
                        if compressed.len() < total_raw_size {
                            arch.add_part(buffer.stream_id, &compressed, total_raw_size as u64)
                                .context("Failed to write compressed final partial pack")?;
                        } else {
                            arch.add_part(buffer.stream_id, &packed_data, 0)
                                .context("Failed to write uncompressed final partial pack")?;
                        }
                    }

                    if self.config.verbosity > 1 {
                        eprintln!(
                            "  Wrote final partial pack for group {} with {} deltas",
                            buffer.group_id, buffer.pending_deltas.len()
                        );
                    }

                    buffer.pending_deltas.clear();
                    buffer.pending_delta_ids.clear();
                }
            }

            if self.config.verbosity > 0 {
                eprintln!("Flushed {} segment groups", num_groups);
            }
        }

        if self.config.verbosity > 0 {
            eprintln!("Writing metadata...");
        }

        // Get total sample count for metadata writing
        let num_samples = {
            let coll = self.collection.lock().unwrap();
            coll.get_no_samples()
        };

        // Write collection metadata to archive
        {
            let mut archive = self.archive.lock().unwrap();
            let mut collection = self.collection.lock().unwrap();

            // DEFERRED METADATA WRITES (C++ AGC compatibility)
            // C++ AGC writes metadata streams AFTER segment data, in this order:
            // 1. params
            // 2. splitters
            // 3. segment-splitters
            // 4. collection metadata (samples, contigs, details)
            // 5. file_type_info
            let (params_stream_id, params_data) = &self.deferred_params;
            archive.add_part(*params_stream_id, params_data, 0)
                .context("Failed to write params")?;

            let (splitters_stream_id, splitters_data) = &self.deferred_splitters;
            archive.add_part(*splitters_stream_id, splitters_data, 0)
                .context("Failed to write splitters")?;

            let (seg_splitters_stream_id, seg_splitters_data) = &self.deferred_segment_splitters;
            archive.add_part(*seg_splitters_stream_id, seg_splitters_data, 0)
                .context("Failed to write segment-splitters")?;

            // Write sample names
            collection
                .store_batch_sample_names(&mut archive)
                .context("Failed to write sample names")?;

            // Write contig names and segment details
            collection
                .store_contig_batch(&mut archive, 0, num_samples)
                .context("Failed to write contig batch")?;

            // Write file_type_info LAST (matches C++ AGC store_file_type_info order)
            let (file_type_info_stream_id, file_type_info_data) = &self.deferred_file_type_info;
            archive.add_part(*file_type_info_stream_id, file_type_info_data, 7)
                .context("Failed to write file_type_info")?;

            if self.config.verbosity > 0 {
                eprintln!("Collection metadata written successfully");
            }

            // Close archive (writes footer)
            archive.close().context("Failed to close archive")?;
        }

        if self.config.verbosity > 0 {
            eprintln!("Compression complete!");
        }

        Ok(())
    }

    /// Get current queue statistics
    pub fn queue_stats(&self) -> QueueStats {
        QueueStats {
            current_size_bytes: self.queue.current_size(),
            current_items: self.queue.len(),
            capacity_bytes: self.queue.capacity(),
            is_closed: self.queue.is_closed(),
        }
    }
}

/// Queue statistics
#[derive(Debug, Clone)]
pub struct QueueStats {
    pub current_size_bytes: usize,
    pub current_items: usize,
    pub capacity_bytes: usize,
    pub is_closed: bool,
}

/// Flush a complete pack of segments (compress, LZ encode, write to archive)
fn flush_pack(
    buffer: &mut SegmentGroupBuffer,
    collection: &Arc<Mutex<CollectionV3>>,
    archive: &Arc<Mutex<Archive>>,
    config: &StreamingQueueConfig,
    reference_segments: &Arc<Mutex<HashMap<u32, Vec<u8>>>>,
) -> Result<()> {
    use crate::segment_compression::{compress_reference_segment, compress_segment_configured};

    // Skip if no segments to write (but still write reference if present)
    if buffer.segments.is_empty() && buffer.ref_written {
        return Ok(());
    }

    let use_lz_encoding = buffer.group_id >= NO_RAW_GROUPS;

    // DO NOT SORT: With file-by-file batch processing (matching C++ AGC's pack_cardinality batching),
    // segments arrive in the correct order already (file order = batch order).
    // C++ AGC sorts WITHIN each batch BEFORE distributing to groups,
    // so segments are already sorted when they arrive at each group.
    // Sorting here would reorder segments from DIFFERENT batches alphabetically, which is wrong!

    // Write reference segment if not already written (first pack for this group)
    // Extract reference from sorted segments (matching C++ AGC: first segment after sort becomes reference)
    // NOTE: Raw groups (0-15) do NOT have a reference - all segments stored raw
    if use_lz_encoding && !buffer.ref_written && !buffer.segments.is_empty() {
        // Remove first segment (alphabetically first) to use as reference
        let ref_seg = buffer.segments.remove(0);

        if std::env::var("RAGC_DEBUG_REF_WRITE").is_ok() {
            eprintln!(
                "DEBUG_REF_WRITE: group={} sample={} contig={} seg={} data_len={} segments_remaining={}",
                buffer.group_id, ref_seg.sample_name, ref_seg.contig_name,
                ref_seg.seg_part_no, ref_seg.data.len(), buffer.segments.len()
            );
        }

        if config.verbosity > 1 {
            eprintln!(
                "  Flushing group {}: reference from {} (chosen from {} sorted segments)",
                buffer.group_id, ref_seg.sample_name, buffer.segments.len() + 1
            );
        }

        // Compress reference using adaptive compression
        let (mut compressed, marker) = compress_reference_segment(&ref_seg.data)
            .context("Failed to compress reference")?;
        compressed.push(marker);

        // Metadata stores the uncompressed size
        let ref_size = ref_seg.data.len() as u64;

        // CRITICAL: Check if compression helped (matching C++ AGC segment.h lines 179, 204)
        // C++ AGC: if(packed_size + 1u < (uint32_t) data.size())
        // If compression didn't help, write UNCOMPRESSED raw data with metadata=0
        {
            let mut arch = archive.lock().unwrap();
            if compressed.len() < ref_seg.data.len() {
                // Compression helped - write compressed data with metadata=original_size
                arch.add_part(buffer.ref_stream_id, &compressed, ref_size)
                    .context("Failed to write reference")?;
            } else {
                // Compression didn't help - write UNCOMPRESSED data with metadata=0
                arch.add_part(buffer.ref_stream_id, &ref_seg.data, 0)
                    .context("Failed to write uncompressed reference")?;
            }
        }

        // Register reference in collection with in_group_id = 0
        {
            let mut coll = collection.lock().unwrap();
            coll.add_segment_placed(
                &ref_seg.sample_name,
                &ref_seg.contig_name,
                ref_seg.seg_part_no,
                buffer.group_id,
                0, // Reference is always at position 0
                ref_seg.is_rev_comp,
                ref_seg.data.len() as u32,
            )
            .context("Failed to register reference")?;
        }

        buffer.ref_written = true;
        buffer.reference_segment = Some(ref_seg.clone()); // Store for LZ encoding

        // Store reference in global map for split cost validation
        // This matches C++ AGC's v_segments which keeps all segment data in memory
        {
            let mut ref_segs = reference_segments.lock().unwrap();
            ref_segs.insert(buffer.group_id, ref_seg.data.clone());
        }

        // Prepare LZ encoder with reference (matching C++ AGC segment.cpp line 43: lz_diff->Prepare(s))
        // This is done ONCE when the reference is written, then reused for all subsequent segments
        if use_lz_encoding {
            let mut lz = LZDiff::new(config.min_match_len as u32);
            lz.prepare(&ref_seg.data);
            buffer.lz_diff = Some(lz);
        }
    }

    // CRITICAL: Sort segments by (sample_name, contig_name, seg_part_no) before processing
    // This matches C++ AGC's sort_known() in agc_compressor.h line 236-238
    // Without this sorting, segments are processed in file-system order which differs from
    // C++ AGC's lexicographic order, causing different in_group_id assignments and worse compression
    buffer.segments.sort();

    // Pack segments together with delta deduplication (matching C++ AGC segment.cpp lines 66-74)
    // Note: segments do NOT include the reference - it's stored separately
    //
    // CRITICAL FIX: Partial packs must persist across flush_pack calls to ensure pack boundaries
    // align with decompression expectations. Pack N must contain entries for in_group_ids
    // (N*50)+1 to (N+1)*50. Only write a pack when it has exactly 50 entries (or at finalization).
    // Use buffer.pending_deltas and buffer.pending_delta_ids to persist partial packs.

    let mut segment_in_group_ids: Vec<(usize, u32)> = Vec::new(); // (segment_index, in_group_id) for each segment

    // Helper closure to write a complete pack (exactly 50 entries)
    let write_complete_pack = |deltas: &[Vec<u8>], needs_raw_placeholder: bool, archive: &Arc<Mutex<Archive>>, stream_id: usize, config: &StreamingQueueConfig| -> Result<()> {
        if deltas.is_empty() && !needs_raw_placeholder {
            return Ok(());
        }

        let mut packed_data = Vec::new();

        // CRITICAL: Raw groups need a placeholder segment at position 0
        if needs_raw_placeholder {
            packed_data.push(0x7f);
            packed_data.push(CONTIG_SEPARATOR);
        }

        for delta in deltas.iter() {
            packed_data.extend_from_slice(delta);
            packed_data.push(CONTIG_SEPARATOR);
        }

        if packed_data.is_empty() {
            return Ok(());
        }

        let total_raw_size = packed_data.len();
        let mut compressed = compress_segment_configured(&packed_data, config.compression_level)
            .context("Failed to compress pack")?;
        compressed.push(0); // Marker 0 = plain ZSTD

        let mut arch = archive.lock().unwrap();
        if compressed.len() < total_raw_size {
            arch.add_part(stream_id, &compressed, total_raw_size as u64)
                .context("Failed to write compressed pack")?;
        } else {
            arch.add_part(stream_id, &packed_data, 0)
                .context("Failed to write uncompressed pack")?;
        }

        Ok(())
    };

    for (seg_idx, seg) in buffer.segments.iter().enumerate() {
        let contig_data = if !use_lz_encoding || buffer.reference_segment.is_none() {
            // Raw segment: groups 0-15 OR groups without reference
            seg.data.clone()
        } else {
            // LZ-encoded segment (groups >= 16 with reference)
            // DEBUG: Log sizes before encoding
            if let Some(ref_seg) = &buffer.reference_segment {
                if config.verbosity > 1 {
                    eprintln!("  LZ encoding: group={} ref_len={} target_len={} sample={} contig={} part={}",
                        buffer.group_id, ref_seg.data.len(), seg.data.len(),
                        seg.sample_name, seg.contig_name, seg.seg_part_no);
                }
            }
            // Reuse prepared lz_diff (matching C++ AGC segment.cpp line 59: lz_diff->Encode(s, delta))
            let ragc_encoded = buffer.lz_diff.as_mut()
                .expect("lz_diff should be prepared when reference is written")
                .encode(&seg.data);

            // Compare with C++ AGC encode
            #[cfg(feature = "cpp_agc")]
            if config.verbosity > 0 {
                if let Some(ref_seg) = &buffer.reference_segment {
                    if let Some(cpp_encoded) = crate::ragc_ffi::lzdiff_v2_encode(
                        &ref_seg.data,
                        &seg.data,
                        config.min_match_len as u32,
                    ) {
                        if ragc_encoded != cpp_encoded {
                            eprintln!("ENCODE_MISMATCH: ragc_len={} cpp_len={} ref_len={} seg_len={} group_id={}",
                                ragc_encoded.len(), cpp_encoded.len(), ref_seg.data.len(), seg.data.len(), buffer.group_id);
                            // Show context around first difference
                            for (i, (r, c)) in ragc_encoded.iter().zip(cpp_encoded.iter()).enumerate() {
                                if r != c {
                                    let start = if i > 10 { i - 10 } else { 0 };
                                    let ragc_ctx: Vec<_> = ragc_encoded[start..].iter().take(30).map(|&b| if b >= 32 && b < 127 { b as char } else { '?' }).collect();
                                    let cpp_ctx: Vec<_> = cpp_encoded[start..].iter().take(30).map(|&b| if b >= 32 && b < 127 { b as char } else { '?' }).collect();
                                    eprintln!("  first diff at byte {}: ragc='{}' cpp='{}'", i, *r as char, *c as char);
                                    eprintln!("  ragc context: {:?}", String::from_iter(ragc_ctx));
                                    eprintln!("  cpp  context: {:?}", String::from_iter(cpp_ctx));
                                    // Show the extra bytes in RAGC
                                    if ragc_encoded.len() > cpp_encoded.len() {
                                        let extra: Vec<_> = ragc_encoded[cpp_encoded.len()..].iter().take(20).map(|&b| if b >= 32 && b < 127 { b as char } else { '?' }).collect();
                                        eprintln!("  ragc extra bytes: {:?}", String::from_iter(extra));
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            ragc_encoded
        };

        // Handle LZ groups with IMPROVED_LZ_ENCODING: empty delta means same as reference
        // (matching C++ AGC segment.cpp lines 62-63)
        if use_lz_encoding && contig_data.is_empty() {
            // Same as reference - use in_group_id = 0
            segment_in_group_ids.push((seg_idx, 0));
            continue;
        }

        // Check if this delta already exists in pending pack (matching C++ AGC segment.cpp line 66)
        // Note: deduplication is per-pack, not global
        if let Some(existing_idx) = buffer.pending_deltas.iter().position(|d| d == &contig_data) {
            // Reuse existing delta's in_group_id (matching C++ AGC segment.cpp line 69)
            let reused_id = buffer.pending_delta_ids[existing_idx];
            segment_in_group_ids.push((seg_idx, reused_id));
        } else {
            // New unique delta - assign next in_group_id (matching C++ AGC segment.cpp lines 74, 77)
            // FIX: Apply .max(1) BEFORE using segments_written to ensure unique IDs when no reference
            // Bug was: max(0,1)=1, increment to 1 → max(1,1)=1 (COLLISION!)
            // Fixed: max(0,1)=1, id=1, increment to 2 → id=2 (UNIQUE!)
            buffer.segments_written = buffer.segments_written.max(1);
            let in_group_id = buffer.segments_written;
            buffer.segments_written += 1;
            buffer.pending_delta_ids.push(in_group_id);
            segment_in_group_ids.push((seg_idx, in_group_id));
            buffer.pending_deltas.push(contig_data);

            // CRITICAL: Flush when pending_deltas reaches PACK_CARDINALITY (matching C++ AGC segment.cpp lines 51-54)
            // C++ AGC: if (v_lzp.size() == contigs_in_pack) { store_in_archive(v_lzp); v_lzp.clear(); }
            if buffer.pending_deltas.len() == PACK_CARDINALITY {
                // Write complete pack with raw placeholder only if this is the first pack of a raw group
                let needs_placeholder = !use_lz_encoding && !buffer.raw_placeholder_written;
                write_complete_pack(&buffer.pending_deltas, needs_placeholder, archive, buffer.stream_id, config)?;
                buffer.raw_placeholder_written = true;

                // Clear for next pack - deduplication starts fresh
                buffer.pending_deltas.clear();
                buffer.pending_delta_ids.clear();
            }
        }
    }

    // DO NOT write partial pack here - leave it in buffer.pending_deltas for next flush_pack call
    // Partial packs are only written in finalize() to ensure pack boundaries align with decompression

    // Register ALL segments in collection with their in_group_ids
    for &(seg_idx, in_group_id) in segment_in_group_ids.iter() {
        let seg = &buffer.segments[seg_idx];
        let mut coll = collection.lock().unwrap();
        coll.add_segment_placed(
            &seg.sample_name,
            &seg.contig_name,
            seg.seg_part_no,
            buffer.group_id,
            in_group_id,
            seg.is_rev_comp,
            seg.data.len() as u32,
        )
        .context("Failed to register segment")?;
    }

    // Clear segments for next batch (but keep pending_deltas!)
    buffer.segments.clear();

    Ok(())
}

/// Write reference segment immediately when first segment arrives in group
/// (Matches C++ AGC segment.cpp lines 41-48: if (no_seqs == 0) writes reference right away)
/// This ensures LZ encoding works correctly for subsequent segments
fn write_reference_immediately(
    segment: &BufferedSegment,
    buffer: &mut SegmentGroupBuffer,
    collection: &Arc<Mutex<CollectionV3>>,
    archive: &Arc<Mutex<Archive>>,
    reference_segments: &Arc<Mutex<HashMap<u32, Vec<u8>>>>,
    reference_orientations: &Arc<Mutex<HashMap<u32, bool>>>,
    config: &StreamingQueueConfig,
) -> Result<()> {
    use crate::segment_compression::compress_reference_segment;

    if std::env::var("RAGC_DEBUG_REF_WRITE").is_ok() {
        eprintln!(
            "DEBUG_REF_IMMEDIATE: group={} sample={} contig={} seg={} data_len={}",
            buffer.group_id, segment.sample_name, segment.contig_name,
            segment.seg_part_no, segment.data.len()
        );
    }

    if config.verbosity > 1 {
        eprintln!(
            "  Writing immediate reference for group {}: {} {}:{} (part {})",
            buffer.group_id, segment.sample_name, segment.contig_name,
            segment.seg_part_no, segment.seg_part_no
        );
    }

    // 1. Compress reference using adaptive compression (matching flush_pack lines 635-637)
    let (mut compressed, marker) = compress_reference_segment(&segment.data)
        .context("Failed to compress reference")?;
    compressed.push(marker);

    let ref_size = segment.data.len();

    // 2. Write to archive immediately (matching C++ AGC segment.cpp line 43: store_in_archive)
    // CRITICAL: Check if compression helped (matching C++ AGC segment.h line 179)
    {
        let mut arch = archive.lock().unwrap();
        if compressed.len() < ref_size {
            // Compression helped - write compressed data with metadata=original_size
            arch.add_part(buffer.ref_stream_id, &compressed, ref_size as u64)
                .context("Failed to write compressed reference")?;
        } else {
            // Compression didn't help - write UNCOMPRESSED data with metadata=0
            arch.add_part(buffer.ref_stream_id, &segment.data, 0)
                .context("Failed to write uncompressed reference")?;
        }
    }

    // 3. Register reference in collection with in_group_id = 0 (matching flush_pack lines 650-661)
    {
        let mut coll = collection.lock().unwrap();
        coll.add_segment_placed(
            &segment.sample_name,
            &segment.contig_name,
            segment.seg_part_no,
            buffer.group_id,
            0, // Reference is always at position 0
            segment.is_rev_comp,
            segment.data.len() as u32,
        )
        .context("Failed to register immediate reference")?;
    }

    // 4. Mark reference as written and store for LZ encoding (matching flush_pack lines 663-664)
    buffer.ref_written = true;
    buffer.reference_segment = Some(segment.clone());
    // CRITICAL: Mark that in_group_id=0 is taken, so subsequent segments start from 1
    buffer.segments_written = 1;

    // 4b. Store reference data persistently (matching C++ AGC v_segments)
    // This enables LZ cost estimation for subsequent samples even after flush
    {
        let mut ref_segs = reference_segments.lock().unwrap();
        ref_segs.insert(buffer.group_id, segment.data.clone());
    }

    // 4c. Store reference orientation for ZERO_MATCH bug fix
    // When a delta segment joins this group later, it MUST use the same orientation
    // as the reference to ensure LZ encoding works correctly
    {
        let mut ref_orients = reference_orientations.lock().unwrap();
        ref_orients.insert(buffer.group_id, segment.is_rev_comp);
    }

    // 5. Prepare LZ encoder with reference (matching C++ AGC segment.cpp line 43: lz_diff->Prepare(s))
    // This is done ONCE when the reference is written, then reused for all subsequent segments
    let use_lz_encoding = buffer.group_id >= NO_RAW_GROUPS;
    if use_lz_encoding {
        let mut lz = LZDiff::new(config.min_match_len as u32);
        lz.prepare(&segment.data);
        buffer.lz_diff = Some(lz);
    }

    Ok(())
}

/// Compute reverse complement of a sequence
fn reverse_complement_sequence(seq: &[u8]) -> Vec<u8> {
    use crate::kmer::reverse_complement;
    seq.iter()
        .rev()
        .map(|&base| reverse_complement(base as u64) as u8)
        .collect()
}

/// Find best existing group for a segment with only one k-mer present
/// (Implements C++ AGC's find_cand_segment_with_one_splitter logic from lines 1659-1745)
fn find_group_with_one_kmer(
    kmer: u64,
    kmer_is_dir: bool,
    segment_data: &[u8],      // Segment data in forward orientation
    segment_data_rc: &[u8],   // Segment data in reverse complement
    map_segments_terminators: &Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    map_segments: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    reference_segments: &Arc<Mutex<HashMap<u32, Vec<u8>>>>,
    config: &StreamingQueueConfig,
) -> (u64, u64, bool) {
    let segment_len = segment_data.len();
    use crate::segment::MISSING_KMER;

    // Look up kmer in terminators map to find connected k-mers
    let connected_kmers = {
        let terminators = map_segments_terminators.lock().unwrap();
        match terminators.get(&kmer) {
            Some(vec) => vec.clone(),
            None => {
                // No connections found - create new group with MISSING
                // Match C++ AGC lines 1671-1679: check is_dir_oriented()
                // Debug: log entry to no-connection path
                if std::env::var("RAGC_DEBUG_IS_DIR").is_ok() {
                    eprintln!("RAGC_FIND_GROUP_NO_CONN: kmer={} kmer_is_dir={}", kmer, kmer_is_dir);
                }
                if kmer_is_dir {
                    // Dir-oriented: (kmer, MISSING) with rc=false
                    if config.verbosity > 1 {
                        #[cfg(feature = "verbose_debug")]
                        eprintln!("RAGC_CASE3_NO_CONNECTION: kmer={} is_dir=true -> ({}, MISSING) rc=false", kmer, kmer);
                    }
                    return (kmer, MISSING_KMER, false);
                } else {
                    // NOT dir-oriented: (MISSING, kmer) with rc=true
                    if config.verbosity > 1 {
                        #[cfg(feature = "verbose_debug")]
                        eprintln!("RAGC_CASE3_NO_CONNECTION: kmer={} is_dir=false -> (MISSING, {}) rc=true", kmer, kmer);
                    }
                    return (MISSING_KMER, kmer, true);
                }
            }
        }
    };

    if config.verbosity > 1 {
        #[cfg(feature = "verbose_debug")]
        eprintln!("RAGC_CASE3_FOUND_CONNECTIONS: kmer={} connections={}",
            kmer, connected_kmers.len());
    }
    // Debug: log connections found
    if std::env::var("RAGC_DEBUG_IS_DIR").is_ok() {
        eprintln!("RAGC_FIND_GROUP_FOUND_CONN: kmer={} kmer_is_dir={} connections={:?}",
            kmer, kmer_is_dir, connected_kmers);
    }

    // Build list of candidate groups
    // Each candidate: (key_front, key_back, needs_rc, ref_segment_size)
    let mut candidates: Vec<(u64, u64, bool, usize)> = Vec::new();

    {
        let groups = segment_groups.lock().unwrap();
        let seg_map = map_segments.lock().unwrap();
        let ref_segs = reference_segments.lock().unwrap();

        for &cand_kmer in &connected_kmers {
            // Create candidate group key normalized (smaller, larger)
            // C++ AGC lines 1691-1704
            //
            // IMPORTANT: When cand_kmer is MISSING, we need to try BOTH orderings!
            // Groups with MISSING k-mers can be stored as either (MISSING, kmer) or (kmer, MISSING)
            // depending on kmer_is_dir when they were created. We must match the actual stored key.
            let candidate_orderings: Vec<(u64, u64, bool)> = if cand_kmer == MISSING_KMER {
                // MISSING is involved - try both orderings to find the group
                vec![
                    (MISSING_KMER, kmer, true),   // (MISSING, kmer) with RC
                    (kmer, MISSING_KMER, false),  // (kmer, MISSING) without RC
                ]
            } else if cand_kmer < kmer {
                // cand_kmer is smaller - it goes first
                // This means we need to RC (C++ AGC line 1696: get<2>(ck) = true)
                vec![(cand_kmer, kmer, true)]
            } else {
                // kmer is smaller - it goes first
                // No RC needed (C++ AGC line 1703: get<2>(ck) = false)
                vec![(kmer, cand_kmer, false)]
            };

            for (key_front, key_back, needs_rc) in candidate_orderings {
                let cand_key = SegmentGroupKey {
                    kmer_front: key_front,
                    kmer_back: key_back,
                };

                // Check if this group exists in global registry OR batch-local buffer
                // CRITICAL FIX: Check BOTH seg_map (global) AND groups (batch-local)
                // This matches C++ AGC behavior where groups created in same batch are visible
                // Debug: trace candidate lookup
                if std::env::var("RAGC_DEBUG_IS_DIR").is_ok() {
                    eprintln!("RAGC_FIND_GROUP_CAND_CHECK: cand_key=({},{}) exists_in_seg_map={} exists_in_groups={}",
                        key_front, key_back, seg_map.contains_key(&cand_key), groups.contains_key(&cand_key));
                }

                // Check global registry first
                let in_global = seg_map.get(&cand_key);
                // Also check batch-local buffer
                let in_batch_local = groups.get(&cand_key);

                if in_global.is_some() || in_batch_local.is_some() {
                    // Get reference segment size from buffer OR persistent storage
                    let ref_size = if let Some(group_buffer) = in_batch_local {
                        // Group in buffer - get size from buffer
                        if let Some(ref_seg) = &group_buffer.reference_segment {
                            ref_seg.data.len()
                        } else {
                            segment_len // No reference yet, use current segment size
                        }
                    } else if let Some(&group_id) = in_global {
                        // Group flushed - get size from persistent storage
                        if let Some(ref_data) = ref_segs.get(&group_id) {
                            ref_data.len()
                        } else {
                            // Group exists but no reference data yet
                            segment_len
                        }
                    } else {
                        segment_len
                    };

                    candidates.push((key_front, key_back, needs_rc, ref_size));
                    break; // Found a match for this cand_kmer, move to next
                }
            }
        }
    }

    if candidates.is_empty() {
        // No existing groups found - create new with MISSING
        // Must match C++ AGC is_dir_oriented logic (same as no-connections case above)
        if std::env::var("RAGC_DEBUG_IS_DIR").is_ok() {
            if kmer_is_dir {
                eprintln!("RAGC_FIND_GROUP_NO_CAND: kmer={} kmer_is_dir={} -> returning ({},MISSING,false)",
                    kmer, kmer_is_dir, kmer);
            } else {
                eprintln!("RAGC_FIND_GROUP_NO_CAND: kmer={} kmer_is_dir={} -> returning (MISSING,{},true)",
                    kmer, kmer_is_dir, kmer);
            }
        }
        if kmer_is_dir {
            // Dir-oriented: (kmer, MISSING) with rc=false
            if config.verbosity > 1 {
                #[cfg(feature = "verbose_debug")]
                eprintln!("RAGC_CASE3_NO_CANDIDATES: kmer={} is_dir=true -> ({}, MISSING) rc=false", kmer, kmer);
            }
            return (kmer, MISSING_KMER, false);
        } else {
            // NOT dir-oriented: (MISSING, kmer) with rc=true
            if config.verbosity > 1 {
                #[cfg(feature = "verbose_debug")]
                eprintln!("RAGC_CASE3_NO_CANDIDATES: kmer={} is_dir=false -> (MISSING, {}) rc=true", kmer, kmer);
            }
            return (MISSING_KMER, kmer, true);
        }
    }

    // Sort candidates by reference segment size (C++ AGC lines 1710-1719)
    // Prefer candidates with ref size closest to our segment size
    candidates.sort_by(|a, b| {
        let a_diff = (a.3 as i64 - segment_len as i64).abs();
        let b_diff = (b.3 as i64 - segment_len as i64).abs();

        if a_diff != b_diff {
            a_diff.cmp(&b_diff)
        } else {
            a.3.cmp(&b.3) // If equal distance, prefer smaller ref size
        }
    });

    // Debug: Print sorted candidates before evaluation
    if config.verbosity > 2 {
        eprintln!(
            "RAGC_CASE3_SORTED_CANDIDATES: kmer={} segment_len={} n_candidates={}",
            kmer, segment_len, candidates.len()
        );
        for (i, &(kf, kb, rc, rs)) in candidates.iter().enumerate() {
            let size_diff = (rs as i64 - segment_len as i64).abs();
            eprintln!(
                "  CAND[{}]: ({},{}) rc={} ref_size={} size_diff={}",
                i, kf, kb, rc, rs, size_diff
            );
        }
    }

    // Test compression for each candidate (C++ AGC lines 1726-1788)
    // Match C++ AGC's TWO-PASS approach:
    //   Pass 1: Compute all estimates, track minimum (lines 1726-1732)
    //   Pass 2: Pick candidate with minimum estimate (lines 1775-1787)
    //
    // CRITICAL: Initialize best_pk to (~0ull, ~0ull) like C++ AGC (line 1628)
    let mut best_key_front = u64::MAX;  // ~0ull in C++
    let mut best_key_back = u64::MAX;   // ~0ull in C++
    let mut best_needs_rc = false;
    let mut best_estim_size = if segment_len < 16 {
        segment_len
    } else {
        segment_len - 16
    };

    // Pass 1: Compute estimates and find minimum
    // Store estimates alongside candidates: Vec<(front, back, needs_rc, ref_size, estim_size)>
    let mut candidate_estimates: Vec<(u64, u64, bool, usize, usize)> = Vec::new();

    {
        let groups = segment_groups.lock().unwrap();
        let seg_map = map_segments.lock().unwrap();
        let ref_segs = reference_segments.lock().unwrap();

        for &(key_front, key_back, needs_rc, ref_size) in &candidates {
            let cand_key = SegmentGroupKey {
                kmer_front: key_front,
                kmer_back: key_back,
            };

            // Get the reference segment for this candidate from buffer OR persistent storage
            let (ref_data_opt, ref_source): (Option<&[u8]>, &str) = if let Some(group_buffer) = groups.get(&cand_key) {
                if config.verbosity > 2 && key_front == 1244212049458757632 && key_back == 1244212049458757632 {
                    let ref_seg = group_buffer.reference_segment.as_ref();
                    let ref_len = ref_seg.map(|s| s.data.len()).unwrap_or(0);
                    let ref_first5: Vec<u8> = ref_seg.map(|s| s.data.iter().take(5).cloned().collect()).unwrap_or_default();
                    eprintln!("RAGC_REF_LOOKUP_BUFFER: degenerate key ({},{}) buffer ref_len={} ref[0..5]={:?}",
                        key_front, key_back, ref_len, ref_first5);
                }
                (group_buffer.reference_segment.as_ref().map(|seg| seg.data.as_slice()), "buffer")
            } else if let Some(&group_id) = seg_map.get(&cand_key) {
                if config.verbosity > 2 && key_front == 1244212049458757632 && key_back == 1244212049458757632 {
                    let ref_data = ref_segs.get(&group_id);
                    let ref_len = ref_data.map(|d| d.len()).unwrap_or(0);
                    let ref_first5: Vec<u8> = ref_data.map(|d| d.iter().take(5).cloned().collect()).unwrap_or_default();
                    eprintln!("RAGC_REF_LOOKUP_PERSISTENT: degenerate key ({},{}) -> group_id={} ref_len={} ref[0..5]={:?}",
                        key_front, key_back, group_id, ref_len, ref_first5);
                }
                (ref_segs.get(&group_id).map(|data| data.as_slice()), "persistent")
            } else {
                (None, "none")
            };
            let ref_data_opt = ref_data_opt;

            if let Some(ref_data) = ref_data_opt {
                // Test LZ encoding against this reference (C++ AGC line 1728: estimate())
                let target_data = if needs_rc {
                    segment_data_rc
                } else {
                    segment_data
                };

                // Compute estimate - compare both RAGC native and C++ FFI when verbose
                let estim_size = {
                    let mut lz = LZDiff::new(config.min_match_len as u32);
                    lz.prepare(&ref_data.to_vec());
                    // Use estimate() which matches C++ CLZDiff_V2::Estimate exactly
                    lz.estimate(&target_data.to_vec(), best_estim_size as u32) as usize
                };

                // Also compute with C++ FFI and compare
                #[cfg(feature = "cpp_agc")]
                let cpp_estim_size = crate::ragc_ffi::lzdiff_v2_estimate(
                    ref_data,
                    target_data,
                    config.min_match_len as u32,
                    best_estim_size as u32,
                ) as usize;

                #[cfg(feature = "cpp_agc")]
                if estim_size != cpp_estim_size && config.verbosity > 0 {
                    eprintln!("ESTIMATE_MISMATCH: ragc={} cpp={} ref_len={} tgt_len={} bound={}",
                        estim_size, cpp_estim_size, ref_data.len(), target_data.len(), best_estim_size);
                }

                // DEBUG: Also compute estimate with initial threshold to check if tie would occur
                #[cfg(not(feature = "cpp_agc"))]
                let estim_no_bound = if config.verbosity > 2 {
                    let mut lz2 = LZDiff::new(config.min_match_len as u32);
                    lz2.prepare(&ref_data.to_vec());
                    lz2.estimate(&target_data.to_vec(), (segment_len - 16) as u32) as usize
                } else { 0 };

                if config.verbosity > 2 {
                    // Print detailed debug info including bound and first/last bytes
                    let ref_first: Vec<u8> = ref_data.iter().take(5).cloned().collect();
                    let ref_last: Vec<u8> = ref_data.iter().rev().take(5).cloned().collect();
                    let tgt_first: Vec<u8> = target_data.iter().take(5).cloned().collect();
                    let tgt_last: Vec<u8> = target_data.iter().rev().take(5).cloned().collect();
                    #[cfg(not(feature = "cpp_agc"))]
                    eprintln!(
                        "RAGC_CASE3_ESTIMATE: kmer={} cand=({},{}) rc={} ref_len={} target_len={} bound={} estim={} estim_nobound={} ref[0..5]={:?} ref[-5..]={:?} tgt[0..5]={:?} tgt[-5..]={:?}",
                        kmer, key_front, key_back, needs_rc, ref_data.len(), target_data.len(), best_estim_size, estim_size, estim_no_bound, ref_first, ref_last, tgt_first, tgt_last
                    );
                    #[cfg(feature = "cpp_agc")]
                    eprintln!(
                        "RAGC_CASE3_ESTIMATE: kmer={} cand=({},{}) rc={} ref_len={} target_len={} bound={} estim={} ref[0..5]={:?} ref[-5..]={:?} tgt[0..5]={:?} tgt[-5..]={:?}",
                        kmer, key_front, key_back, needs_rc, ref_data.len(), target_data.len(), best_estim_size, estim_size, ref_first, ref_last, tgt_first, tgt_last
                    );
                }

                // Track minimum estim_size (C++ AGC lines 1730-1732)
                if estim_size < best_estim_size {
                    best_estim_size = estim_size;
                }

                candidate_estimates.push((key_front, key_back, needs_rc, ref_size, estim_size));
            }
        }
    }

    // Pass 2: Pick candidate with minimum estimate among ALL candidates, using tie-breakers
    // (C++ AGC lines 1775-1788)
    //
    // CRITICAL FIX: C++ AGC only picks candidates that BEAT the initial threshold (segment_size - 16).
    // If no candidate beats the threshold, best_pk stays at (~0ull, ~0ull) and fallback MISSING is used.
    //
    // The previous bug was unconditionally picking the first candidate (first_candidate = true).
    // This caused RAGC to always pick the first candidate even when its estimate was worse than threshold,
    // preventing fallback to existing MISSING groups.
    //
    // C++ AGC's selection logic (lines 1780-1787):
    //   if (v_estim_size[i] < best_estim_size || ...)
    // This only updates best_pk if estimate is BETTER than current best (initially threshold).
    for &(key_front, key_back, needs_rc, _ref_size, estim_size) in &candidate_estimates {
        let cand_pk = (key_front, key_back);
        let best_pk = (best_key_front, best_key_back);

        // Match C++ AGC's selection logic exactly (lines 1780-1787):
        // Only pick candidate if:
        // - Smaller estimate than current best (initially threshold), OR
        // - Same estimate with lexicographically smaller pk, OR
        // - Same estimate+pk with better RC (prefers forward orientation)
        if estim_size < best_estim_size
            || (estim_size == best_estim_size && cand_pk < best_pk)
            || (estim_size == best_estim_size && cand_pk == best_pk && !needs_rc)
        {
            best_estim_size = estim_size;
            best_key_front = key_front;
            best_key_back = key_back;
            best_needs_rc = needs_rc;
        }
    }

    // Debug: Print Pass 2 results
    if config.verbosity > 2 && !candidate_estimates.is_empty() {
        let threshold = if segment_len < 16 { segment_len } else { segment_len - 16 };
        eprintln!(
            "RAGC_CASE3_PASS2_RESULTS: threshold={} best=({},{}) best_estim={}",
            threshold, best_key_front, best_key_back, best_estim_size
        );
        for (i, &(kf, kb, rc, rs, es)) in candidate_estimates.iter().enumerate() {
            let is_winner = kf == best_key_front && kb == best_key_back;
            let marker = if is_winner { "*WINNER*" } else { "" };
            eprintln!(
                "  RESULT[{}]: ({},{}) rc={} ref_size={} estimate={} {}",
                i, kf, kb, rc, rs, es, marker
            );
        }
    }

    // If no candidate was selected (best_pk is still (~0ull, ~0ull)), create MISSING key
    // This matches C++ AGC lines 1791-1799: fallback to (kmer, MISSING) or (MISSING, kmer)
    if best_key_front == u64::MAX && best_key_back == u64::MAX {
        if kmer_is_dir {
            // Dir-oriented: (kmer, MISSING) with rc=false
            if config.verbosity > 1 {
                #[cfg(feature = "verbose_debug")]
                eprintln!("RAGC_CASE3_NO_WINNER: kmer={} is_dir=true -> ({}, MISSING) rc=false", kmer, kmer);
            }
            return (kmer, MISSING_KMER, false);
        } else {
            // NOT dir-oriented: (MISSING, kmer) with rc=true
            if config.verbosity > 1 {
                #[cfg(feature = "verbose_debug")]
                eprintln!("RAGC_CASE3_NO_WINNER: kmer={} is_dir=false -> (MISSING, {}) rc=true", kmer, kmer);
            }
            return (MISSING_KMER, kmer, true);
        }
    }

    if config.verbosity > 1 {
        #[cfg(feature = "verbose_debug")]
        eprintln!("RAGC_CASE3_PICKED: kmer={} best=({},{}) rc={} estim_size={} segment_size={}",
            kmer, best_key_front, best_key_back, best_needs_rc, best_estim_size, segment_len);
    }

    (best_key_front, best_key_back, best_needs_rc)
}

/// Find candidate segment using fallback minimizers
/// Matches C++ AGC's find_cand_segment_using_fallback_minimizers (lines 1807-1958)
///
/// This function is called when Case 3 (one k-mer present) fails to find a good match.
/// It scans the segment for k-mers that pass the fallback filter, looks them up in
/// the fallback minimizers map, and finds candidate groups with shared k-mers.
///
/// # Arguments
/// * `segment_data` - The segment data to search
/// * `k` - K-mer length
/// * `min_shared_kmers` - Minimum number of shared k-mers to consider a candidate
/// * `fallback_filter` - Filter to select which k-mers to check
/// * `map_fallback_minimizers` - Map from k-mer to candidate group keys
/// * `map_segments` - Map from group key to group ID
/// * `segment_groups` - Buffer of segment groups
/// * `reference_segments` - Stored reference segments
/// * `config` - Compression configuration
///
/// # Returns
/// (key_front, key_back, should_reverse) if a candidate is found, or (MISSING, MISSING, false) if none
#[allow(clippy::too_many_arguments)]
fn find_cand_segment_using_fallback_minimizers(
    segment_data: &[u8],
    segment_data_rc: &[u8],
    k: usize,
    min_shared_kmers: u64,
    fallback_filter: &FallbackFilter,
    map_fallback_minimizers: &Arc<Mutex<HashMap<u64, Vec<(u64, u64)>>>>,
    map_segments: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    reference_segments: &Arc<Mutex<HashMap<u32, Vec<u8>>>>,
    config: &StreamingQueueConfig,
) -> (u64, u64, bool) {
    use crate::segment::MISSING_KMER;

    const MAX_NUM_TO_ESTIMATE: usize = 10;
    let short_segments = config.segment_size <= 10000;
    let segment_len = segment_data.len();

    if !fallback_filter.is_enabled() {
        return (MISSING_KMER, MISSING_KMER, false);
    }

    // Scan segment for k-mers and count candidates
    // Map from candidate group key to list of shared k-mers
    let mut cand_seg_counts: HashMap<(u64, u64), Vec<u64>> = HashMap::new();

    // K-mer scanning state (matches C++ AGC CKmer behavior)
    let mut kmer_data: u64 = 0;
    let mut kmer_rc: u64 = 0;
    let mut kmer_len: usize = 0;
    let mask: u64 = (1u64 << (2 * k)) - 1;

    // Scan segment for k-mers
    for &base in segment_data {
        if base > 3 {
            // Non-ACGT character - reset k-mer
            kmer_data = 0;
            kmer_rc = 0;
            kmer_len = 0;
            continue;
        }

        // Add base to forward k-mer (shift left, add at LSB)
        kmer_data = ((kmer_data << 2) | (base as u64)) & mask;

        // Add complement to reverse k-mer (shift right, add at MSB)
        let comp = 3 - base; // A<->T, C<->G
        kmer_rc = (kmer_rc >> 2) | ((comp as u64) << (2 * (k - 1)));

        kmer_len += 1;

        if kmer_len >= k {
            // Use canonical k-mer (smaller of forward and reverse)
            let canonical = kmer_data.min(kmer_rc);
            let is_dir_oriented = kmer_data <= kmer_rc;

            // Check if k-mer passes fallback filter and is not symmetric
            if fallback_filter.passes(canonical) && kmer_data != kmer_rc {
                // Look up in fallback minimizers map
                let fb_map = map_fallback_minimizers.lock().unwrap();
                if let Some(candidates) = fb_map.get(&canonical) {
                    for &(key1, key2) in candidates {
                        // Skip MISSING keys
                        if key1 == MISSING_KMER || key2 == MISSING_KMER {
                            continue;
                        }

                        // Normalize based on orientation
                        let cand_key = if !is_dir_oriented {
                            (key2, key1)
                        } else {
                            (key1, key2)
                        };

                        cand_seg_counts.entry(cand_key)
                            .or_insert_with(Vec::new)
                            .push(canonical);
                    }
                }
            }
        }
    }

    // Prune candidates to those with >= min_shared_kmers unique k-mers
    let mut pruned_candidates: Vec<(u64, (u64, u64))> = Vec::new();
    for (key, mut kmers) in cand_seg_counts {
        kmers.sort_unstable();
        kmers.dedup();
        let unique_count = kmers.len() as u64;
        if unique_count >= min_shared_kmers {
            pruned_candidates.push((unique_count, key));
        }
    }

    if pruned_candidates.is_empty() {
        if config.verbosity > 1 {
            #[cfg(feature = "verbose_debug")]
            eprintln!("RAGC_FALLBACK_NO_CANDIDATES: min_shared={}", min_shared_kmers);
        }
        return (MISSING_KMER, MISSING_KMER, false);
    }

    // Sort by count (descending) and take top MAX_NUM_TO_ESTIMATE
    pruned_candidates.sort_by(|a, b| b.0.cmp(&a.0));
    if pruned_candidates.len() > MAX_NUM_TO_ESTIMATE {
        pruned_candidates.truncate(MAX_NUM_TO_ESTIMATE);
    }

    // Avoid trying poor candidates (less than half the best count)
    let best_count = pruned_candidates[0].0;
    pruned_candidates.retain(|c| c.0 * 2 >= best_count);

    if config.verbosity > 1 {
        #[cfg(feature = "verbose_debug")]
        eprintln!("RAGC_FALLBACK_CANDIDATES: count={} best_shared={} min_shared={}",
            pruned_candidates.len(), best_count, min_shared_kmers);
    }

    // For short segments, use fast decision based on shared k-mer count
    if short_segments {
        let (count, (key_front, key_back)) = pruned_candidates[0];
        if config.verbosity > 1 {
            #[cfg(feature = "verbose_debug")]
            eprintln!("RAGC_FALLBACK_SHORT_SEGMENT: key=({},{}) shared_kmers={}", key_front, key_back, count);
        }
        // Normalize: ensure front <= back
        if key_front <= key_back {
            return (key_front, key_back, false);
        } else {
            return (key_back, key_front, true);
        }
    }

    // For longer segments, estimate compression cost for each candidate
    let mut best_key: Option<(u64, u64)> = None;
    let mut best_estimate: usize = segment_len;
    let mut _best_is_rc = false;

    {
        let groups = segment_groups.lock().unwrap();
        let seg_map = map_segments.lock().unwrap();
        let ref_segs = reference_segments.lock().unwrap();

        for &(_count, (key_front, key_back)) in &pruned_candidates {
            // Normalize key
            let (norm_front, norm_back, is_seg_rc) = if key_front <= key_back {
                (key_front, key_back, false)
            } else {
                (key_back, key_front, true)
            };

            let cand_key = SegmentGroupKey {
                kmer_front: norm_front,
                kmer_back: norm_back,
            };

            // Get reference segment for this candidate
            let ref_data_opt: Option<&[u8]> = if let Some(group_buffer) = groups.get(&cand_key) {
                group_buffer.reference_segment.as_ref().map(|seg| seg.data.as_slice())
            } else if let Some(&group_id) = seg_map.get(&cand_key) {
                ref_segs.get(&group_id).map(|data| data.as_slice())
            } else {
                None
            };

            if let Some(ref_data) = ref_data_opt {
                let target_data = if is_seg_rc { segment_data_rc } else { segment_data };

                // Estimate compression cost
                #[cfg(feature = "cpp_agc")]
                let estimate = crate::ragc_ffi::lzdiff_v2_estimate(
                    ref_data,
                    target_data,
                    config.min_match_len as u32,
                    best_estimate as u32,
                ) as usize;

                #[cfg(not(feature = "cpp_agc"))]
                let estimate = {
                    let mut lz = LZDiff::new(config.min_match_len as u32);
                    lz.prepare(&ref_data.to_vec());
                    // Use estimate() which matches C++ CLZDiff_V2::Estimate exactly
                    lz.estimate(&target_data.to_vec(), best_estimate as u32) as usize
                };

                if config.verbosity > 2 {
                    #[cfg(feature = "verbose_debug")]
                    eprintln!("RAGC_FALLBACK_ESTIMATE: key=({},{}) rc={} estimate={}",
                        norm_front, norm_back, is_seg_rc, estimate);
                }

                // Track best (lowest estimate)
                if estimate > 0 && estimate < best_estimate {
                    best_estimate = estimate;
                    best_key = Some((norm_front, norm_back));
                    _best_is_rc = is_seg_rc;
                }
            }
        }
    }

    // In adaptive mode, check if result is worth using
    if config.adaptive_mode {
        let threshold = if short_segments {
            (segment_len as f64 * 0.9) as usize
        } else {
            (segment_len as f64 * 0.2) as usize
        };

        if best_estimate >= threshold {
            if config.verbosity > 1 {
                #[cfg(feature = "verbose_debug")]
                eprintln!("RAGC_FALLBACK_ADAPTIVE_REJECT: estimate={} threshold={}", best_estimate, threshold);
            }
            return (MISSING_KMER, MISSING_KMER, false);
        }
    }

    match best_key {
        Some((front, back)) => {
            // Normalize: ensure front <= back
            if front <= back {
                if config.verbosity > 1 {
                    #[cfg(feature = "verbose_debug")]
                    eprintln!("RAGC_FALLBACK_PICKED: key=({},{}) rc=false estimate={}", front, back, best_estimate);
                }
                (front, back, false)
            } else {
                if config.verbosity > 1 {
                    #[cfg(feature = "verbose_debug")]
                    eprintln!("RAGC_FALLBACK_PICKED: key=({},{}) rc=true estimate={}", back, front, best_estimate);
                }
                (back, front, true)
            }
        }
        None => {
            if config.verbosity > 1 {
                #[cfg(feature = "verbose_debug")]
                eprintln!("RAGC_FALLBACK_NO_WINNER: no candidate beat threshold");
            }
            (MISSING_KMER, MISSING_KMER, false)
        }
    }
}

/// Add fallback mapping for a segment's k-mers
/// Matches C++ AGC's add_fallback_mapping (lines 1961-1989)
///
/// Called when a segment is assigned to a group to populate the fallback minimizers map.
fn add_fallback_mapping(
    segment_data: &[u8],
    k: usize,
    splitter1: u64,
    splitter2: u64,
    fallback_filter: &FallbackFilter,
    map_fallback_minimizers: &Arc<Mutex<HashMap<u64, Vec<(u64, u64)>>>>,
) {
    use crate::segment::MISSING_KMER;

    if !fallback_filter.is_enabled() {
        return;
    }

    // Skip if splitters are MISSING
    if splitter1 == MISSING_KMER || splitter2 == MISSING_KMER {
        return;
    }

    let splitter_dir = (splitter1, splitter2);
    let splitter_rev = (splitter2, splitter1);
    let mask: u64 = (1u64 << (2 * k)) - 1;

    // K-mer scanning state
    let mut kmer_data: u64 = 0;
    let mut kmer_rc: u64 = 0;
    let mut kmer_len: usize = 0;

    let mut fb_map = map_fallback_minimizers.lock().unwrap();

    for &base in segment_data {
        if base > 3 {
            kmer_data = 0;
            kmer_rc = 0;
            kmer_len = 0;
            continue;
        }

        kmer_data = ((kmer_data << 2) | (base as u64)) & mask;
        let comp = 3 - base;
        kmer_rc = (kmer_rc >> 2) | ((comp as u64) << (2 * (k - 1)));
        kmer_len += 1;

        if kmer_len >= k {
            let canonical = kmer_data.min(kmer_rc);
            let is_dir_oriented = kmer_data <= kmer_rc;

            // Check filter and skip symmetric k-mers
            if fallback_filter.passes(canonical) && kmer_data != kmer_rc {
                let to_add = if is_dir_oriented { splitter_dir } else { splitter_rev };
                let entry = fb_map.entry(canonical).or_insert_with(Vec::new);

                // Only add if not already present
                if !entry.contains(&to_add) {
                    entry.push(to_add);
                }
            }
        }
    }
}

/// Flush batch-local groups to global state (matches C++ AGC batch boundary)
/// This updates the global map_segments registry with batch-local groups,
/// then clears the batch-local state (like C++ AGC destroying m_kmers at batch end)
fn flush_batch(
    _segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    _pending_batch_segments: &Arc<Mutex<Vec<PendingSegment>>>,
    batch_local_groups: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    batch_local_terminators: &Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    map_segments: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    _group_counter: &Arc<AtomicU32>,
    map_segments_terminators: &Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    _archive: &Arc<Mutex<Archive>>,
    _collection: &Arc<Mutex<CollectionV3>>,
    _reference_segments: &Arc<Mutex<HashMap<u32, Vec<u8>>>>,
    #[cfg(feature = "cpp_agc")]
    _grouping_engine: &Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,
    config: &StreamingQueueConfig,
) -> Result<()> {
    // Get batch-local groups
    let batch_map = batch_local_groups.lock().unwrap();
    let batch_terms = batch_local_terminators.lock().unwrap();

    if batch_map.is_empty() && batch_terms.is_empty() {
        #[cfg(feature = "verbose_debug")]
        if config.verbosity > 0 {
            eprintln!("FLUSH_BATCH: No pending groups to flush");
        }
        return Ok(());
    }

    #[cfg(feature = "verbose_debug")]
    if config.verbosity > 0 {
        eprintln!("FLUSH_BATCH: Updating global registry with {} batch-local groups, {} terminator keys",
            batch_map.len(), batch_terms.len());
    }

    // Update global registry with this batch's new groups
    {
        let mut global_map = map_segments.lock().unwrap();
        for (key, group_id) in batch_map.iter() {
            // Only insert if not already present (shouldn't happen, but defensive)
            global_map.entry(key.clone()).or_insert(*group_id);
        }
    }

    // CRITICAL: Merge batch-local terminators into global terminators
    // This is where C++ AGC makes terminators visible for find_middle in subsequent samples
    {
        let mut global_terms = map_segments_terminators.lock().unwrap();
        for (kmer, connections) in batch_terms.iter() {
            let entry = global_terms.entry(*kmer).or_insert_with(Vec::new);
            entry.extend(connections.iter().cloned());
            entry.sort_unstable();
            entry.dedup();
        }
    }

    // Clear batch-local state (like C++ AGC destroying m_kmers)
    drop(batch_map);
    drop(batch_terms);
    batch_local_groups.lock().unwrap().clear();
    batch_local_terminators.lock().unwrap().clear();

    #[cfg(feature = "verbose_debug")]
    if config.verbosity > 0 {
        eprintln!("FLUSH_BATCH: Batch flush complete, batch-local state cleared");
    }

    Ok(())
}

/// Helper function to fix orientation of segment data to match reference orientation.
/// Returns (fixed_data, fixed_is_rev_comp) tuple.
/// Used for both normal segments and split segments to ensure consistent orientation within groups.
fn fix_orientation_for_group(
    data: &[u8],
    should_reverse: bool,
    key: &SegmentGroupKey,
    map_segments: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    batch_local_groups: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    reference_orientations: &Arc<Mutex<HashMap<u32, bool>>>,
) -> (Vec<u8>, bool) {
    // Look up group_id in both global (map_segments) and batch-local registries
    let group_id_opt = {
        let seg_map = map_segments.lock().unwrap();
        if let Some(&gid) = seg_map.get(key) {
            Some(gid)
        } else {
            drop(seg_map);
            let batch_map = batch_local_groups.lock().unwrap();
            batch_map.get(key).copied()
        }
    };

    if let Some(group_id) = group_id_opt {
        let ref_orients = reference_orientations.lock().unwrap();
        if let Some(&ref_rc) = ref_orients.get(&group_id) {
            // Group has a stored reference orientation
            if ref_rc != should_reverse {
                // Orientation mismatch - transform data to match reference
                // The data is already transformed based on should_reverse, so we need to
                // apply RC if ref_rc is true (to get RC), or use original if ref_rc is false (to get forward)
                // But we're given the ALREADY-TRANSFORMED data, so we need to invert if mismatch
                let fixed_data: Vec<u8> = data.iter().rev().map(|&base| {
                    match base {
                        0 => 3, // A -> T
                        1 => 2, // C -> G
                        2 => 1, // G -> C
                        3 => 0, // T -> A
                        _ => base, // N or other non-ACGT
                    }
                }).collect();
                return (fixed_data, ref_rc);
            }
        }
    }
    // No fix needed - return as-is
    (data.to_vec(), should_reverse)
}

/// Worker thread that pulls from queue and compresses
fn worker_thread(
    worker_id: usize,
    queue: Arc<MemoryBoundedQueue<ContigTask>>,
    collection: Arc<Mutex<CollectionV3>>,
    splitters: Arc<HashSet<u64>>,
    archive: Arc<Mutex<Archive>>,
    segment_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    group_counter: Arc<AtomicU32>,
    raw_group_counter: Arc<AtomicU32>,
    reference_sample_name: Arc<Mutex<Option<String>>>,
    map_segments: Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    reference_segments: Arc<Mutex<HashMap<u32, Vec<u8>>>>,
    reference_orientations: Arc<Mutex<HashMap<u32, bool>>>,
    split_offsets: Arc<Mutex<HashMap<(String, String, usize), usize>>>,
    #[cfg(feature = "cpp_agc")]
    grouping_engine: Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,
    current_batch_sample: Arc<Mutex<Option<String>>>,
    batch_local_groups: Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    batch_local_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    pending_batch_segments: Arc<Mutex<Vec<PendingSegment>>>,
    map_fallback_minimizers: Arc<Mutex<HashMap<u64, Vec<(u64, u64)>>>>,
    config: StreamingQueueConfig,
) -> Result<()> {
    let mut processed_count = 0;

    // Create fallback filter from config
    let fallback_filter = FallbackFilter::new(config.fallback_frac);

    loop {
        // Pull from queue (blocks if empty, returns None when closed)
        let Some(task) = queue.pull() else {
            // Queue is closed and empty - flush any pending batch before exiting
            if config.verbosity > 0 {
                eprintln!("Worker {} flushing final batch before exit", worker_id);
            }
            flush_batch(
                &segment_groups,
                &pending_batch_segments,
                &batch_local_groups,
                &batch_local_terminators,
                &map_segments,
                &group_counter,
                &map_segments_terminators,
                &archive,
                &collection,
                &reference_segments,
                #[cfg(feature = "cpp_agc")]
                &grouping_engine,
                &config,
            ).ok(); // Ignore errors on final flush

            if config.verbosity > 1 {
                eprintln!(
                    "Worker {} finished ({} contigs processed)",
                    worker_id, processed_count
                );
            }
            break;
        };

        // Check if we're starting a new sample (batch boundary)
        {
            let mut current = current_batch_sample.lock().unwrap();
            if current.as_ref() != Some(&task.sample_name) {
                if current.is_some() {
                    // Batch boundary detected - flush previous batch
                    #[cfg(feature = "verbose_debug")]
                    if config.verbosity > 0 {
                        eprintln!("BATCH_BOUNDARY: Flushing batch for sample {:?}, starting new batch for {}",
                            current, task.sample_name);
                    }
                    drop(current); // Release lock before flush
                    flush_batch(
                        &segment_groups,
                        &pending_batch_segments,
                        &batch_local_groups,
                        &batch_local_terminators,
                        &map_segments,
                        &group_counter,
                        &map_segments_terminators,
                        &archive,
                        &collection,
                        &reference_segments,
                        #[cfg(feature = "cpp_agc")]
                        &grouping_engine,
                        &config,
                    )?;
                    current = current_batch_sample.lock().unwrap();
                }
                *current = Some(task.sample_name.clone());
            }
        }

        // Split into segments
        let segments =
            split_at_splitters_with_size(&task.data, &splitters, config.k, config.segment_size);

        if config.verbosity > 2 {
            eprintln!(
                "Worker {} processing {} (split into {} segments)",
                worker_id,
                task.contig_name,
                segments.len()
            );
        }

        // Buffer segments for packing (matching batch mode)
        for (original_place, segment) in segments.iter().enumerate() {
            // Calculate adjusted place based on prior splits in this contig
            // (matches C++ AGC lines 2033-2036: increment seg_part_no twice when split occurs)
            let place = {
                let offsets = split_offsets.lock().unwrap();
                let mut adjusted = original_place;
                // Count how many splits occurred before this position
                for pos in 0..original_place {
                    if offsets.contains_key(&(task.sample_name.clone(), task.contig_name.clone(), pos)) {
                        adjusted += 1;
                    }
                }
                adjusted
            };

            // DEBUG: Output every segment for comparison with C++ AGC
            #[cfg(feature = "verbose_debug")]
            eprintln!("RAGC_SEGMENT: sample={} contig={} part={} len={} front={} back={}",
                task.sample_name, task.contig_name, place, segment.data.len(),
                segment.front_kmer, segment.back_kmer);

            // Match C++ AGC Case 2: Normalize segment group key by ensuring front <= back
            // (agc_compressor.cpp lines 1306-1327)
            use crate::segment::MISSING_KMER;

            // Precompute reverse complement for all cases that might need it
            // Segment data uses numeric encoding: 0=A, 1=C, 2=G, 3=T
            let segment_data_rc: Vec<u8> = segment.data.iter().rev().map(|&base| {
                match base {
                    0 => 3, // A -> T
                    1 => 2, // C -> G
                    2 => 1, // G -> C
                    3 => 0, // T -> A
                    _ => base, // N or other non-ACGT
                }
            }).collect();

            let (key_front, key_back, should_reverse) =
                if segment.front_kmer != MISSING_KMER && segment.back_kmer != MISSING_KMER {
                    // Both k-mers present
                    // C++ AGC uses `<` not `<=`, which means degenerate k-mers (front == back)
                    // go to the else branch and get store_rc=true (lines 1306-1313)
                    if segment.front_kmer < segment.back_kmer {
                        // Already normalized - keep original orientation
                        if config.verbosity > 2 {
                            #[cfg(feature = "verbose_debug")]
                            eprintln!(
                                "RAGC_CASE2_KEEP: sample={} front={} back={} len={}",
                                task.sample_name, segment.front_kmer, segment.back_kmer, segment.data.len()
                            );
                        }
                        (segment.front_kmer, segment.back_kmer, false)
                    } else {
                        // Swap k-mers and reverse complement data
                        if config.verbosity > 2 {
                            #[cfg(feature = "verbose_debug")]
                            eprintln!(
                                "RAGC_CASE2_SWAP: sample={} front={} back={} -> key=({},{}) len={}",
                                task.sample_name, segment.front_kmer, segment.back_kmer,
                                segment.back_kmer, segment.front_kmer, segment.data.len()
                            );
                        }
                        (segment.back_kmer, segment.front_kmer, true)
                    }
                } else if segment.front_kmer != MISSING_KMER {
                    // Case 3a: Only front k-mer present, back is MISSING (terminator)
                    // Match C++ AGC lines 1315-1336: reverse complement and find candidate with one splitter
                    // Use the actual is_dir_oriented value from segment detection
                    #[cfg(feature = "verbose_debug")]
                    eprintln!("RAGC_CASE3A_TERMINATOR: sample={} front={} front_is_dir={} back=MISSING -> finding best group",
                        task.sample_name, segment.front_kmer, segment.front_kmer_is_dir);
                    // Debug: trace is_dir value before find_group call
                    if std::env::var("RAGC_DEBUG_IS_DIR").is_ok() {
                        eprintln!("RAGC_CASE3A_CALL: contig={} seg_part={} front_kmer={} front_kmer_is_dir={}",
                            task.contig_name, place, segment.front_kmer, segment.front_kmer_is_dir);
                    }
                    let (mut kf, mut kb, mut sr) = find_group_with_one_kmer(
                        segment.front_kmer,
                        segment.front_kmer_is_dir, // Use actual orientation from segment detection
                        &segment.data,
                        &segment_data_rc,
                        &map_segments_terminators,
                        &map_segments,
                        &segment_groups,
                        &reference_segments,
                        &config,
                    );

                    // Fallback: If Case 3a returned MISSING, try fallback minimizers (C++ AGC lines 1322-1334)
                    if (kf == MISSING_KMER || kb == MISSING_KMER) && fallback_filter.is_enabled() {
                        let (fb_kf, fb_kb, fb_sr) = find_cand_segment_using_fallback_minimizers(
                            &segment.data,
                            &segment_data_rc,
                            config.k,
                            5, // min_shared_kmers = 5 for Case 3 (matches C++ AGC)
                            &fallback_filter,
                            &map_fallback_minimizers,
                            &map_segments,
                            &segment_groups,
                            &reference_segments,
                            &config,
                        );
                        if fb_kf != MISSING_KMER && fb_kb != MISSING_KMER {
                            if config.verbosity > 1 {
                                #[cfg(feature = "verbose_debug")]
                                eprintln!("RAGC_CASE3A_FALLBACK: found ({},{}) rc={}", fb_kf, fb_kb, fb_sr);
                            }
                            kf = fb_kf;
                            kb = fb_kb;
                            sr = fb_sr;
                        }
                    }
                    (kf, kb, sr)
                } else if segment.back_kmer != MISSING_KMER {
                    // Case 3b: Only back k-mer present, front is MISSING (terminator)
                    // Match C++ AGC lines 1337-1360: swap_dir_rc() inverts is_dir_oriented()
                    //
                    // C++ AGC calls kmer.swap_dir_rc() which swaps kmer_dir and kmer_rc fields,
                    // effectively inverting is_dir_oriented() (which checks kmer_dir <= kmer_rc).
                    // So if back_kmer was originally dir-oriented, after swap it becomes NOT dir-oriented.
                    let kmer_is_dir_after_swap = !segment.back_kmer_is_dir;
                    #[cfg(feature = "verbose_debug")]
                    eprintln!("RAGC_CASE3B_TERMINATOR: sample={} front=MISSING back={} back_is_dir={} -> kmer_is_dir_after_swap={}",
                        task.sample_name, segment.back_kmer, segment.back_kmer_is_dir, kmer_is_dir_after_swap);

                    // C++ AGC line 1344 passes (segment_rc, segment) to find_cand_segment_with_one_splitter
                    // and then inverts the result: store_rc = !store_dir
                    // So we swap the segment parameters here AND invert sr below
                    let (mut kf, mut kb, mut sr) = find_group_with_one_kmer(
                        segment.back_kmer, // Use original k-mer value
                        kmer_is_dir_after_swap, // Inverted due to swap_dir_rc()
                        &segment_data_rc,  // SWAPPED: RC first (matches C++ AGC segment_rc param)
                        &segment.data,     // SWAPPED: Original second (matches C++ AGC segment param)
                        &map_segments_terminators,
                        &map_segments,
                        &segment_groups,
                        &reference_segments,
                        &config,
                    );
                    // Invert sr to match C++ AGC's store_rc = !store_dir
                    sr = !sr;

                    // Fallback: If Case 3b returned MISSING, try fallback minimizers (C++ AGC lines 1347-1359)
                    // Note: C++ AGC uses segment_rc for fallback in Case 3b
                    if (kf == MISSING_KMER || kb == MISSING_KMER) && fallback_filter.is_enabled() {
                        let (fb_kf, fb_kb, fb_sr) = find_cand_segment_using_fallback_minimizers(
                            &segment_data_rc, // Use RC for Case 3b (matches C++ AGC)
                            &segment.data,
                            config.k,
                            5, // min_shared_kmers = 5 for Case 3 (matches C++ AGC)
                            &fallback_filter,
                            &map_fallback_minimizers,
                            &map_segments,
                            &segment_groups,
                            &reference_segments,
                            &config,
                        );
                        if fb_kf != MISSING_KMER && fb_kb != MISSING_KMER {
                            if config.verbosity > 1 {
                                #[cfg(feature = "verbose_debug")]
                                eprintln!("RAGC_CASE3B_FALLBACK: found ({},{}) rc={}", fb_kf, fb_kb, !fb_sr);
                            }
                            kf = fb_kf;
                            kb = fb_kb;
                            sr = !fb_sr; // C++ AGC: store_rc = !store_dir_alt
                        }
                    }
                    (kf, kb, sr)
                } else {
                    // Case 1: Both MISSING - try fallback minimizers (C++ AGC lines 1286-1298)
                    let mut kf = MISSING_KMER;
                    let mut kb = MISSING_KMER;
                    let mut sr = false;

                    if fallback_filter.is_enabled() {
                        let (fb_kf, fb_kb, fb_sr) = find_cand_segment_using_fallback_minimizers(
                            &segment.data,
                            &segment_data_rc,
                            config.k,
                            1, // min_shared_kmers = 1 for Case 1 (matches C++ AGC line 1293)
                            &fallback_filter,
                            &map_fallback_minimizers,
                            &map_segments,
                            &segment_groups,
                            &reference_segments,
                            &config,
                        );
                        if fb_kf != MISSING_KMER && fb_kb != MISSING_KMER {
                            if config.verbosity > 1 {
                                #[cfg(feature = "verbose_debug")]
                                eprintln!("RAGC_CASE1_FALLBACK: sample={} found ({},{}) rc={} len={}",
                                    task.sample_name, fb_kf, fb_kb, fb_sr, segment.data.len());
                            }
                            kf = fb_kf;
                            kb = fb_kb;
                            sr = fb_sr;
                        }
                    }

                    (kf, kb, sr)
                };

            // Create grouping key from normalized k-mers
            // For raw segments (both k-mers MISSING), use round-robin distribution
            // across groups 0-15 to match C++ AGC behavior
            let key = if key_front == MISSING_KMER && key_back == MISSING_KMER {
                let raw_id = raw_group_counter.fetch_add(1, Ordering::SeqCst);
                let raw_group = raw_id % NO_RAW_GROUPS;
                SegmentGroupKey {
                    kmer_front: raw_group as u64,
                    kmer_back: MISSING_KMER,
                }
            } else {
                SegmentGroupKey {
                    kmer_front: key_front,
                    kmer_back: key_back,
                }
            };

            // Reverse complement data if needed (matching C++ AGC lines 1315-1316, 1320-1321)
            let segment_data = if should_reverse {
                segment.data.iter().rev().map(|&base| {
                    match base {
                        0 => 3, // A -> T
                        1 => 2, // C -> G
                        2 => 1, // G -> C
                        3 => 0, // T -> A
                        _ => base, // N or other non-ACGT
                    }
                }).collect()
            } else {
                segment.data.clone()
            };

            // Lock segment groups and add segment to buffer
            // NOTE: Split check must happen BEFORE creating BufferedSegment
            // to avoid moving segment_data prematurely
            {
                let mut groups = segment_groups.lock().unwrap();

                // Phase 1: Check if group already exists
                // (matches C++ AGC: seg_map_mtx.lock() then find at line 1020)
                let key_exists = {
                    let seg_map = map_segments.lock().unwrap();
                    seg_map.contains_key(&key)
                };

                // Phase 2: Try to split
                // C++ AGC only attempts splits when key doesn't exist (agc_compressor.cpp:1367)
                // This is the condition: p == map_segments.end() && both k-mers valid && both in terminators
                // Set RAGC_SPLIT_ALL=1 to try splitting even when key exists (experimental)
                // CRITICAL: C++ AGC lines 1374-1378 skip segment splitting when front == back!
                // When front == back, it just sets store_rc based on orientation, does NOT call
                // find_cand_segment_with_missing_middle_splitter. We must do the same.
                let split_allowed = if std::env::var("RAGC_SPLIT_ALL").as_deref() == Ok("1") { true } else { !key_exists };

                // Debug: trace split decision
                if std::env::var("RAGC_DEBUG_SPLIT").is_ok() && task.contig_name.contains("chrVII") && place >= 2 && place <= 5 {
                    eprintln!("RAGC_SPLIT_CHECK: contig={} seg={} key=({},{}) key_exists={} split_allowed={} front_missing={} back_missing={} front==back={}",
                        task.contig_name, place, key_front, key_back, key_exists, split_allowed,
                        key_front == MISSING_KMER, key_back == MISSING_KMER, key_front == key_back);
                }

                if split_allowed && key_front != MISSING_KMER && key_back != MISSING_KMER && key_front != key_back {
                    // CRITICAL: First attempt to find middle splitter
                    // Use ONLY global terminators (not batch-local) to match C++ AGC behavior
                    // C++ AGC only sees terminators from previous batches, not the current one
                    let middle_kmer_opt = {
                        let terminators = map_segments_terminators.lock().unwrap();
                        let result = find_middle_splitter(key_front, key_back, &terminators);
                        // Debug: trace middle splitter result
                        if std::env::var("RAGC_DEBUG_SPLIT").is_ok() && task.contig_name.contains("chrVII") && place >= 2 && place <= 5 {
                            let front_conn = terminators.get(&key_front).map(|v| v.len()).unwrap_or(0);
                            let back_conn = terminators.get(&key_back).map(|v| v.len()).unwrap_or(0);
                            eprintln!("RAGC_SPLIT_MIDDLE: contig={} seg={} key=({},{}) middle={:?} front_conn={} back_conn={}",
                                task.contig_name, place, key_front, key_back, result, front_conn, back_conn);
                        }
                        result
                    };

                    #[cfg(feature = "verbose_debug")]
                    if config.verbosity > 0 {
                        if middle_kmer_opt.is_some() {
                            eprintln!("DEBUG_SPLIT: Found middle k-mer for ({},{}) sample={}",
                                key_front, key_back, task.sample_name);
                        } else if config.verbosity > 1 {
                            eprintln!(
                                "SPLIT_NO_MIDDLE: ({},{}) sample={} place={} should_reverse={}",
                                key_front, key_back, task.sample_name, place, should_reverse
                            );
                        }
                    }

                    if let Some(middle_kmer) = middle_kmer_opt {
                        // Found potential middle k-mer
                        // Now check if BOTH split groups already exist in map_segments
                        // (This is the key difference from just checking terminators!)

                        // Debug: trace middle found
                        if std::env::var("RAGC_DEBUG_SPLIT").is_ok() {
                            eprintln!("RAGC_SPLIT_FOUND_MIDDLE: contig={} seg={} middle={}", task.contig_name, place, middle_kmer);
                        }

                        let left_key = if key_front <= middle_kmer {
                            SegmentGroupKey {
                                kmer_front: key_front,
                                kmer_back: middle_kmer,
                            }
                        } else {
                            SegmentGroupKey {
                                kmer_front: middle_kmer,
                                kmer_back: key_front,
                            }
                        };

                        let right_key = if middle_kmer <= key_back {
                            SegmentGroupKey {
                                kmer_front: middle_kmer,
                                kmer_back: key_back,
                            }
                        } else {
                            SegmentGroupKey {
                                kmer_front: key_back,
                                kmer_back: middle_kmer,
                            }
                        };

                        // CRITICAL: C++ AGC requires BOTH target groups to exist in map_segments
                        // at split decision time (agc_compressor.cpp lines 1472, 1486 use .at() which throws)
                        // If either group doesn't exist, C++ AGC aborts the split.
                        // We must check map_segments (global), not batch_local_groups, to match C++ behavior.
                        let (left_exists, right_exists) = {
                            let global_map = map_segments.lock().unwrap();
                            (global_map.contains_key(&left_key), global_map.contains_key(&right_key))
                        };

                        if !left_exists || !right_exists {
                            // Skip split - one or both target groups don't exist yet
                            // This matches C++ AGC behavior where .at() would throw
                            if config.verbosity > 1 {
                                eprintln!(
                                    "SPLIT_SKIP_NO_GROUP: left_key=({},{}) exists={} right_key=({},{}) exists={}",
                                    left_key.kmer_front, left_key.kmer_back, left_exists,
                                    right_key.kmer_front, right_key.kmer_back, right_exists
                                );
                            }
                            if std::env::var("RAGC_DEBUG_SPLIT").is_ok() {
                                eprintln!(
                                    "RAGC_SPLIT_SKIP_NO_GROUP: left=({},{}) exists={} right=({},{}) exists={}",
                                    left_key.kmer_front, left_key.kmer_back, left_exists,
                                    right_key.kmer_front, right_key.kmer_back, right_exists
                                );
                            }
                            // Don't attempt split - fall through to normal segment processing
                        } else {
                        // Both groups exist - proceed with split cost calculation
                        #[cfg(feature = "verbose_debug")]
                        if config.verbosity > 0 {
                            eprintln!("DEBUG_SPLIT: Attempting cost-based split for ({},{}) sample={}",
                                key_front, key_back, task.sample_name);
                        }

                        let split_result = try_split_segment_with_cost(
                            &segment_data,
                            key_front,
                            key_back,
                            middle_kmer,
                            &left_key,
                            &right_key,
                            &map_segments,
                            &map_segments_terminators,
                            &reference_segments,
                            &config,
                            should_reverse,
                        );

                        if let Some((left_data, right_data, _mid)) = split_result {
                            // Check if this is a degenerate split (one side empty)
                            let is_degenerate_left = left_data.is_empty();
                            let is_degenerate_right = right_data.is_empty();

                            if config.verbosity > 1 {
                                if is_degenerate_right {
                                    eprintln!(
                                        "SPLIT_DEGENERATE_RIGHT: ({},{}) -> left_only=({},{})",
                                        key_front, key_back, left_key.kmer_front, left_key.kmer_back
                                    );
                                } else if is_degenerate_left {
                                    eprintln!(
                                        "SPLIT_DEGENERATE_LEFT: ({},{}) -> right_only=({},{})",
                                        key_front, key_back, right_key.kmer_front, right_key.kmer_back
                                    );
                                } else {
                                    eprintln!(
                                        "SPLIT: original=({},{}) -> left=({},{}) right=({},{})",
                                        key_front, key_back, left_key.kmer_front, left_key.kmer_back,
                                        right_key.kmer_front, right_key.kmer_back
                                    );
                                }
                            }

                            // Determine emission order. By default match C++ logic:
                            // - Normal orientation (should_reverse=false): emit left then right
                            // - Reversed orientation (should_reverse=true): emit right then left
                            // Allow env override for diagnostics:
                            //   RAGC_EMIT_ORDER=left  -> force left-first
                            //   RAGC_EMIT_ORDER=right -> force right-first
                            //   RAGC_EMIT_ORDER=flip  -> invert default
                            //   RAGC_EMIT_ORDER=auto  -> default behavior (or if unset)
                            let emit_left_first = match std::env::var("RAGC_EMIT_ORDER") {
                                Ok(val) => match val.to_ascii_lowercase().as_str() {
                                    "left" | "left-first" => true,
                                    "right" | "right-first" => false,
                                    "flip" => should_reverse, // invert default (!should_reverse)
                                    _ => !should_reverse,      // auto/default
                                },
                                Err(_) => !should_reverse,
                            };
                            if config.verbosity > 1 {
                                eprintln!(
                                    "EMIT_ORDER: should_reverse={} -> emit_left_first={} (env RAGC_EMIT_ORDER)",
                                    should_reverse, emit_left_first
                                );
                            }

                            // Optional targeted split trace for a specific (sample, contig, index)
                            if let (Ok(ts), Ok(tc), Ok(ti)) = (
                                std::env::var("RAGC_TRACE_SAMPLE"),
                                std::env::var("RAGC_TRACE_CONTIG"),
                                std::env::var("RAGC_TRACE_INDEX").and_then(|s| s.parse::<usize>().map_err(|e| std::env::VarError::NotPresent)),
                            ) {
                                if ts == task.sample_name && tc == task.contig_name && ti == place {
                                    // Derive seg2_start from lengths (robust for both FFI and local mapping)
                                    let seg_len = segment_data.len();
                                    let right_len = right_data.len();
                                    let left_len = left_data.len();
                                    let seg2_start_derived = seg_len.saturating_sub(right_len);
                                    let left_end_derived = seg2_start_derived.saturating_add(config.k).min(seg_len);
                                    eprintln!(
                                        "TRACE_SPLIT: {}/{} idx={} rev={} emit_left_first={} degL={} degR={} seg2_start={} left_end={} left_len={} right_len={}",
                                        task.sample_name, task.contig_name, place, should_reverse, emit_left_first,
                                        is_degenerate_left, is_degenerate_right, seg2_start_derived, left_end_derived, left_len, right_len
                                    );
                                }
                            }

                            // Emit in correct contig order
                            if emit_left_first {
                                // left first
                                if !is_degenerate_left {
                                    let left_buffer = groups.entry(left_key.clone()).or_insert_with(|| {
                                        // BATCH-LOCAL: Check global first, then batch-local
                                        let group_id = {
                                            let global_map = map_segments.lock().unwrap();
                                            if let Some(&existing_id) = global_map.get(&left_key) {
                                                drop(global_map);
                                                existing_id
                                            } else {
                                                drop(global_map);
                                                let mut batch_map = batch_local_groups.lock().unwrap();
                                                *batch_map.entry(left_key.clone()).or_insert_with(|| {
                                                    group_counter.fetch_add(1, Ordering::SeqCst)
                                                })
                                            }
                                        };
                                        // Register with FFI engine
                                        #[cfg(feature = "cpp_agc")]
                                        if left_key.kmer_front != MISSING_KMER && left_key.kmer_back != MISSING_KMER {
                                            let mut eng = grouping_engine.lock().unwrap();
                                            eng.register_group(left_key.kmer_front, left_key.kmer_back, group_id);
                                        }
                                        // Update GLOBAL terminators map IMMEDIATELY (matches C++ AGC)
                                        if left_key.kmer_front != MISSING_KMER && left_key.kmer_back != MISSING_KMER {
                                            let mut term_map = map_segments_terminators.lock().unwrap();
                                            term_map.entry(left_key.kmer_front).or_insert_with(Vec::new).push(left_key.kmer_back);
                                            if left_key.kmer_front != left_key.kmer_back {
                                                term_map.entry(left_key.kmer_back).or_insert_with(Vec::new).push(left_key.kmer_front);
                                            }
                                            if let Some(front_vec) = term_map.get_mut(&left_key.kmer_front) { front_vec.sort_unstable(); front_vec.dedup(); }
                                            if left_key.kmer_front != left_key.kmer_back {
                                                if let Some(back_vec) = term_map.get_mut(&left_key.kmer_back) { back_vec.sort_unstable(); back_vec.dedup(); }
                                            }
                                        }
                                        // Register streams for this group
                                        let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
                                        let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
                                        let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);
                                        let mut arch = archive.lock().unwrap();
                                        let stream_id = arch.register_stream(&delta_stream_name);
                                        let ref_stream_id = arch.register_stream(&ref_stream_name);
                                        drop(arch);
                                        SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
                                    });
                                    let (fixed_left_data, fixed_left_rc) = fix_orientation_for_group(&left_data, should_reverse, &left_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let left_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: place, data: fixed_left_data, is_rev_comp: fixed_left_rc };
                                    left_buffer.segments.push(left_buffered);
                                    if left_buffer.segments.len() >= PACK_CARDINALITY + 1 { flush_pack(left_buffer, &collection, &archive, &config, &reference_segments).context("Failed to flush left pack")?; }
                                }
                                if !is_degenerate_right {
                                    let right_buffer = groups.entry(right_key.clone()).or_insert_with(|| {
                                        // BATCH-LOCAL: Check global first, then batch-local
                                        let group_id = {
                                            let global_map = map_segments.lock().unwrap();
                                            if let Some(&existing_id) = global_map.get(&right_key) {
                                                drop(global_map);
                                                existing_id
                                            } else {
                                                drop(global_map);
                                                let mut batch_map = batch_local_groups.lock().unwrap();
                                                *batch_map.entry(right_key.clone()).or_insert_with(|| {
                                                    group_counter.fetch_add(1, Ordering::SeqCst)
                                                })
                                            }
                                        };
                                        // Register with FFI engine
                                        #[cfg(feature = "cpp_agc")]
                                        if right_key.kmer_front != MISSING_KMER && right_key.kmer_back != MISSING_KMER {
                                            let mut eng = grouping_engine.lock().unwrap();
                                            eng.register_group(right_key.kmer_front, right_key.kmer_back, group_id);
                                        }
                                        // Update GLOBAL terminators map IMMEDIATELY (matches C++ AGC)
                                        if right_key.kmer_front != MISSING_KMER && right_key.kmer_back != MISSING_KMER {
                                            let mut term_map = map_segments_terminators.lock().unwrap();
                                            term_map.entry(right_key.kmer_front).or_insert_with(Vec::new).push(right_key.kmer_back);
                                            if right_key.kmer_front != right_key.kmer_back {
                                                term_map.entry(right_key.kmer_back).or_insert_with(Vec::new).push(right_key.kmer_front);
                                            }
                                            if let Some(front_vec) = term_map.get_mut(&right_key.kmer_front) { front_vec.sort_unstable(); front_vec.dedup(); }
                                            if right_key.kmer_front != right_key.kmer_back {
                                                if let Some(back_vec) = term_map.get_mut(&right_key.kmer_back) { back_vec.sort_unstable(); back_vec.dedup(); }
                                            }
                                        }
                                        // Register streams for this group
                                        let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
                                        let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
                                        let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);
                                        let mut arch = archive.lock().unwrap();
                                        let stream_id = arch.register_stream(&delta_stream_name);
                                        let ref_stream_id = arch.register_stream(&ref_stream_name);
                                        drop(arch);
                                        SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
                                    });
                                    let seg_part = if is_degenerate_left { place } else { place + 1 };
                                    let (fixed_right_data, fixed_right_rc) = fix_orientation_for_group(&right_data, should_reverse, &right_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let right_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: seg_part, data: fixed_right_data, is_rev_comp: fixed_right_rc };
                                    right_buffer.segments.push(right_buffered);
                                    if right_buffer.segments.len() >= PACK_CARDINALITY + 1 { flush_pack(right_buffer, &collection, &archive, &config, &reference_segments).context("Failed to flush right pack")?; }
                                }
                            } else {
                                // reversed: right first
                                if !is_degenerate_right {
                                    let right_buffer = groups.entry(right_key.clone()).or_insert_with(|| {
                                        // BATCH-LOCAL: Check global first, then batch-local (group must exist from earlier)
                                        let group_id = {
                                            let global_map = map_segments.lock().unwrap();
                                            if let Some(&id) = global_map.get(&right_key) {
                                                id
                                            } else {
                                                drop(global_map);
                                                let batch_map = batch_local_groups.lock().unwrap();
                                                *batch_map.get(&right_key).expect("Split right group must exist in batch_local_groups or map_segments")
                                            }
                                        };
                                        let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
                                        let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
                                        let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);
                                        let mut arch = archive.lock().unwrap();
                                        let stream_id = arch.register_stream(&delta_stream_name);
                                        let ref_stream_id = arch.register_stream(&ref_stream_name);
                                        drop(arch);
                                        SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
                                    });
                                    let (fixed_right_data, fixed_right_rc) = fix_orientation_for_group(&right_data, should_reverse, &right_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let right_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: place, data: fixed_right_data, is_rev_comp: fixed_right_rc };
                                    right_buffer.segments.push(right_buffered);
                                    if right_buffer.segments.len() >= PACK_CARDINALITY + 1 { flush_pack(right_buffer, &collection, &archive, &config, &reference_segments).context("Failed to flush right pack")?; }
                                }
                                if !is_degenerate_left {
                                    let left_buffer = groups.entry(left_key.clone()).or_insert_with(|| {
                                        // BATCH-LOCAL: Check global first, then batch-local (group must exist from earlier)
                                        let group_id = {
                                            let global_map = map_segments.lock().unwrap();
                                            if let Some(&id) = global_map.get(&left_key) {
                                                id
                                            } else {
                                                drop(global_map);
                                                let batch_map = batch_local_groups.lock().unwrap();
                                                *batch_map.get(&left_key).expect("Split left group must exist in batch_local_groups or map_segments")
                                            }
                                        };
                                        let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
                                        let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
                                        let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);
                                        let mut arch = archive.lock().unwrap();
                                        let stream_id = arch.register_stream(&delta_stream_name);
                                        let ref_stream_id = arch.register_stream(&ref_stream_name);
                                        drop(arch);
                                        SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
                                    });
                                    let seg_part = if is_degenerate_right { place } else { place + 1 };
                                    let (fixed_left_data, fixed_left_rc) = fix_orientation_for_group(&left_data, should_reverse, &left_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let left_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: seg_part, data: fixed_left_data, is_rev_comp: fixed_left_rc };
                                    left_buffer.segments.push(left_buffered);
                                    if left_buffer.segments.len() >= PACK_CARDINALITY + 1 { flush_pack(left_buffer, &collection, &archive, &config, &reference_segments).context("Failed to flush left pack")?; }
                                }
                            }

                            // Optional: assert lengths vs C++ archive if provided
                            if let Ok(assert_path) = std::env::var("RAGC_ASSERT_CPP_ARCHIVE") {
                                use crate::{Decompressor, DecompressorConfig};
                                let mut dec = match Decompressor::open(&assert_path, DecompressorConfig{ verbosity: 0 }) {
                                    Ok(d) => d, Err(_) => {
                                        if config.verbosity > 1 { eprintln!("ASSERT_SKIP: cannot open {}", assert_path); }
                                        return Ok(());
                                    }
                                };
                                if let Ok(all) = dec.get_all_segments() {
                                    if let Some((_, _, segs)) = all.into_iter().find(|(s,c,_)| *s == task.sample_name && *c == task.contig_name) {
                                        // Compute our emitted lens and expected lens at indices
                                        let mut checks: Vec<(usize, usize)> = Vec::new();
                                        if emit_left_first {
                                            if !is_degenerate_left { checks.push((place, left_data.len())); }
                                            if !is_degenerate_right { checks.push((if is_degenerate_left { place } else { place + 1 }, right_data.len())); }
                                        } else {
                                            if !is_degenerate_right { checks.push((place, right_data.len())); }
                                            if !is_degenerate_left { checks.push((if is_degenerate_right { place } else { place + 1 }, left_data.len())); }
                                        }

                                        // Derive segmentation geometry for detailed diagnostics
                                        let seg_len = segment_data.len();
                                        let right_len = right_data.len();
                                        let left_len = left_data.len();
                                        let seg2_start_derived = seg_len.saturating_sub(right_len);
                                        let left_end_derived = seg2_start_derived.saturating_add(config.k).min(seg_len);
                                        let emit_idx_left = if emit_left_first { place } else { if is_degenerate_right { place } else { place + 1 } };
                                        let emit_idx_right = if emit_left_first { if is_degenerate_left { place } else { place + 1 } } else { place };

                                        for (idx, got) in checks {
                                            if idx < segs.len() {
                                                let exp = segs[idx].raw_length as usize;
                                                if exp != got {
                                                    eprintln!("ASSERT_LEN_MISMATCH: {}/{} idx={} got={} exp={} keys L=({:#x},{:#x}) R=({:#x},{:#x})",
                                                        task.sample_name, task.contig_name, idx, got, exp,
                                                        left_key.kmer_front, left_key.kmer_back,
                                                        right_key.kmer_front, right_key.kmer_back);
                                                    // Extended context (guarded by env to limit noise)
                                                    if std::env::var("RAGC_ASSERT_VERBOSE").is_ok() {
                                                        eprintln!("  CONTEXT: place={} orig_place={} emit_left_first={} should_reverse={}",
                                                            place, original_place, emit_left_first, should_reverse);
                                                        eprintln!("  GEOM: seg_len={} left_len={} right_len={} seg2_start={} left_end={}",
                                                            seg_len, left_len, right_len, seg2_start_derived, left_end_derived);
                                                        eprintln!("  EMIT_IDX: left_at={} right_at={}", emit_idx_left, emit_idx_right);
                                                    }
                                                }
                                            } else {
                                                eprintln!("ASSERT_IDX_OOB: {}/{} idx={} (segs={})", task.sample_name, task.contig_name, idx, segs.len());
                                            }
                                        }
                                    }
                                }
                            }

                            // Record this split so subsequent segments from this contig get shifted
                            // (matches C++ AGC lines 2033-2036: ++seg_part_no twice when split)
                            // For degenerate splits, only increment once (no actual split)
                            if !is_degenerate_left && !is_degenerate_right {
                                let mut offsets = split_offsets.lock().unwrap();
                                offsets.insert((task.sample_name.clone(), task.contig_name.clone(), original_place), 1);
                            }

                            // Skip adding original segment - we've added the split/reclassified segment
                            continue;
                        }
                        // If split_result was None, fall through to normal path
                        } // end of else { both groups exist }
                    }
                }

                // Phase 2.5: Secondary fallback attempt (C++ AGC lines 1477-1494)
                // If the group doesn't exist yet, try fallback minimizers one more time with min_shared=2
                // This helps segments find existing groups that share internal k-mers
                let (key, key_front, key_back, should_reverse) = {
                    // Re-check if key exists (may have changed since split logic ran)
                    let key_exists_now = {
                        let seg_map = map_segments.lock().unwrap();
                        seg_map.contains_key(&key)
                    };

                    // Debug: count how many segments could be eligible for secondary fallback
                    if std::env::var("RAGC_DEBUG_FALLBACK2").as_deref() == Ok("1") {
                        if !key_exists_now && key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                            eprintln!("SECONDARY_FB_CANDIDATE: sample={} contig={} place={} key=({},{})",
                                task.sample_name, task.contig_name, place, key.kmer_front, key.kmer_back);
                        }
                    }

                    if !key_exists_now
                        && key.kmer_front != MISSING_KMER
                        && key.kmer_back != MISSING_KMER
                        && fallback_filter.is_enabled()
                    {
                        // Generate reverse complement for fallback lookup
                        let segment_data_rc_fb: Vec<u8> = segment_data.iter().rev().map(|&b| {
                            if b > 3 { b } else { 3 - b }
                        }).collect();

                        let (fb_kf, fb_kb, fb_sr) = find_cand_segment_using_fallback_minimizers(
                            &segment_data,
                            &segment_data_rc_fb,
                            config.k,
                            2, // min_shared_kmers = 2 for secondary fallback (C++ AGC line 1482)
                            &fallback_filter,
                            &map_fallback_minimizers,
                            &map_segments,
                            &segment_groups,
                            &reference_segments,
                            &config,
                        );

                        if std::env::var("RAGC_DEBUG_FALLBACK2").as_deref() == Ok("1") {
                            if fb_kf == MISSING_KMER || fb_kb == MISSING_KMER {
                                eprintln!("SECONDARY_FB_NO_MATCH: orig_key=({},{})", key.kmer_front, key.kmer_back);
                            } else {
                                eprintln!("SECONDARY_FB_FOUND: orig=({},{}) found=({},{}) rc={}",
                                    key.kmer_front, key.kmer_back, fb_kf, fb_kb, fb_sr);
                            }
                        }

                        if fb_kf != MISSING_KMER && fb_kb != MISSING_KMER {
                            // Verify the found group actually exists
                            let found_key = SegmentGroupKey {
                                kmer_front: fb_kf,
                                kmer_back: fb_kb,
                            };
                            let found_exists = {
                                let seg_map = map_segments.lock().unwrap();
                                seg_map.contains_key(&found_key)
                            };

                            if found_exists {
                                if config.verbosity > 1 {
                                    eprintln!("SECONDARY_FALLBACK_SUCCESS: ({},{}) -> ({},{}) sr={}->{}",
                                        key_front, key_back, fb_kf, fb_kb, should_reverse, fb_sr);
                                }
                                (found_key, fb_kf, fb_kb, fb_sr)
                            } else {
                                // Fallback found k-mers but group doesn't exist - keep original
                                (key, key_front, key_back, should_reverse)
                            }
                        } else {
                            // Fallback didn't find anything - keep original
                            (key, key_front, key_back, should_reverse)
                        }
                    } else {
                        // Group exists or not eligible for fallback - keep original
                        (key, key_front, key_back, should_reverse)
                    }
                };

                // Phase 3: Normal path - add segment to group as-is (group exists, or split failed/impossible)

                // FIX: When joining an existing group, use the reference's orientation
                // to ensure LZ encoding works correctly (fixes ZERO_MATCH bug in Case 3 terminators)
                let (final_should_reverse, final_segment_data) = {
                    // Check if this group already has a reference with stored orientation
                    // Check BOTH map_segments (previous batches) AND batch_local_groups (current batch)
                    let group_id_opt = {
                        let seg_map = map_segments.lock().unwrap();
                        if let Some(&gid) = seg_map.get(&key) {
                            Some(gid)
                        } else {
                            drop(seg_map);
                            let batch_map = batch_local_groups.lock().unwrap();
                            batch_map.get(&key).copied()
                        }
                    };

                    if let Some(group_id) = group_id_opt {
                        let ref_orients = reference_orientations.lock().unwrap();
                        if let Some(&ref_rc) = ref_orients.get(&group_id) {
                            // Group has a stored reference orientation
                            if ref_rc != should_reverse {
                                // Orientation mismatch - recompute segment_data with reference orientation
                                if config.verbosity > 1 {
                                    eprintln!("ORIENTATION_FIX: group={} sample={} contig={}:{} computed_rc={} ref_rc={}",
                                        group_id, task.sample_name, task.contig_name, place, should_reverse, ref_rc);
                                }
                                // Recompute segment_data with the reference's orientation
                                // Note: We need to reverse the ORIGINAL data, not segment_data which may already be reversed
                                let fixed_data = if ref_rc {
                                    // Need RC: reverse complement the original segment data
                                    segment.data.iter().rev().map(|&base| {
                                        match base {
                                            0 => 3, // A -> T
                                            1 => 2, // C -> G
                                            2 => 1, // G -> C
                                            3 => 0, // T -> A
                                            _ => base, // N or other non-ACGT
                                        }
                                    }).collect()
                                } else {
                                    // Need forward: use original data as-is
                                    segment.data.clone()
                                };
                                (ref_rc, fixed_data)
                            } else {
                                // Orientations match - use computed values
                                (should_reverse, segment_data)
                            }
                        } else {
                            // No stored orientation yet (this will be the reference)
                            (should_reverse, segment_data)
                        }
                    } else {
                        // New group - use computed orientation
                        (should_reverse, segment_data)
                    }
                };

                // Create buffered segment with corrected orientation
                let buffered = BufferedSegment {
                    sample_name: task.sample_name.clone(),
                    contig_name: task.contig_name.clone(),
                    seg_part_no: place,
                    data: final_segment_data,
                    is_rev_comp: final_should_reverse,
                };

                // DEBUG: Print data for degenerate k-mer reference segments
                if key_front == 1244212049458757632 && key_back == 1244212049458757632 && config.verbosity > 1 {
                    let data_first5: Vec<u8> = buffered.data.iter().take(5).cloned().collect();
                    eprintln!("RAGC_DEGENERATE_BUFFERED: sample={} contig={} part={} key=({},{}) should_reverse={} len={} data[0..5]={:?}",
                        task.sample_name, task.contig_name, place, key_front, key_back, should_reverse, buffered.data.len(), data_first5);
                }

                // Get or create buffer for this group
                // If group doesn't exist yet, allocate group_id using BATCH-LOCAL logic
                let buffer = groups.entry(key.clone()).or_insert_with(|| {
                    // BATCH-LOCAL GROUP ASSIGNMENT (matches C++ AGC m_kmers per-batch behavior)
                    // 1. Check GLOBAL registry first (for groups from previous batches)
                    // 2. If not found globally, check/create in BATCH-LOCAL map
                    let group_id = {
                        let global_map = map_segments.lock().unwrap();
                        if let Some(&existing_id) = global_map.get(&key) {
                            // Group exists from previous batch - reuse it
                            drop(global_map);
                            existing_id
                        } else {
                            // Not in global registry - check batch-local map
                            drop(global_map);
                            let mut batch_map = batch_local_groups.lock().unwrap();
                            *batch_map.entry(key.clone()).or_insert_with(|| {
                                // New group for this batch - allocate ID
                                if key.kmer_back == MISSING_KMER && key.kmer_front < NO_RAW_GROUPS as u64 {
                                    key.kmer_front as u32  // Raw groups use ID from key
                                } else {
                                    group_counter.fetch_add(1, Ordering::SeqCst)  // LZ groups use counter
                                }
                            })
                        }
                    };

                    // Register with FFI engine
                    #[cfg(feature = "cpp_agc")]
                    if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                        let mut eng = grouping_engine.lock().unwrap();
                        eng.register_group(key.kmer_front, key.kmer_back, group_id);
                    }

                    // Update BATCH-LOCAL terminators (NOT global) - C++ AGC agc_compressor.cpp:1015-1024
                    //
                    // CRITICAL: C++ AGC updates terminators in store_segments (called at batch END),
                    // NOT during segment processing. This means segments in the SAME batch cannot
                    // see each other through terminators - they only see groups from PREVIOUS batches.
                    //
                    // Previous code incorrectly updated GLOBAL terminators immediately, causing
                    // later segments in a batch to find and use groups created earlier in the same
                    // batch. This produced different grouping than C++ AGC.
                    //
                    // Now we update batch_local_terminators, which is merged to global at batch end
                    // via flush_batch() - matching C++ AGC's timing exactly.
                    //
                    // C++ AGC also requires BOTH k-mers to be non-MISSING (line 1015):
                    //   if (kmer1 != ~0ull && kmer2 != ~0ull)
                    if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                        let mut term_map = batch_local_terminators.lock().unwrap();

                        // Add bidirectional edge: front -> back
                        term_map.entry(key.kmer_front)
                            .or_insert_with(Vec::new)
                            .push(key.kmer_back);

                        // Add bidirectional edge: back -> front (if different from front)
                        if key.kmer_front != key.kmer_back {
                            term_map.entry(key.kmer_back)
                                .or_insert_with(Vec::new)
                                .push(key.kmer_front);
                        }

                        // Sort to maintain sorted order for set_intersection
                        if let Some(front_vec) = term_map.get_mut(&key.kmer_front) {
                            front_vec.sort_unstable();
                            front_vec.dedup();
                        }
                        if key.kmer_front != key.kmer_back {
                            if let Some(back_vec) = term_map.get_mut(&key.kmer_back) {
                                back_vec.sort_unstable();
                                back_vec.dedup();
                            }
                        }
                    }

                    if config.verbosity > 1 {
                        eprintln!("NEW_GROUP: group_id={} front={} back={} sample={}",
                            group_id, key_front, key_back, task.sample_name);
                    }

                    // Register streams for this group
                    let archive_version =
                        ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
                    let delta_stream_name =
                        ragc_common::stream_delta_name(archive_version, group_id);
                    let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);

                    let mut arch = archive.lock().unwrap();
                    let stream_id = arch.register_stream(&delta_stream_name);
                    let ref_stream_id = arch.register_stream(&ref_stream_name);
                    drop(arch);

                    // NOTE: Terminator updates already done in batch_local_terminators above
                    // Removed redundant duplicate block here

                    SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
                });

                // CSV logging for ALL segment grouping decisions (enabled via RAGC_GROUP_LOG=1)
                // This logs every segment, including those joining existing batch-local groups
                if std::env::var("RAGC_GROUP_LOG").as_deref() == Ok("1") {
                    let source = if key_exists { "global" } else { "batch_or_new" };
                    eprintln!("SEGMENT_GROUP,{},{},{},{},{},{},{},{}",
                        task.sample_name, task.contig_name, place,
                        key.kmer_front, key.kmer_back, final_should_reverse, source, buffer.group_id);
                }

                // Add fallback mapping for this segment (matches C++ AGC add_fallback_mapping)
                // This populates the fallback minimizers map for use by find_cand_segment_using_fallback_minimizers
                add_fallback_mapping(
                    &buffered.data,
                    config.k,
                    key.kmer_front,
                    key.kmer_back,
                    &fallback_filter,
                    &map_fallback_minimizers,
                );

                if key_exists && config.verbosity > 1 {
                    eprintln!("REUSE_GROUP: group_id={} front={} back={} sample={}",
                        buffer.group_id, key_front, key_back, task.sample_name);
                }

                // CRITICAL: Match C++ AGC behavior - write reference IMMEDIATELY for LZ groups
                // (C++ AGC segment.cpp lines 41-48: if (no_seqs == 0) writes reference right away)
                // BUT: Raw groups (0-15) don't have references - they skip position 0
                let is_raw_group = buffer.group_id < NO_RAW_GROUPS;
                if !is_raw_group && buffer.reference_segment.is_none() && buffer.segments.is_empty() {
                    // This is the FIRST segment in an LZ group - make it the reference NOW
                    // (matches C++ AGC: lz_diff->Prepare(s); store_in_archive(s, zstd_cctx);)
                    if let Err(e) = write_reference_immediately(&buffered, buffer, &collection, &archive, &reference_segments, &reference_orientations, &config) {
                        eprintln!("ERROR: Failed to write immediate reference: {}", e);
                        // Fall back to buffering
                        buffer.segments.push(buffered);
                    }
                } else {
                    // Buffer segment: either raw group, or subsequent LZ segment
                    buffer.segments.push(buffered);
                }

                // Flush pack if buffer is full
                if buffer.segments.len() >= PACK_CARDINALITY + 1 {
                    flush_pack(buffer, &collection, &archive, &config, &reference_segments)
                        .context("Failed to flush pack")?;
                }
            }
        }

        processed_count += 1;
    }

    Ok(())
}

// ========== SEGMENT SPLITTING HELPER FUNCTIONS ==========
// (Phase 3-6 implementation)

/// Phase 3: Find a k-mer that connects both front and back
/// Returns the first k-mer that appears in the terminator lists of BOTH front and back
/// (matches C++ AGC find_cand_segment_with_missing_middle_splitter lines 1531-1554)
fn find_middle_splitter(
    front_kmer: u64,
    back_kmer: u64,
    terminators: &HashMap<u64, Vec<u64>>,
) -> Option<u64> {
    let front_connections = terminators.get(&front_kmer)?;
    let back_connections = terminators.get(&back_kmer)?;

    #[cfg(feature = "cpp_agc")]
    {
        if let Some(m) = crate::ragc_ffi::find_middle(front_connections, back_connections) {
            return Some(m);
        }
        if std::env::var("RAGC_DEBUG_SPLIT_FIND").is_ok() {
            eprintln!(
                "DEBUG_FIND_MIDDLE_MISS: front={} back={} front_conn={} back_conn={} shared=0",
                front_kmer, back_kmer, front_connections.len(), back_connections.len()
            );
        }
        None
    }

    #[cfg(not(feature = "cpp_agc"))]
    {
        // Fallback: local set_intersection
        let mut i = 0;
        let mut j = 0;
        while i < front_connections.len() && j < back_connections.len() {
            let a = front_connections[i];
            let b = back_connections[j];
            if a == b {
                if a != MISSING_KMER { return Some(a); }
                i += 1; j += 1;
            } else if a < b {
                i += 1;
            } else {
                j += 1;
            }
        }
        if std::env::var("RAGC_DEBUG_SPLIT_FIND").is_ok() {
            eprintln!(
                "DEBUG_FIND_MIDDLE_MISS: front={} back={} front_conn={} back_conn={} shared=0",
                front_kmer, back_kmer, front_connections.len(), back_connections.len()
            );
            eprintln!("  front_connections: {:?}", &front_connections[..front_connections.len().min(5)]);
            eprintln!("  back_connections: {:?}", &back_connections[..back_connections.len().min(5)]);
        }
        None
    }
}

/// Phase 4: Find split position by scanning for middle k-mer
/// Scans the segment to find where the middle k-mer actually occurs
/// Returns the split position (in bytes) at the END of the middle k-mer
fn find_split_position(segment_data: &[u8], middle_kmer: u64, segment_len: usize, k: usize) -> Option<usize> {
    use crate::kmer::{Kmer, KmerMode};

    // Ensure we don't split too close to the ends
    // Need at least k+1 bytes on each side for valid segments
    if segment_len < 2 * (k + 1) {
        return None;
    }

    // Scan segment to find where middle_kmer occurs
    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);

    for (pos, &base) in segment_data.iter().enumerate() {
        kmer.insert(base as u64);

        if kmer.is_full() {
            let current_kmer = kmer.data_canonical();
            if current_kmer == middle_kmer {
                // Found the middle k-mer! Position is at the end of the k-mer
                let split_pos = pos + 1;

                // Validate: ensure we have enough space on both sides
                let left_size = split_pos;
                let right_size = segment_len - split_pos + k;

                if left_size >= k + 1 && right_size >= k + 1 {
                    return Some(split_pos);
                }
            }
        }
    }

    // Middle k-mer not found in segment - shouldn't happen but handle gracefully
    None
}

/// Phase 5: Split segment into two overlapping segments
/// Returns (left_segment, right_segment) with k-mer overlap
/// (matches C++ AGC lines 1461-1464)
fn split_segment_at_position(
    segment_data: &[u8],
    split_pos: usize,
    k: usize,
) -> (Vec<u8>, Vec<u8>) {
    // C++ AGC creates overlap of k bytes (not k/2!):
    //   seg2_start_pos = left_size - ceil(kmer_length / 2)
    //   segment2 starts at seg2_start_pos
    //   segment ends at seg2_start_pos + kmer_length
    // This creates k bytes of overlap: [split_pos - k/2 .. split_pos + k/2]
    let half_ceil = (k + 1) / 2;
    let seg2_start_pos = split_pos.saturating_sub(half_ceil);

    // Right segment: [seg2_start_pos .. end]
    let right = segment_data[seg2_start_pos..].to_vec();

    // Left segment: [0 .. seg2_start_pos + k]
    let left_end = seg2_start_pos + k;
    let left = segment_data[..left_end].to_vec();

    (left, right)
}

/// Split using seg2_start byte index (start of right segment) matching C++ layout
fn split_segment_from_start(segment_data: &[u8], seg2_start: usize, k: usize) -> (Vec<u8>, Vec<u8>) {
    let seg2_start_pos = seg2_start.min(segment_data.len());
    let right = segment_data[seg2_start_pos..].to_vec();
    let left_end = seg2_start_pos.saturating_add(k).min(segment_data.len());
    let left = segment_data[..left_end].to_vec();
    (left, right)
}

/// Phase 6: Attempt to split using compression cost heuristic (EXACT C++ AGC algorithm)
/// Matches agc_compressor.cpp lines 1387-1503 and 1531-1663
/// Returns Some((left_data, right_data, middle_kmer)) if split is beneficial
/// Returns None if split would be degenerate (creates segments too small)
fn try_split_segment_with_cost(
    segment_data: &Contig,
    front_kmer: u64,
    back_kmer: u64,
    middle_kmer: u64,
    left_key: &SegmentGroupKey,
    right_key: &SegmentGroupKey,
    map_segments: &Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: &Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    reference_segments: &Arc<Mutex<HashMap<u32, Vec<u8>>>>,
    config: &StreamingQueueConfig,
    should_reverse: bool,
) -> Option<(Vec<u8>, Vec<u8>, u64)> {
    if config.verbosity > 1 {
        eprintln!(
            "SPLIT_ATTEMPT: front={} back={} middle={}",
            front_kmer, back_kmer, middle_kmer
        );
    }

    // Debug: trace split attempt
    if std::env::var("RAGC_DEBUG_SPLIT").is_ok() {
        eprintln!("RAGC_SPLIT_TRY: front={} back={} middle={} left_key=({},{}) right_key=({},{})",
            front_kmer, back_kmer, middle_kmer,
            left_key.kmer_front, left_key.kmer_back,
            right_key.kmer_front, right_key.kmer_back);
    }

    // Prepare LZDiff for both groups from persistent storage
    // C++ AGC uses global v_segments[segment_id] (agc_compressor.cpp:1535-1536)
    // RAGC uses reference_segments HashMap - ALWAYS prepare on-demand
    // Don't require groups to be in local buffer (other workers may have created them)

    // Helper to prepare LZDiff from global reference_segments
    let prepare_on_demand = |key: &SegmentGroupKey, label: &str| -> Option<LZDiff> {
        let map_segments_locked = map_segments.lock().unwrap();
        let ref_segments_locked = reference_segments.lock().unwrap();

        // C++ AGC uses map_segments[key] which returns 0 (default) if key doesn't exist.
        // v_segments[0] is a raw group initialized with empty_ctg = { 0x7f } (1 byte).
        // This gives maximum LZ cost (no compression) for non-existent groups.
        // To match C++ AGC behavior, use the actual reference if available,
        // otherwise use empty reference (gives max cost like C++ AGC's v_segments[0]).
        let segment_id = map_segments_locked.get(key).copied();

        if let Some(ref_data) = segment_id.and_then(|id| ref_segments_locked.get(&id)) {
            // Reference exists! Prepare LZDiff on-demand
            if std::env::var("RAGC_DEBUG_SPLIT_REF").is_ok() {
                eprintln!(
                    "RAGC_SPLIT_REF: {}_key=({},{}) segment_id={:?} ref_size={} (ACTUAL)",
                    label, key.kmer_front, key.kmer_back, segment_id, ref_data.len()
                );
            }

            let mut lz = LZDiff::new(config.min_match_len as u32);
            lz.prepare(ref_data);
            return Some(lz);
        } else {
            // No reference data available for this group
            // C++ AGC uses v_segments[0] which is initialized with empty_ctg = { 0x7f } (1 byte)
            // This gives maximum LZ cost (no compression matches possible)
            // Return LZDiff prepared with empty reference to match C++ AGC behavior
            if std::env::var("RAGC_DEBUG_SPLIT_REF").is_ok() {
                eprintln!(
                    "RAGC_SPLIT_REF: {}_key=({},{}) segment_id={:?} ref_size=1 (EMPTY FALLBACK)",
                    label, key.kmer_front, key.kmer_back, segment_id
                );
            }

            // Use 1-byte dummy reference like C++ AGC's empty_ctg = { 0x7f }
            let empty_ref: Vec<u8> = vec![0x7f];
            let mut lz = LZDiff::new(config.min_match_len as u32);
            lz.prepare(&empty_ref);
            return Some(lz);
        }
    };

    // Build segment in both orientations once
    let segment_dir = segment_data; // &Vec<u8>
    // Reverse-complement once
    let segment_rc_vec: Vec<u8> = reverse_complement_sequence(segment_data);

    // Calculate compression costs and best split position using C++ FFI if enabled
    // Falls back to Rust implementation otherwise
    let mut maybe_best: Option<(usize, usize)> = None; // (best_pos, seg2_start)
    #[cfg(feature = "cpp_agc")]
    {
        // Inspect availability of left/right references and log keys
        let (left_seg_id_opt, right_seg_id_opt) = {
            let map_segments_locked = map_segments.lock().unwrap();
            (map_segments_locked.get(left_key).copied(), map_segments_locked.get(right_key).copied())
        };
        let (left_have_ref, right_have_ref) = {
            let ref_segments_locked = reference_segments.lock().unwrap();
            (left_seg_id_opt.and_then(|id| ref_segments_locked.get(&id)).is_some(),
             right_seg_id_opt.and_then(|id| ref_segments_locked.get(&id)).is_some())
        };

        if config.verbosity > 1 {
            eprintln!(
                "SPLIT_KEYS: left=({:#x},{:#x}) right=({:#x},{:#x}) left_seg_id={:?} right_seg_id={:?} have_left_ref={} have_right_ref={}",
                left_key.kmer_front, middle_kmer, middle_kmer, right_key.kmer_back,
                left_seg_id_opt, right_seg_id_opt, left_have_ref, right_have_ref
            );
        }

        // Prepare neighbor lists for FFI decision
        let (front_neighbors, back_neighbors) = {
            let term_map = map_segments_terminators.lock().unwrap();
            (term_map.get(&front_kmer).cloned().unwrap_or_default(),
             term_map.get(&back_kmer).cloned().unwrap_or_default())
        };

        // Always attempt FFI decision; if refs are missing, C++ will decide no-split
        let (ref_left_opt, ref_right_opt) = {
            let ref_segments_locked = reference_segments.lock().unwrap();
            let l = left_seg_id_opt.and_then(|id| ref_segments_locked.get(&id).cloned());
            let r = right_seg_id_opt.and_then(|id| ref_segments_locked.get(&id).cloned());
            (l, r)
        };
        let empty: Vec<u8> = Vec::new();
        let ref_left = ref_left_opt.as_ref().unwrap_or(&empty);
        let ref_right = ref_right_opt.as_ref().unwrap_or(&empty);

        if let Some((has_mid, mid, bp, s2, should)) = crate::ragc_ffi::decide_split(
            &front_neighbors, &back_neighbors,
            ref_left, ref_right,
            segment_dir,
            front_kmer, back_kmer,
            config.min_match_len as u32,
            config.k as u32,
            should_reverse,
        ) {
            if config.verbosity > 1 {
                eprintln!("FFI_DECIDE: has_middle={} middle={:#x} best_pos={} seg2_start={} should_split={} refs L={} R={}", has_mid, mid, bp, s2, should, ref_left.len(), ref_right.len());
            }
            if !has_mid || !should { return None; }
            maybe_best = Some((bp, s2));
        } else if config.verbosity > 1 {
            eprintln!("FFI_DECIDE: unavailable (decide_split returned None)");
        }
    }

    // If FFI provided best position, use it; otherwise compute costs in Rust
    let mut v_costs1 = if maybe_best.is_none() {
        if let Some(mut lz_left) = prepare_on_demand(left_key, "left") {
        #[cfg(feature = "cpp_agc")]
        {
            // Unused path when FFI returns best split; kept for completeness
            let ref_left = {
                let map_segments_locked = map_segments.lock().unwrap();
                let ref_segments_locked = reference_segments.lock().unwrap();
                let seg_id = map_segments_locked.get(left_key).copied().unwrap_or(0);
                ref_segments_locked.get(&seg_id).cloned()
            };
            if let Some(ref_data) = ref_left {
                if front_kmer < middle_kmer {
                    crate::ragc_ffi::cost_vector(true, &ref_data, segment_dir, config.min_match_len as u32)
                } else {
                    let mut v = crate::ragc_ffi::cost_vector(false, &ref_data, &segment_rc_vec, config.min_match_len as u32);
                    v.reverse(); v
                }
            } else {
                if config.verbosity > 1 { eprintln!("SPLIT_SKIP: left group has no reference yet"); }
                return None;
            }
        }
        #[cfg(not(feature = "cpp_agc"))]
        {
            if front_kmer < middle_kmer {
                lz_left.get_coding_cost_vector(segment_dir, true)
            } else {
                let mut v = lz_left.get_coding_cost_vector(&segment_rc_vec, false);
                v.reverse(); v
            }
        }
        } else {
        if config.verbosity > 1 { eprintln!("SPLIT_SKIP: left group has no reference yet"); }
        if std::env::var("RAGC_DEBUG_SPLIT").is_ok() {
            eprintln!("RAGC_SPLIT_SKIP_LEFT: left_key=({},{}) has no reference", left_key.kmer_front, left_key.kmer_back);
        }
        return None;
        }
    } else { Vec::new() };

    // Cumulative sum forward for v_costs1
    let mut sum = 0u32;
    for cost in v_costs1.iter_mut() {
        sum = sum.saturating_add(*cost);
        *cost = sum;
    }

    let mut v_costs2 = if maybe_best.is_none() {
        if let Some(mut lz_right) = prepare_on_demand(right_key, "right") {
        #[cfg(feature = "cpp_agc")]
        {
            let ref_right = {
                let map_segments_locked = map_segments.lock().unwrap();
                let ref_segments_locked = reference_segments.lock().unwrap();
                let seg_id = map_segments_locked.get(right_key).copied().unwrap_or(0);
                ref_segments_locked.get(&seg_id).cloned()
            };
            if let Some(ref_data) = ref_right {
                let mut v = if middle_kmer < back_kmer {
                    // Suffix placement, cumulative sum right-to-left
                    crate::ragc_ffi::cost_vector(false, &ref_data, segment_dir, config.min_match_len as u32)
                } else {
                    // RC + prefix placement; cumulative sum left-to-right then reverse
                    crate::ragc_ffi::cost_vector(true, &ref_data, &segment_rc_vec, config.min_match_len as u32)
                };
                if middle_kmer < back_kmer {
                    // Reverse cumulative sum
                    let mut acc = 0u32;
                    for cost in v.iter_mut().rev() {
                        acc = acc.saturating_add(*cost);
                        *cost = acc;
                    }
                    v
                } else {
                    // Forward cumulative then reverse
                    let mut acc = 0u32;
                    for cost in v.iter_mut() {
                        acc = acc.saturating_add(*cost);
                        *cost = acc;
                    }
                    v.reverse();
                    v
                }
            } else {
                if config.verbosity > 1 { eprintln!("SPLIT_SKIP: right group has no reference yet"); }
                return None;
            }
        }
        #[cfg(not(feature = "cpp_agc"))]
        {
            if middle_kmer < back_kmer {
                let mut v = lz_right.get_coding_cost_vector(segment_dir, false);
                let mut acc = 0u32; for cost in v.iter_mut().rev() { acc = acc.saturating_add(*cost); *cost = acc; } v
            } else {
                let mut v = lz_right.get_coding_cost_vector(&segment_rc_vec, true);
                let mut acc = 0u32; for cost in v.iter_mut() { acc = acc.saturating_add(*cost); *cost = acc; } v.reverse(); v
            }
        }
        } else {
        if config.verbosity > 1 { eprintln!("SPLIT_SKIP: right group has no reference yet"); }
        return None;
        }
    } else { Vec::new() };

    if maybe_best.is_none() && (v_costs1.is_empty() || v_costs2.is_empty()) {
        if config.verbosity > 1 {
            eprintln!("SPLIT_SKIP: cost vectors empty");
        }
        return None;
    }

    if maybe_best.is_none() && v_costs1.len() != v_costs2.len() {
        if config.verbosity > 1 {
            eprintln!("SPLIT_SKIP: cost vector length mismatch");
        }
        return None;
    }

    // Find position with minimum combined cost
    // Matches C++ AGC agc_compressor.cpp:1663-1674
    let mut best_pos = if let Some((p, _)) = maybe_best { p } else {
        let mut best_sum = u32::MAX;
        let mut pos = 0usize;
        for i in 0..v_costs1.len() {
            let cs = v_costs1[i].saturating_add(v_costs2[i]);
            if cs < best_sum {
                best_sum = cs;
                pos = i;
            }
        }
        pos
    };

    #[cfg(feature = "verbose_debug")]
    if std::env::var("RAGC_DEBUG_SPLIT_MAP").is_ok() && maybe_best.is_none() {
        let start = best_pos.saturating_sub(3);
        let end = (best_pos + 4).min(v_costs1.len());
        eprintln!("RAGC_COST_WINDOW: len={} best_pos={}", v_costs1.len(), best_pos);
        for i in start..end {
            eprintln!(
                "  i={} Lcum={} Rcum={} Sum={}{}",
                i,
                v_costs1[i],
                v_costs2[i],
                v_costs1[i].saturating_add(v_costs2[i]),
                if i == best_pos { "  <--" } else { "" }
            );
        }
    }

    // Apply degenerate position rules only when computing best_pos locally.
    // When FFI is used, these rules have already been applied in C++.
    if maybe_best.is_none() {
        let k = config.k;
        let original_best_pos = best_pos;  // Save for logging
        if best_pos < k + 1 {
            best_pos = 0; // Too close to start
        }
        if best_pos + k + 1 > v_costs1.len() {
            best_pos = v_costs1.len(); // Too close to end
        }

        if config.verbosity > 1 && original_best_pos != best_pos {
            eprintln!(
                "COST_CALC: original_best_pos={} forced_to={} (len={}, k+1={})",
                original_best_pos, best_pos, v_costs1.len(), k + 1
            );
        }
    }

    // Check if split is degenerate (C++ AGC agc_compressor.cpp:1400-1415)
    // C++ AGC ACCEPTS degenerate splits and assigns whole segment to one group
    // First compute sizes with exact best_pos; map to bytes afterward.
    let left_size_pre = best_pos;
    let right_size_pre = segment_data.len().saturating_sub(best_pos);

    if left_size_pre == 0 {
        // Degenerate: whole segment matches RIGHT group
        // Return empty left, full segment as right (C++ AGC line 1400-1407)
        if config.verbosity > 1 {
            eprintln!(
                "SPLIT_DEGENERATE_RIGHT: best_pos=0, assigning whole segment to RIGHT group"
            );
        }
        return Some((Vec::new(), segment_data.to_vec(), middle_kmer));
    }

    if right_size_pre == 0 {
        // Degenerate: whole segment matches LEFT group
        // Return full segment as left, empty right (C++ AGC line 1408-1415)
        if config.verbosity > 1 {
            eprintln!(
                "SPLIT_DEGENERATE_LEFT: best_pos=len, assigning whole segment to LEFT group"
            );
        }
        return Some((segment_data.to_vec(), Vec::new(), middle_kmer));
    }

    // Non-degenerate split: use FFI seg2_start directly (it already accounts for orientation)
    let (left_data, right_data) = if let Some((bp, s2)) = maybe_best {
        if config.verbosity > 1 {
            eprintln!(
                "SPLIT_GEOM_SELECT(FFI): best_pos={} seg2_start={} should_reverse={}",
                bp, s2, should_reverse
            );
        }
        split_segment_from_start(segment_data.as_slice(), s2, config.k)
    } else {
        let half = if should_reverse { (config.k + 1) / 2 } else { config.k / 2 };
        let seg2_start = best_pos.saturating_sub(half);
        if config.verbosity > 1 {
            eprintln!(
                "SPLIT_GEOM_SELECT(local): best_pos={} k={} half={} seg2_start={} should_reverse={}",
                best_pos, config.k, half, seg2_start, should_reverse
            );
        }
        split_segment_from_start(segment_data.as_slice(), seg2_start, config.k)
    };

    if config.verbosity > 1 {
        eprintln!(
            "SPLIT_SUCCESS: best_pos={} cost={} left_len={} right_len={}",
            best_pos,
            0u32, // best_sum not available under FFI path; placeholder
            left_data.len(),
            right_data.len()
        );
    }

    Some((left_data, right_data, middle_kmer))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_compressor() {
        let config = StreamingQueueConfig::default();
        let splitters = HashSet::new();
        let compressor =
            StreamingQueueCompressor::with_splitters("/tmp/test_stream.agc", config, splitters);
        assert!(compressor.is_ok());
    }

    #[test]
    fn test_queue_stats() {
        let config = StreamingQueueConfig::default();
        let splitters = HashSet::new();
        let compressor =
            StreamingQueueCompressor::with_splitters("/tmp/test_stats.agc", config, splitters)
                .unwrap();

        let stats = compressor.queue_stats();
        assert_eq!(stats.current_size_bytes, 0);
        assert_eq!(stats.current_items, 0);
        assert_eq!(stats.capacity_bytes, 2 * 1024 * 1024 * 1024);
        assert!(!stats.is_closed);
    }

    #[test]
    fn test_push_and_finalize() {
        let config = StreamingQueueConfig {
            verbosity: 0, // Quiet for tests
            ..Default::default()
        };
        let splitters = HashSet::new();
        let mut compressor =
            StreamingQueueCompressor::with_splitters("/tmp/test_push.agc", config, splitters)
                .unwrap();

        // Push a small contig
        let data = vec![b'A'; 1000];
        compressor
            .push("sample1".to_string(), "chr1".to_string(), data)
            .unwrap();

        // Finalize
        compressor.finalize().unwrap();
    }
}
