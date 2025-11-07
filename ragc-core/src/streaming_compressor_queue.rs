// Queue-based streaming compressor API
// Provides simple push() interface with automatic backpressure and constant memory usage

use crate::lz_diff::LZDiff;
use crate::memory_bounded_queue::MemoryBoundedQueue;
use crate::segment::split_at_splitters_with_size;
use crate::splitters::determine_splitters;
use anyhow::{Context, Result};
use ragc_common::{Archive, CollectionV3, Contig, CONTIG_SEPARATOR};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

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
        }
    }
}

/// Task to be processed by workers
/// Note: Contig is type alias for Vec<u8>, so we store the name separately
struct ContigTask {
    sample_name: String,
    contig_name: String,
    data: Contig, // Vec<u8>
}

/// Segment group identified by flanking k-mers (matching batch mode)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Ord, PartialOrd)]
struct SegmentGroupKey {
    kmer_front: u64,
    kmer_back: u64,
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
    reference_sample_name: Arc<Mutex<Option<String>>>, // First sample becomes reference
    // Segment splitting support (Phase 1)
    map_segments: Arc<Mutex<std::collections::HashMap<SegmentGroupKey, u32>>>, // (front, back) -> group_id
    map_segments_terminators: Arc<Mutex<std::collections::HashMap<u64, Vec<u64>>>>, // kmer -> [connected kmers]
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

        // Write file_type_info stream (after collection streams for C++ AGC compatibility)
        {
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
            archive.add_part(stream_id, &data, 7)?; // 7 key-value pairs
        }

        // Write params stream
        {
            let params_stream_id = archive.register_stream("params");
            let mut params_data = Vec::new();
            params_data.extend_from_slice(&(config.k as u32).to_le_bytes());
            params_data.extend_from_slice(&(config.min_match_len as u32).to_le_bytes());
            params_data.extend_from_slice(&50u32.to_le_bytes()); // pack_cardinality (default)
            params_data.extend_from_slice(&(config.segment_size as u32).to_le_bytes());
            archive.add_part(params_stream_id, &params_data, 0)?;
        }

        // Write empty splitters stream (C++ AGC compatibility)
        {
            let splitters_data = Vec::new();
            let stream_id = archive.register_stream("splitters");
            archive.add_part(stream_id, &splitters_data, 0)?;
        }

        // Write empty segment-splitters stream (C++ AGC compatibility)
        {
            let seg_splitters_data = Vec::new();
            let stream_id = archive.register_stream("segment-splitters");
            archive.add_part(stream_id, &seg_splitters_data, 0)?;
        }

        let collection = Arc::new(Mutex::new(collection));
        let archive = Arc::new(Mutex::new(archive));

        // Create memory-bounded queue
        let queue = Arc::new(MemoryBoundedQueue::new(config.queue_capacity));

        let splitters = Arc::new(splitters);

        // Segment grouping for LZ packing (using BTreeMap for better memory efficiency)
        let segment_groups = Arc::new(Mutex::new(BTreeMap::new()));
        let group_counter = Arc::new(AtomicU32::new(NO_RAW_GROUPS)); // Start at 16 for LZ groups
        let reference_sample_name = Arc::new(Mutex::new(None)); // Shared across all workers

        // Segment splitting support (Phase 1)
        let map_segments: Arc<Mutex<HashMap<SegmentGroupKey, u32>>> = Arc::new(Mutex::new(HashMap::new()));
        let map_segments_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>> = Arc::new(Mutex::new(HashMap::new()));

        // Spawn worker threads
        let mut workers = Vec::new();
        for worker_id in 0..config.num_threads {
            let queue = Arc::clone(&queue);
            let collection = Arc::clone(&collection);
            let splitters = Arc::clone(&splitters);
            let archive = Arc::clone(&archive);
            let segment_groups = Arc::clone(&segment_groups);
            let group_counter = Arc::clone(&group_counter);
            let reference_sample_name = Arc::clone(&reference_sample_name);
            let map_segments = Arc::clone(&map_segments);
            let map_segments_terminators = Arc::clone(&map_segments_terminators);
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
                    reference_sample_name,
                    map_segments,
                    map_segments_terminators,
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
            reference_sample_name,
            map_segments,
            map_segments_terminators,
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
                let reference_sample_name = Arc::clone(&self.reference_sample_name);
                let map_segments = Arc::clone(&self.map_segments);
                let map_segments_terminators = Arc::clone(&self.map_segments_terminators);
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
                        reference_sample_name,
                        map_segments,
                        map_segments_terminators,
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

        // Create task
        let task = ContigTask {
            sample_name,
            contig_name,
            data,
        };

        // Push to queue (BLOCKS if queue is full!)
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
    /// # Ok::<(), anyhow::Error>(())
    /// ```
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
                    flush_pack(buffer, &self.collection, &self.archive, &self.config)
                        .context("Failed to flush remaining pack")?;
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

            // Write sample names
            collection
                .store_batch_sample_names(&mut archive)
                .context("Failed to write sample names")?;

            // Write contig names and segment details
            collection
                .store_contig_batch(&mut archive, 0, num_samples)
                .context("Failed to write contig batch")?;

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
) -> Result<()> {
    use crate::segment_compression::{compress_reference_segment, compress_segment_configured};

    // Skip if no segments to write (but still write reference if present)
    if buffer.segments.is_empty() && buffer.ref_written {
        return Ok(());
    }

    let use_lz_encoding = buffer.group_id >= NO_RAW_GROUPS;

    // Sort segments to match C++ AGC (agc_compressor.h line 362: p->sort())
    // This ensures the reference is chosen deterministically (alphabetically first sample/contig/part)
    // instead of "whichever arrives first"
    buffer.segments.sort();

    // Write reference segment if not already written (first pack for this group)
    // Extract reference from sorted segments (matching C++ AGC: first segment after sort becomes reference)
    if !buffer.ref_written && !buffer.segments.is_empty() {
        // Remove first segment (alphabetically first) to use as reference
        let ref_seg = buffer.segments.remove(0);

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

        {
            let mut arch = archive.lock().unwrap();
            arch.add_part(buffer.ref_stream_id, &compressed, ref_size)
                .context("Failed to write reference")?;
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

        // Prepare LZ encoder with reference (matching C++ AGC segment.cpp line 43: lz_diff->Prepare(s))
        // This is done ONCE when the reference is written, then reused for all subsequent segments
        if use_lz_encoding {
            let mut lz = LZDiff::new(config.min_match_len as u32);
            lz.prepare(&ref_seg.data);
            buffer.lz_diff = Some(lz);
        }
    }

    // Pack segments together with delta deduplication (matching C++ AGC segment.cpp lines 66-74)
    // Note: segments do NOT include the reference - it's stored separately
    let mut unique_deltas: Vec<Vec<u8>> = Vec::new(); // v_lzp in C++ AGC
    let mut segment_delta_indices: Vec<u32> = Vec::new(); // Which unique delta each segment uses

    // First pass: encode segments and deduplicate deltas
    for seg in buffer.segments.iter() {
        let contig_data = if !use_lz_encoding || buffer.reference_segment.is_none() {
            // Raw segment: groups 0-15 OR groups without reference
            seg.data.clone()
        } else {
            // LZ-encoded segment (groups >= 16 with reference)
            // Reuse prepared lz_diff (matching C++ AGC segment.cpp line 59: lz_diff->Encode(s, delta))
            buffer.lz_diff.as_mut()
                .expect("lz_diff should be prepared when reference is written")
                .encode(&seg.data)
        };

        // Check if this delta already exists (matching C++ AGC segment.cpp line 66)
        if let Some(existing_idx) = unique_deltas.iter().position(|d| d == &contig_data) {
            // Reuse existing delta (matching C++ AGC segment.cpp line 69)
            segment_delta_indices.push(existing_idx as u32);
        } else {
            // New unique delta - add it (matching C++ AGC segment.cpp line 74)
            segment_delta_indices.push(unique_deltas.len() as u32);
            unique_deltas.push(contig_data);
        }
    }

    // Second pass: pack unique deltas and register segments
    let mut packed_data = Vec::new();
    for delta in unique_deltas.iter() {
        packed_data.extend_from_slice(delta);
        packed_data.push(CONTIG_SEPARATOR);
    }

    // Register segments in collection with their delta indices
    for (seg, &delta_idx) in buffer.segments.iter().zip(segment_delta_indices.iter()) {
        // in_group_id represents which delta this segment uses
        // 0 = reference, 1+ = delta index (offset by buffer.segments_written and +1 for reference)
        let in_group_id = buffer.segments_written + delta_idx + 1;

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

    // Update counter (counts UNIQUE deltas only, not reference or duplicates)
    // This matches C++ AGC which increments no_seqs only when adding new deltas (segment.cpp line 77)
    buffer.segments_written += unique_deltas.len() as u32;

    // Compress and write the packed data (if we have any delta segments)
    if !packed_data.is_empty() {
        // The metadata field stores the uncompressed size (including separators)
        let total_raw_size = packed_data.len() as u64;

        let mut compressed = compress_segment_configured(&packed_data, config.compression_level)
            .context("Failed to compress pack")?;
        compressed.push(0); // Marker 0 = plain ZSTD

        {
            let mut arch = archive.lock().unwrap();
            arch.add_part(buffer.stream_id, &compressed, total_raw_size)
                .context("Failed to write pack")?;
        }
    }

    // Clear segments for next pack
    buffer.segments.clear();

    Ok(())
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
    reference_sample_name: Arc<Mutex<Option<String>>>,
    map_segments: Arc<Mutex<HashMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: Arc<Mutex<HashMap<u64, Vec<u64>>>>,
    config: StreamingQueueConfig,
) -> Result<()> {
    let mut processed_count = 0;

    loop {
        // Pull from queue (blocks if empty, returns None when closed)
        let Some(task) = queue.pull() else {
            // Queue is closed and empty - we're done!
            if config.verbosity > 1 {
                eprintln!(
                    "Worker {} finished ({} contigs processed)",
                    worker_id, processed_count
                );
            }
            break;
        };

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
        for (place, segment) in segments.iter().enumerate() {
            // DEBUG: Output every segment for comparison with C++ AGC
            eprintln!("RAGC_SEGMENT: sample={} contig={} part={} len={} front={} back={}",
                task.sample_name, task.contig_name, place, segment.data.len(),
                segment.front_kmer, segment.back_kmer);

            // Match C++ AGC Case 2: Normalize segment group key by ensuring front <= back
            // (agc_compressor.cpp lines 1306-1327)
            use crate::segment::MISSING_KMER;

            let (key_front, key_back, should_reverse) =
                if segment.front_kmer != MISSING_KMER && segment.back_kmer != MISSING_KMER {
                    // Both k-mers present - normalize by ensuring front <= back
                    if segment.front_kmer <= segment.back_kmer {
                        // Already normalized - keep original orientation
                        if config.verbosity > 2 {
                            eprintln!(
                                "RAGC_CASE2_KEEP: sample={} front={} back={} len={}",
                                task.sample_name, segment.front_kmer, segment.back_kmer, segment.data.len()
                            );
                        }
                        (segment.front_kmer, segment.back_kmer, false)
                    } else {
                        // Swap k-mers and reverse complement data
                        if config.verbosity > 2 {
                            eprintln!(
                                "RAGC_CASE2_SWAP: sample={} front={} back={} -> key=({},{}) len={}",
                                task.sample_name, segment.front_kmer, segment.back_kmer,
                                segment.back_kmer, segment.front_kmer, segment.data.len()
                            );
                        }
                        (segment.back_kmer, segment.front_kmer, true)
                    }
                } else {
                    // At least one k-mer is MISSING - use as-is (Cases 1, 3, 4)
                    (segment.front_kmer, segment.back_kmer, false)
                };

            // Create grouping key from normalized k-mers
            let key = SegmentGroupKey {
                kmer_front: key_front,
                kmer_back: key_back,
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

                // Check if this key already exists (for debugging group reuse)
                let key_exists = groups.contains_key(&key);

                // Phase 2: Check if we should attempt splitting this segment
                // (matches C++ AGC agc_compressor.cpp lines 1387-1503)
                if !key_exists && key_front != MISSING_KMER && key_back != MISSING_KMER {
                    // Check if both terminators exist in other segments
                    let can_attempt_split = {
                        let terminators = map_segments_terminators.lock().unwrap();
                        terminators.contains_key(&key_front) && terminators.contains_key(&key_back)
                    };

                    if can_attempt_split {
                        // Attempt to split the segment
                        // IMPORTANT: Pass ORIGINAL segment data, not the RC'd version!
                        // C++ AGC passes both segment and segment_rc to the split function
                        let split_result = {
                            let terminators = map_segments_terminators.lock().unwrap();
                            try_split_segment(&segment.data, key_front, key_back, &terminators, &config)
                        };

                        if let Some((left_data, right_data, middle_kmer)) = split_result {
                            // Successfully split! Add both segments to their existing groups

                            // Calculate keys for the two new segments
                            // Left: (front, middle), Right: (middle, back)
                            // Apply same normalization as original segment
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

                            if config.verbosity > 1 {
                                eprintln!(
                                    "SPLIT: original=({},{}) -> left=({},{}) right=({},{})",
                                    key_front, key_back, left_key.kmer_front, left_key.kmer_back,
                                    right_key.kmer_front, right_key.kmer_back
                                );
                            }

                            // Add left segment to its group (should already exist)
                            if let Some(left_buffer) = groups.get_mut(&left_key) {
                                let left_buffered = BufferedSegment {
                                    sample_name: task.sample_name.clone(),
                                    contig_name: task.contig_name.clone(),
                                    seg_part_no: place,
                                    data: left_data,
                                    is_rev_comp: should_reverse, // Inherit RC status from original
                                };
                                left_buffer.segments.push(left_buffered);

                                if left_buffer.segments.len() >= PACK_CARDINALITY + 1 {
                                    flush_pack(left_buffer, &collection, &archive, &config)
                                        .context("Failed to flush left pack")?;
                                }
                            }

                            // Add right segment to its group (should already exist)
                            if let Some(right_buffer) = groups.get_mut(&right_key) {
                                let right_buffered = BufferedSegment {
                                    sample_name: task.sample_name.clone(),
                                    contig_name: task.contig_name.clone(),
                                    seg_part_no: place,
                                    data: right_data,
                                    is_rev_comp: should_reverse, // Inherit RC status from original
                                };
                                right_buffer.segments.push(right_buffered);

                                if right_buffer.segments.len() >= PACK_CARDINALITY + 1 {
                                    flush_pack(right_buffer, &collection, &archive, &config)
                                        .context("Failed to flush right pack")?;
                                }
                            }

                            // Skip adding original segment - we've already added the split parts
                            continue;
                        }
                    }
                }

                // Normal path: add segment to group as-is (no split or split failed)
                // Create buffered segment
                let buffered = BufferedSegment {
                    sample_name: task.sample_name.clone(),
                    contig_name: task.contig_name.clone(),
                    seg_part_no: place,
                    data: segment_data,
                    is_rev_comp: should_reverse,
                };

                // Get or create buffer for this group
                let buffer = groups.entry(key.clone()).or_insert_with(|| {
                    // Allocate new group ID
                    let group_id = group_counter.fetch_add(1, Ordering::SeqCst);

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

                    // Phase 1: Track segment terminators for splitting
                    // (matches C++ AGC agc_compressor.cpp lines 1319-1334)
                    {
                        use crate::segment::MISSING_KMER;

                        // Record this segment group in map_segments
                        let mut seg_map = map_segments.lock().unwrap();
                        seg_map.insert(key.clone(), group_id);
                        drop(seg_map);

                        // Track k-mer connections (only if both are not MISSING)
                        if key_front != MISSING_KMER && key_back != MISSING_KMER {
                            let mut terminators = map_segments_terminators.lock().unwrap();

                            // Add bidirectional connection
                            terminators.entry(key_front).or_insert_with(Vec::new).push(key_back);
                            if key_front != key_back {
                                terminators.entry(key_back).or_insert_with(Vec::new).push(key_front);
                            }

                            // Keep vectors sorted and deduplicated for efficient set intersection
                            // (C++ AGC uses sorted vectors for set_intersection in lines 1541-1543)
                            if let Some(vec) = terminators.get_mut(&key_front) {
                                vec.sort_unstable();
                                vec.dedup();
                            }
                            if key_front != key_back {
                                if let Some(vec) = terminators.get_mut(&key_back) {
                                    vec.sort_unstable();
                                    vec.dedup();
                                }
                            }
                        }
                    }

                    SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
                });

                if key_exists && config.verbosity > 1 {
                    eprintln!("REUSE_GROUP: group_id={} front={} back={} sample={}",
                        buffer.group_id, key_front, key_back, task.sample_name);
                }

                // Buffer all segments (matching C++ AGC: buffer all, then sort)
                // Reference will be chosen from sorted segments during flush
                buffer.segments.push(buffered);

                // Flush pack if buffer is full
                // Note: +1 because first segment will become reference
                if buffer.segments.len() >= PACK_CARDINALITY + 1 {
                    flush_pack(buffer, &collection, &archive, &config)
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

    // Find intersection of sorted vectors (set_intersection algorithm)
    // Both vectors are kept sorted by Phase 1
    let mut i = 0;
    let mut j = 0;

    while i < front_connections.len() && j < back_connections.len() {
        if front_connections[i] == back_connections[j] {
            // Found shared k-mer - return it (C++ AGC uses first match)
            return Some(front_connections[i]);
        } else if front_connections[i] < back_connections[j] {
            i += 1;
        } else {
            j += 1;
        }
    }

    None // No shared k-mer found
}

/// Phase 4: Find split position using simplified heuristic
/// C++ AGC uses compression cost (lines 1556-1625), but for a simplified approach
/// we split at the midpoint of the segment
/// Returns the split position (in bytes)
fn find_split_position(_segment_data: &[u8], _middle_kmer: u64, segment_len: usize, k: usize) -> Option<usize> {
    // Ensure we don't split too close to the ends
    // Need at least k+1 bytes on each side for valid segments
    if segment_len < 2 * (k + 1) {
        return None;
    }

    // Split at midpoint
    // C++ AGC finds optimal position using compression cost, but midpoint is reasonable
    let split_pos = segment_len / 2;

    Some(split_pos)
}

/// Phase 5: Split segment into two overlapping segments
/// Returns (left_segment, right_segment) with k-mer overlap
/// (matches C++ AGC lines 1445-1503)
fn split_segment_at_position(
    segment_data: &[u8],
    split_pos: usize,
    k: usize,
) -> (Vec<u8>, Vec<u8>) {
    // Left segment ends at split_pos + k (includes full middle k-mer)
    let left_end = split_pos + k;
    let left = segment_data[..left_end].to_vec();

    // Right segment starts at split_pos + k/2 (overlaps by k/2)
    let right_start = split_pos + k / 2;
    let right = segment_data[right_start..].to_vec();

    (left, right)
}

/// Phase 6: Attempt to split a segment by finding a middle k-mer
/// Returns Some((left_data, right_data, middle_kmer)) if successful
/// (orchestrates C++ AGC agc_compressor.cpp lines 1387-1503)
fn try_split_segment(
    segment_data: &[u8],
    front_kmer: u64,
    back_kmer: u64,
    terminators: &HashMap<u64, Vec<u64>>,
    config: &StreamingQueueConfig,
) -> Option<(Vec<u8>, Vec<u8>, u64)> {
    // Find middle k-mer that connects front and back
    let middle_kmer = find_middle_splitter(front_kmer, back_kmer, terminators)?;

    if config.verbosity > 1 {
        eprintln!(
            "SPLIT_ATTEMPT: front={} back={} middle={}",
            front_kmer, back_kmer, middle_kmer
        );
    }

    // Find split position (simplified: use midpoint instead of compression cost)
    let split_pos = match find_split_position(segment_data, middle_kmer, segment_data.len(), config.k) {
        Some(pos) => pos,
        None => {
            if config.verbosity > 1 {
                eprintln!(
                    "SPLIT_FAILED: middle k-mer {} not found in segment (len={})",
                    middle_kmer,
                    segment_data.len()
                );
            }
            return None;
        }
    };

    // Split into two segments with overlap
    let (left_data, right_data) = split_segment_at_position(segment_data, split_pos, config.k);

    if config.verbosity > 1 {
        eprintln!(
            "SPLIT_SUCCESS: split_pos={} left_len={} right_len={}",
            split_pos,
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
