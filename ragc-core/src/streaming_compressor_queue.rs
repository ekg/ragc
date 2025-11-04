// Queue-based streaming compressor API
// Provides simple push() interface with automatic backpressure and constant memory usage

use crate::memory_bounded_queue::MemoryBoundedQueue;
use crate::segment::split_at_splitters_with_size;
use crate::segment_compression::{compress_reference_segment, compress_segment_configured};
use crate::splitters::determine_splitters;
use anyhow::{Context, Result};
use ragc_common::{Archive, CollectionV3, Contig};
use std::collections::HashSet;
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
    segment_counter: Arc<AtomicU32>,
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
            append_str(&mut data, &format!("RAGC v.{}.{}", ragc_common::AGC_FILE_MAJOR, ragc_common::AGC_FILE_MINOR));

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
        // Start at 16 for LZ groups (< 16 are raw groups, >= 16 are LZ groups with ref streams)
        // Each segment becomes its own reference (no LZ differential encoding yet)
        let segment_counter = Arc::new(AtomicU32::new(16));

        // Spawn worker threads
        let mut workers = Vec::new();
        for worker_id in 0..config.num_threads {
            let queue = Arc::clone(&queue);
            let collection = Arc::clone(&collection);
            let splitters = Arc::clone(&splitters);
            let archive = Arc::clone(&archive);
            let segment_counter = Arc::clone(&segment_counter);
            let config = config.clone();

            let handle = thread::spawn(move || {
                worker_thread(worker_id, queue, collection, splitters, archive, segment_counter, config)
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
            segment_counter,
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
    pub fn push(
        &mut self,
        sample_name: String,
        contig_name: String,
        data: Contig,
    ) -> Result<()> {
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
                let segment_counter = Arc::clone(&self.segment_counter);
                let config = self.config.clone();

                let handle = thread::spawn(move || {
                    worker_thread(worker_id, queue, collection, splitters, archive, segment_counter, config)
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
    pub fn finalize(mut self) -> Result<()> {
        if self.config.verbosity > 0 {
            eprintln!("Finalizing compression...");
            eprintln!("  Closing queue...");
        }

        // Close queue - no more pushes allowed
        self.queue.close();

        if self.config.verbosity > 0 {
            eprintln!(
                "  Waiting for {} workers to finish...",
                self.workers.len()
            );
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
            collection.store_batch_sample_names(&mut archive)
                .context("Failed to write sample names")?;

            // Write contig names and segment details
            collection.store_contig_batch(&mut archive, 0, num_samples)
                .context("Failed to write contig batch")?;

            if self.config.verbosity > 0 {
                eprintln!("Collection metadata written successfully");
            }

            // Close archive (writes footer)
            archive.close()
                .context("Failed to close archive")?;
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

/// Worker thread that pulls from queue and compresses
fn worker_thread(
    worker_id: usize,
    queue: Arc<MemoryBoundedQueue<ContigTask>>,
    collection: Arc<Mutex<CollectionV3>>,
    splitters: Arc<HashSet<u64>>,
    archive: Arc<Mutex<Archive>>,
    segment_counter: Arc<AtomicU32>,
    config: StreamingQueueConfig,
) -> Result<()> {
    use crate::segment_compression::compress_reference_segment;

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
        let segments = split_at_splitters_with_size(&task.data, &splitters, config.k, config.segment_size);

        if config.verbosity > 2 {
            eprintln!(
                "Worker {} processing {} (split into {} segments)",
                worker_id,
                task.contig_name,
                segments.len()
            );
        }

        // Register and compress each segment
        for (place, segment) in segments.iter().enumerate() {
            // Get a global segment number for this segment
            let seg_num = segment_counter.fetch_add(1, Ordering::SeqCst);

            // Register this segment in the collection
            {
                let mut coll = collection.lock().unwrap();
                coll.add_segment_placed(
                    &task.sample_name,
                    &task.contig_name,
                    place,
                    seg_num,  // group_id (each segment is its own group for now)
                    0,         // in_group_id (0 = reference)
                    false,     // is_rev_comp
                    segment.data.len() as u32,
                )
                .context("Failed to add segment to collection")?;
            }

            // Write segment to archive
            // Store as compressed reference (matching batch mode)
            // TODO: Add LZ differential encoding for subsequent segments in same group
            {
                // Compress segment using adaptive compression (plain ZSTD or tuple packing)
                let (mut compressed, marker) = compress_reference_segment(&segment.data)
                    .context("Failed to compress segment")?;

                // Append marker byte to compressed data
                // Decompressor pops this byte when metadata != 0
                compressed.push(marker);

                let mut arch = archive.lock().unwrap();

                // Register ref stream for this segment (LZ groups use ref streams for references)
                let stream_name = ragc_common::stream_ref_name(3000, seg_num);
                let stream_id = arch.register_stream(&stream_name);

                // Write compressed data with metadata=1 (indicates compressed + marker format)
                arch.add_part(stream_id, &compressed, 1)
                    .context("Failed to write segment to archive")?;
            }
        }

        processed_count += 1;
    }

    Ok(())
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
