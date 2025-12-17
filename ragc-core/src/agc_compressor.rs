// Queue-based streaming compressor API
// Provides simple push() interface with automatic backpressure and constant memory usage

use crate::kmer_extract::{enumerate_kmers, remove_non_singletons};
use crate::lz_diff::LZDiff;
use crate::memory_bounded_queue::MemoryBoundedQueue;
use crate::segment::{split_at_splitters_with_size, MISSING_KMER};
use crate::splitters::{determine_splitters, find_new_splitters_for_contig};
use anyhow::{Context, Result};
use ragc_common::{Archive, CollectionV3, Contig, CONTIG_SEPARATOR};
use ahash::AHashSet;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;
use std::sync::atomic::{AtomicI32, AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock};
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

    /// Batch size: number of samples to accumulate before sorting and distributing
    /// (matches C++ AGC pack_cardinality parameter, default 50)
    /// Segments from batch_size samples are sorted by (sample, contig, seg_part_no)
    /// before distribution to groups, ensuring consistent pack boundaries with C++ AGC.
    pub batch_size: usize,

    /// Pack size: number of segments per pack (matches C++ AGC contigs_in_pack)
    /// When a group reaches this many segments, write a pack immediately
    /// (default: 50, matching PACK_CARDINALITY)
    pub pack_size: usize,

    /// Concatenated genomes mode: if true, send sync tokens every pack_size contigs
    /// If false (multiple input files), only send sync tokens at sample boundaries
    /// Matches C++ AGC's concatenated_genomes behavior
    pub concatenated_genomes: bool,
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
            batch_size: 50,       // Default matches C++ AGC pack_cardinality
            pack_size: 50,        // Default matches C++ AGC contigs_in_pack / PACK_CARDINALITY
            concatenated_genomes: false, // Default: multiple input files (non-concatenated)
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
    cost: usize, // Contig size in bytes (matches C++ AGC cost calculation)
    sequence: u64, // Insertion order within sample - lower = processed first (FASTA order)
    is_sync_token: bool, // True if this is a synchronization token (matches C++ AGC registration tokens)
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
        self.sample_priority == other.sample_priority
            && self.cost == other.cost
            && self.contig_name == other.contig_name
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
        // C++ AGC uses (priority, cost) as the multimap key with PopLarge (rbegin).
        // multimap is sorted by (priority, cost) in ASCENDING order.
        // rbegin() returns LARGEST element, so within same priority, LARGEST cost is popped first.
        //
        // Example: AAA#0 contigs are processed by SIZE (largest first):
        //   chrIV (1.5MB) → chrXV (1.1MB) → chrVII (1.1MB) → ... → chrMT (86KB)
        //
        // This is NOT file order! Instrumentation shows C++ AGC pops by (priority, cost).

        // First compare by sample_priority (higher priority first)
        match self.sample_priority.cmp(&other.sample_priority) {
            std::cmp::Ordering::Equal => {
                // Then by cost (LARGER cost = higher priority, processed first)
                // Match C++ AGC's PopLarge behavior
                match self.cost.cmp(&other.cost) {
                    std::cmp::Ordering::Equal => {
                        // CRITICAL TIE-BREAKER: When sizes are equal, use FASTA order (sequence field)
                        // to ensure deterministic ordering. Without this, the BinaryHeap order is
                        // non-deterministic, causing different segment splitting and 19% size difference.
                        // LOWER sequence = earlier in FASTA = processed first (reverse comparison for max-heap)
                        other.sequence.cmp(&self.sequence)
                    }
                    cost_ord => cost_ord,
                }
            }
            priority_ord => priority_ord,
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
/// Segments are sorted by (sample_name, contig_name, place) to match C++ AGC order
#[derive(Debug, Clone, PartialEq, Eq)]
struct PendingSegment {
    key: SegmentGroupKey,
    segment_data: Vec<u8>,
    should_reverse: bool,
    sample_name: String,
    contig_name: String,
    place: usize,
    sample_priority: i32, // Sample processing order (higher = earlier)
}

// Match C++ AGC sorting order (agc_compressor.h lines 112-119)
// Sort by: sample_name, then contig_name, then place (seg_part_no)
impl PartialOrd for PendingSegment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PendingSegment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Match C++ AGC: pure lexicographic ordering (no sample_priority)
        // Sort by: sample_name, then contig_name, then place (seg_part_no)
        match self.sample_name.cmp(&other.sample_name) {
            std::cmp::Ordering::Equal => {
                // Then by contig_name
                match self.contig_name.cmp(&other.contig_name) {
                    std::cmp::Ordering::Equal => {
                        // Finally by place (seg_part_no)
                        self.place.cmp(&other.place)
                    }
                    other => other,
                }
            }
            other => other,
        }
    }
}

/// Buffered segment waiting to be packed
#[derive(Debug, Clone, PartialEq, Eq)]
struct BufferedSegment {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    data: Contig,
    is_rev_comp: bool,
    sample_priority: i32, // Sample processing order (higher = earlier)
}

// Match C++ AGC sorting order: pure lexicographic (no sample_priority)
// Sort by: sample_name, contig_name, seg_part_no
impl PartialOrd for BufferedSegment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BufferedSegment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Match C++ AGC: pure lexicographic ordering (no sample_priority)
        // Sort by: sample_name, then contig_name, then seg_part_no
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

// =============================================================================
// RAW segment buffering for parallel Phase 1 (BEFORE classification)
// =============================================================================

/// Raw segment data buffered BEFORE k-mer classification.
/// This allows parallel buffering without lock contention from find_group_with_one_kmer.
/// Classification is deferred to Thread 0 at the barrier.
#[derive(Clone)]
struct RawBufferedSegment {
    /// Raw segment data (numeric encoding: 0=A, 1=C, 2=G, 3=T)
    data: Vec<u8>,
    /// Precomputed reverse complement of data
    data_rc: Vec<u8>,
    /// Front k-mer from segment detection
    front_kmer: u64,
    /// Back k-mer from segment detection
    back_kmer: u64,
    /// Is front k-mer in canonical direction?
    front_kmer_is_dir: bool,
    /// Is back k-mer in canonical direction?
    back_kmer_is_dir: bool,
    /// Sample name for sorting and registration
    sample_name: String,
    /// Contig name for sorting and registration
    contig_name: String,
    /// Segment index within contig (before split adjustment)
    original_place: usize,
    /// Sample processing priority (higher = earlier)
    sample_priority: i32,
}

// Implement Ord for deterministic sorting at barrier
impl PartialEq for RawBufferedSegment {
    fn eq(&self, other: &Self) -> bool {
        self.sample_name == other.sample_name
            && self.contig_name == other.contig_name
            && self.original_place == other.original_place
    }
}
impl Eq for RawBufferedSegment {}

impl PartialOrd for RawBufferedSegment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RawBufferedSegment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Sort by: sample_name, contig_name, original_place (matches C++ AGC order)
        match self.sample_name.cmp(&other.sample_name) {
            std::cmp::Ordering::Equal => {
                match self.contig_name.cmp(&other.contig_name) {
                    std::cmp::Ordering::Equal => self.original_place.cmp(&other.original_place),
                    other => other,
                }
            }
            other => other,
        }
    }
}

// =============================================================================
// C++ AGC-style segment buffering for parallel compression (4-phase pattern)
// =============================================================================

/// Per-group segment buffer with its own mutex (C++ AGC: list_seg_part_t)
/// Each group has independent locking to allow parallel writes during Phase 1
struct PerGroupSegments {
    segments: Vec<BufferedSegment>,
}

/// Segment waiting to be assigned a group ID (C++ AGC: kk_seg_part_t)
/// Used during Phase 1 when segment's k-mer pair doesn't exist in map_segments yet
#[derive(Clone)]
struct NewSegment {
    /// K-mer pair (normalized: front <= back)
    kmer_front: u64,
    kmer_back: u64,
    /// Sort key for deterministic processing: (sample_priority, sample_name, contig_name, seg_part_no)
    sample_priority: i32,
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    /// Segment data
    data: Contig,
    should_reverse: bool,
}

// Implement Ord for NewSegment to match C++ AGC BTreeSet ordering
impl PartialEq for NewSegment {
    fn eq(&self, other: &Self) -> bool {
        self.sample_priority == other.sample_priority
            && self.sample_name == other.sample_name
            && self.contig_name == other.contig_name
            && self.seg_part_no == other.seg_part_no
    }
}
impl Eq for NewSegment {}

impl PartialOrd for NewSegment {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NewSegment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Higher sample_priority processes first (descending)
        match other.sample_priority.cmp(&self.sample_priority) {
            std::cmp::Ordering::Equal => {
                // Then by sample_name, contig_name, seg_part_no (ascending)
                match self.sample_name.cmp(&other.sample_name) {
                    std::cmp::Ordering::Equal => {
                        match self.contig_name.cmp(&other.contig_name) {
                            std::cmp::Ordering::Equal => self.seg_part_no.cmp(&other.seg_part_no),
                            other => other,
                        }
                    }
                    other => other,
                }
            }
            other => other,
        }
    }
}

/// Two-tier segment buffering for C++ AGC 4-phase pattern (C++ AGC: CBufferedSegPart)
///
/// Phase 1 (PARALLEL): Workers add segments using add_known() or add_new()
/// Phase 2 (SINGLE): Thread 0 calls process_new() to assign group IDs
/// Phase 3 (PARALLEL): Workers call get_vec_id() + get_part() for atomic work-stealing
/// Phase 4: Thread 0 calls clear() for cleanup
struct BufferedSegPart {
    /// KNOWN segments: indexed by group_id, each has own mutex
    /// RwLock on Vec allows process_new() to resize while add_known() reads
    /// C++ AGC: vector<list_seg_part_t> vl_seg_part
    vl_seg_part: RwLock<Vec<Mutex<PerGroupSegments>>>,

    /// NEW segments: BTreeSet for deterministic iteration
    /// C++ AGC: set<kk_seg_part_t> s_seg_part
    s_seg_part: Mutex<std::collections::BTreeSet<NewSegment>>,

    /// Atomic counter for work distribution (descending from num_groups-1 to -1)
    /// C++ AGC: atomic<int32_t> a_v_part_id
    a_v_part_id: AtomicI32,
}

impl BufferedSegPart {
    fn new(initial_groups: usize) -> Self {
        Self {
            vl_seg_part: RwLock::new(
                (0..initial_groups)
                    .map(|_| Mutex::new(PerGroupSegments { segments: Vec::new() }))
                    .collect()
            ),
            s_seg_part: Mutex::new(std::collections::BTreeSet::new()),
            a_v_part_id: AtomicI32::new(-1),
        }
    }

    /// Add segment to KNOWN group (has group_id)
    /// C++ AGC: add_known() - read lock on Vec, per-group lock on Mutex
    fn add_known(&self, group_id: u32, segment: BufferedSegment) {
        let groups = self.vl_seg_part.read().unwrap();
        if (group_id as usize) < groups.len() {
            groups[group_id as usize]
                .lock()
                .unwrap()
                .segments
                .push(segment);
        }
    }

    /// Add segment with UNKNOWN group (new k-mer pair)
    /// C++ AGC: add_new() - global s_seg_part lock (but brief)
    fn add_new(&self, segment: NewSegment) {
        self.s_seg_part.lock().unwrap().insert(segment);
    }

    /// Process NEW segments, assign group IDs deterministically
    /// C++ AGC: process_new() - ONLY called by thread 0 after barrier
    fn process_new(
        &self,
        map_segments: &mut BTreeMap<SegmentGroupKey, u32>,
        next_group_id: &mut u32,
    ) -> u32 {
        let mut s = self.s_seg_part.lock().unwrap();
        let mut m_kmers: BTreeMap<(u64, u64), u32> = BTreeMap::new();
        let mut new_count = 0u32;

        // First pass: assign group IDs (deterministic - BTreeSet order)
        for seg in s.iter() {
            let key = (seg.kmer_front, seg.kmer_back);
            if !m_kmers.contains_key(&key) && !map_segments.contains_key(&SegmentGroupKey {
                kmer_front: seg.kmer_front,
                kmer_back: seg.kmer_back,
            }) {
                m_kmers.insert(key, *next_group_id);
                *next_group_id += 1;
                new_count += 1;
            }
        }

        // Resize vl_seg_part for new groups (requires write lock)
        {
            let mut groups = self.vl_seg_part.write().unwrap();
            while groups.len() < *next_group_id as usize {
                groups.push(Mutex::new(PerGroupSegments {
                    segments: Vec::new(),
                }));
            }
        }

        // Second pass: move segments to vl_seg_part and update map_segments
        let segments: Vec<NewSegment> = s.iter().cloned().collect();
        s.clear();
        drop(s);

        for seg in segments {
            let key = SegmentGroupKey {
                kmer_front: seg.kmer_front,
                kmer_back: seg.kmer_back,
            };

            // Get group_id from either existing map or newly assigned
            let group_id = if let Some(&id) = map_segments.get(&key) {
                id
            } else if let Some(&id) = m_kmers.get(&(seg.kmer_front, seg.kmer_back)) {
                // Insert into map_segments
                map_segments.insert(key, id);
                id
            } else {
                continue; // Should not happen
            };

            // Add to per-group buffer (uses read lock internally)
            let buffered = BufferedSegment {
                sample_name: seg.sample_name,
                contig_name: seg.contig_name,
                seg_part_no: seg.seg_part_no,
                data: seg.data,
                is_rev_comp: seg.should_reverse,
                sample_priority: seg.sample_priority,
            };
            self.add_known(group_id, buffered);
        }

        new_count
    }

    /// Sort known segments within each group for deterministic output
    fn sort_known(&self) {
        let groups = self.vl_seg_part.read().unwrap();
        for group in groups.iter() {
            group.lock().unwrap().segments.sort();
        }
    }

    /// Reset atomic counter for work distribution
    /// C++ AGC: restart_read_vec()
    fn restart_read_vec(&self) {
        let groups = self.vl_seg_part.read().unwrap();
        self.a_v_part_id.store(
            groups.len() as i32 - 1,
            Ordering::SeqCst,
        );
    }

    /// Get next group_id to process (atomic decrement for work-stealing)
    /// C++ AGC: get_vec_id() - returns -1 when all groups claimed
    fn get_vec_id(&self) -> i32 {
        self.a_v_part_id.fetch_sub(1, Ordering::Relaxed)
    }

    /// Get and remove one segment from group (for store phase)
    /// C++ AGC: get_part()
    fn get_part(&self, group_id: u32) -> Option<BufferedSegment> {
        let groups = self.vl_seg_part.read().unwrap();
        if (group_id as usize) < groups.len() {
            groups[group_id as usize]
                .lock()
                .unwrap()
                .segments
                .pop()
        } else {
            None
        }
    }

    /// Get all segments from a group (for batch processing)
    fn drain_group(&self, group_id: u32) -> Vec<BufferedSegment> {
        let groups = self.vl_seg_part.read().unwrap();
        if (group_id as usize) < groups.len() {
            std::mem::take(&mut groups[group_id as usize].lock().unwrap().segments)
        } else {
            Vec::new()
        }
    }

    /// Clear all buffers after batch
    /// C++ AGC: clear()
    fn clear(&self) {
        let groups = self.vl_seg_part.read().unwrap();
        for group in groups.iter() {
            group.lock().unwrap().segments.clear();
        }
        self.s_seg_part.lock().unwrap().clear();
    }

    /// Check if any segments are buffered
    fn has_segments(&self) -> bool {
        let groups = self.vl_seg_part.read().unwrap();
        for group in groups.iter() {
            if !group.lock().unwrap().segments.is_empty() {
                return true;
            }
        }
        !self.s_seg_part.lock().unwrap().is_empty()
    }

    /// Total number of groups
    fn num_groups(&self) -> usize {
        self.vl_seg_part.read().unwrap().len()
    }
}

// =============================================================================
// Parallel flush coordination for Phase 3 (atomic work-stealing)
// =============================================================================

/// State for coordinating parallel Phase 3 segment storage
/// Workers atomically claim buffers via next_idx, then process independently
struct ParallelFlushState {
    /// Extracted buffers to flush. Each slot has its own Mutex for independent access.
    /// RwLock allows parallel read access to the Vec during Phase 3, avoiding serialization.
    /// Workers only need read access to the Vec to reach their claimed slot's inner Mutex.
    buffers: RwLock<Vec<Mutex<Option<(SegmentGroupKey, SegmentGroupBuffer)>>>>,
    /// Compression results from each buffer (stored by workers, written by Thread 0)
    results: RwLock<Vec<Mutex<Option<FlushPackResult>>>>,
    /// Atomic index for work-stealing (starts at len-1, decrements to -1)
    next_idx: AtomicI32,
}

impl ParallelFlushState {
    fn new() -> Self {
        Self {
            buffers: RwLock::new(Vec::new()),
            results: RwLock::new(Vec::new()),
            next_idx: AtomicI32::new(-1),
        }
    }

    /// Set up buffers to flush and reset atomic counter (called by Thread 0 in Phase 2)
    fn prepare(&self, extracted: Vec<(SegmentGroupKey, SegmentGroupBuffer)>) {
        let len = extracted.len();
        let mut buffers = self.buffers.write().unwrap();
        *buffers = extracted
            .into_iter()
            .map(|(k, b)| Mutex::new(Some((k, b))))
            .collect();
        // Initialize results slots (one per buffer)
        let mut results = self.results.write().unwrap();
        *results = (0..len).map(|_| Mutex::new(None)).collect();
        self.next_idx.store(len as i32 - 1, Ordering::SeqCst);
    }

    /// Claim next buffer index (returns None when all claimed)
    /// Called by ALL workers in Phase 3
    fn claim_next_idx(&self) -> Option<usize> {
        let idx = self.next_idx.fetch_sub(1, Ordering::Relaxed);
        if idx < 0 {
            None
        } else {
            Some(idx as usize)
        }
    }

    /// Get buffer at claimed index (READ lock on Vec, exclusive on slot)
    /// Called by workers after claiming an index in Phase 3
    fn get_buffer_at(&self, idx: usize) -> Option<(SegmentGroupKey, SegmentGroupBuffer)> {
        let buffers = self.buffers.read().unwrap();
        if idx < buffers.len() {
            buffers[idx].lock().unwrap().take()
        } else {
            None
        }
    }

    /// Put buffer back after processing (READ lock on Vec, exclusive on slot)
    fn return_buffer(&self, idx: usize, key: SegmentGroupKey, buffer: SegmentGroupBuffer) {
        let buffers = self.buffers.read().unwrap();
        if idx < buffers.len() {
            *buffers[idx].lock().unwrap() = Some((key, buffer));
        }
    }

    /// Drain all buffers back (called by Thread 0 in Phase 4 - needs WRITE lock)
    fn drain_buffers(&self) -> Vec<(SegmentGroupKey, SegmentGroupBuffer)> {
        let mut buffers = self.buffers.write().unwrap();
        let result: Vec<_> = buffers
            .iter_mut()
            .filter_map(|slot| slot.lock().unwrap().take())
            .collect();
        buffers.clear();
        self.next_idx.store(-1, Ordering::SeqCst);
        result
    }

    /// Store compression result at given index (READ lock on Vec, exclusive on slot)
    fn store_result(&self, idx: usize, result: FlushPackResult) {
        let results = self.results.read().unwrap();
        if idx < results.len() {
            *results[idx].lock().unwrap() = Some(result);
        }
    }

    /// Drain all results sorted by group_id (called by Thread 0 for deterministic writes)
    fn drain_results_sorted(&self) -> Vec<FlushPackResult> {
        let mut results_lock = self.results.write().unwrap();
        let mut all_results: Vec<FlushPackResult> = results_lock
            .iter_mut()
            .filter_map(|slot| slot.lock().unwrap().take())
            .collect();
        results_lock.clear();
        // Sort by group_id for deterministic write order
        all_results.sort_by_key(|r| r.group_id);
        all_results
    }
}

/// Parallel write buffer with per-stream mutexes (C++ AGC pattern: per-segment mutex)
/// Workers operating on different streams don't contend at all.
/// BTreeMap ensures flush writes in sorted stream_id order for determinism.
struct ParallelWriteBuffer {
    /// Per-stream buffers: BTreeMap for sorted iteration, each stream has its own Mutex
    /// RwLock allows concurrent reader access to find the right stream's Mutex
    streams: RwLock<BTreeMap<usize, Mutex<Vec<(Vec<u8>, u64)>>>>,
}

impl ParallelWriteBuffer {
    fn new() -> Self {
        Self {
            streams: RwLock::new(BTreeMap::new()),
        }
    }

    /// Buffer a write for a specific stream (only locks that stream's mutex)
    /// Workers on different streams can call this concurrently without contention
    fn buffer_write(&self, stream_id: usize, data: Vec<u8>, metadata: u64) {
        // First try with read lock - most common case (stream already exists)
        {
            let streams = self.streams.read().unwrap();
            if let Some(stream_mutex) = streams.get(&stream_id) {
                stream_mutex.lock().unwrap().push((data, metadata));
                return;
            }
        }
        // Stream doesn't exist - need write lock to create it
        {
            let mut streams = self.streams.write().unwrap();
            // Double-check (another thread may have created it)
            streams
                .entry(stream_id)
                .or_insert_with(|| Mutex::new(Vec::new()))
                .lock()
                .unwrap()
                .push((data, metadata));
        }
    }

    /// Flush all buffered writes to archive in sorted stream_id order
    /// Called by Thread 0 after barrier - ensures deterministic output
    fn flush_to_archive(&self, archive: &mut Archive) -> Result<()> {
        let streams = self.streams.read().unwrap();
        // BTreeMap iterates in sorted key order (stream_id)
        for (stream_id, stream_mutex) in streams.iter() {
            let parts = stream_mutex.lock().unwrap();
            for (data, metadata) in parts.iter() {
                archive.add_part(*stream_id, data, *metadata)?;
            }
        }
        Ok(())
    }

    /// Clear all buffers (called after flush)
    fn clear(&self) {
        let mut streams = self.streams.write().unwrap();
        for (_, stream_mutex) in streams.iter_mut() {
            stream_mutex.lock().unwrap().clear();
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

    /// Check if this group should write a pack (has >= pack_size segments)
    /// Matches C++ AGC's logic for writing packs when full
    fn should_flush_pack(&self, pack_size: usize) -> bool {
        // Count buffered segments (excluding reference which is handled separately)
        self.segments.len() >= pack_size
    }

    /// Get current segment count (for pack-full detection)
    fn segment_count(&self) -> usize {
        self.segments.len()
    }
}

/// Batch-local state for processing new segments
/// Equivalent to C++ AGC's `m_kmers` local variable in process_new()
/// This is RESET at each sample boundary to match C++ AGC behavior
struct BatchState {
    /// New segments discovered in THIS batch (not found in global registry)
    /// Key: (front_kmer, back_kmer)
    /// Value: Vec of segments with that k-mer pair
    new_segments: BTreeMap<(u64, u64), Vec<PendingSegment>>,

    /// Starting group ID for this batch (continues from global count)
    next_group_id: u32,
}

impl BatchState {
    fn new(starting_group_id: u32) -> Self {
        BatchState {
            new_segments: BTreeMap::new(),
            next_group_id: starting_group_id,
        }
    }

    /// Clear batch state for next sample (resets new_segments map)
    /// next_group_id continues incrementing
    fn clear(&mut self) {
        self.new_segments.clear();
        // next_group_id NOT reset - it continues from where it left off
    }

    /// Add a new segment to this batch
    fn add_segment(&mut self, key: (u64, u64), segment: PendingSegment) {
        self.new_segments.entry(key).or_insert_with(Vec::new).push(segment);
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
    barrier: Arc<std::sync::Barrier>, // Synchronization barrier for batch boundaries (matches C++ AGC bar.arrive_and_wait())
    collection: Arc<Mutex<CollectionV3>>,
    splitters: Arc<AHashSet<u64>>,
    config: StreamingQueueConfig,
    archive: Arc<Mutex<Archive>>,
    segment_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    group_counter: Arc<AtomicU32>, // Starts at 16 for LZ groups
    raw_group_counter: Arc<AtomicU32>, // Round-robin counter for raw groups (0-15)
    reference_sample_name: Arc<Mutex<Option<String>>>, // First sample becomes reference
    // Segment splitting support (Phase 1)
    map_segments: Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>, // (front, back) -> group_id (BTreeMap for deterministic iteration)
    map_segments_terminators: Arc<RwLock<BTreeMap<u64, Vec<u64>>>>, // kmer -> [connected kmers] (BTreeMap for determinism)

    // FFI Grouping Engine - C++ AGC-compatible group assignment
    #[cfg(feature = "cpp_agc")]
    grouping_engine: Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,

    // Persistent reference segment storage (matches C++ AGC v_segments)
    // Stores reference segment data even after groups are flushed, enabling LZ cost estimation
    // for subsequent samples (fixes multi-sample group fragmentation bug)
    reference_segments: Arc<RwLock<BTreeMap<u32, Vec<u8>>>>, // group_id -> reference segment data (BTreeMap for determinism)

    // Reference orientation tracking - stores is_rev_comp for each group's reference segment
    // When a delta segment joins an existing group, it MUST use the same orientation as the reference
    // to ensure LZ encoding works correctly (fixes ZERO_MATCH bug in Case 3 terminator segments)
    reference_orientations: Arc<RwLock<BTreeMap<u32, bool>>>, // group_id -> reference is_rev_comp (BTreeMap for determinism)

    // Track segment splits for renumbering subsequent segments
    // Maps (sample_name, contig_name, original_place) -> number of splits inserted before this position
    split_offsets: Arc<Mutex<BTreeMap<(String, String, usize), usize>>>, // BTreeMap for determinism

    // Priority assignment for interleaved processing (matches C++ AGC)
    // Higher priority = processed first (sample1 > sample2 > sample3...)
    sample_priorities: Arc<RwLock<BTreeMap<String, i32>>>, // sample_name -> priority (BTreeMap for determinism)

    // Track last sample to detect sample boundaries for sync token insertion
    last_sample_name: Arc<Mutex<Option<String>>>, // Last sample that was pushed

    // Batch-local group assignment (matches C++ AGC m_kmers per-batch behavior)
    // When batch_samples reaches batch_size, we flush pending segments and clear batch-local state
    batch_samples: Arc<Mutex<HashSet<String>>>, // Samples in current batch (matches C++ AGC pack_cardinality batch)
    batch_local_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>>, // Batch-local m_kmers equivalent (BTreeMap for deterministic iteration)
    batch_local_terminators: Arc<Mutex<BTreeMap<u64, Vec<u64>>>>, // Batch-local terminators (BTreeMap for determinism)
    pending_batch_segments: Arc<Mutex<Vec<PendingSegment>>>, // Buffer segments until batch boundary
    // Two-tier segment buffering for C++ AGC 4-phase parallel pattern
    buffered_seg_part: Arc<BufferedSegPart>, // Per-group buffers for parallel Phase 1
    // Fallback minimizers map for segments with no terminator match (matches C++ AGC map_fallback_minimizers)
    map_fallback_minimizers: Arc<Mutex<BTreeMap<u64, Vec<(u64, u64)>>>>, // kmer -> [(front, back)] candidate group keys (BTreeMap for determinism)
    next_priority: Arc<Mutex<i32>>, // Decreases for each new sample (starts at i32::MAX)
    next_sequence: Arc<std::sync::atomic::AtomicU64>, // Increases for each contig (FASTA order)
    global_contig_count: Arc<AtomicUsize>, // GLOBAL contig counter for synchronization (C++ AGC: cnt_contigs_in_sample)

    // Deferred metadata streams - written AFTER segment data (C++ AGC compatibility)
    // C++ AGC writes segment data first, then metadata streams at the end
    deferred_file_type_info: (usize, Vec<u8>),    // (stream_id, data)
    deferred_params: (usize, Vec<u8>),            // (stream_id, data)
    deferred_splitters: (usize, Vec<u8>),         // (stream_id, data)
    deferred_segment_splitters: (usize, Vec<u8>), // (stream_id, data)

    // Dynamic splitter discovery for adaptive mode (matches C++ AGC find_new_splitters)
    // Stores reference k-mers to exclude when finding new splitters for non-reference contigs
    ref_singletons: Arc<Vec<u64>>, // Sorted for binary search - reference singleton k-mers (v_candidate_kmers)
    ref_duplicates: Arc<AHashSet<u64>>, // Reference duplicate k-mers (v_duplicated_kmers)

    // Parallel Phase 3 state for atomic work-stealing (matches C++ AGC architecture)
    parallel_state: Arc<ParallelFlushState>,

    // Per-stream write buffer for parallel Phase 3 (C++ AGC pattern: per-segment mutex)
    // Workers on different streams can buffer writes concurrently without contention
    write_buffer: Arc<ParallelWriteBuffer>,

    // RAW segment buffers for deferred classification (parallel Phase 1 optimization)
    // PER-WORKER buffers eliminate contention: each worker pushes to its own buffer
    // Thread 0 drains all buffers at barrier for classification
    raw_segment_buffers: Arc<Vec<Mutex<Vec<RawBufferedSegment>>>>,
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
        splitters: AHashSet<u64>,
    ) -> Result<Self> {
        // Call internal with empty ref data (no dynamic splitter discovery)
        Self::with_splitters_internal(
            output_path,
            config,
            splitters,
            Arc::new(Vec::new()),
            Arc::new(AHashSet::new()),
        )
    }

    /// Internal constructor that accepts all splitter data
    fn with_splitters_internal(
        output_path: impl AsRef<Path>,
        config: StreamingQueueConfig,
        splitters: AHashSet<u64>,
        ref_singletons: Arc<Vec<u64>>,
        ref_duplicates: Arc<AHashSet<u64>>,
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
        // ref_singletons and ref_duplicates are passed as parameters to ensure workers
        // get the same Arc as stored in self (critical for dynamic splitter discovery)

        // Segment grouping for LZ packing (using BTreeMap for better memory efficiency)
        let segment_groups = Arc::new(Mutex::new(BTreeMap::new()));
        let group_counter = Arc::new(AtomicU32::new(NO_RAW_GROUPS)); // Start at 16 (LZ groups), group 0 reserved for orphan segments
        let raw_group_counter = Arc::new(AtomicU32::new(0)); // Round-robin counter for raw groups (0-15)
        let reference_sample_name = Arc::new(Mutex::new(None)); // Shared across all workers

        // Segment splitting support (Phase 1)
        // Initialize map_segments with (MISSING_KMER, MISSING_KMER) → 0
        // This matches C++ AGC line 2396: map_segments[make_pair(~0ull, ~0ull)] = 0
        // All raw segments (both k-mers missing) will map to group 0
        let mut initial_map_segments = BTreeMap::new();
        initial_map_segments.insert(
            SegmentGroupKey {
                kmer_front: MISSING_KMER,
                kmer_back: MISSING_KMER,
            },
            0,
        );
        let map_segments: Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>> = Arc::new(RwLock::new(initial_map_segments));
        let map_segments_terminators: Arc<RwLock<BTreeMap<u64, Vec<u64>>>> = Arc::new(RwLock::new(BTreeMap::new()));
        let split_offsets: Arc<Mutex<BTreeMap<(String, String, usize), usize>>> = Arc::new(Mutex::new(BTreeMap::new()));

        // Persistent reference segment storage (matches C++ AGC v_segments)
        let reference_segments: Arc<RwLock<BTreeMap<u32, Vec<u8>>>> = Arc::new(RwLock::new(BTreeMap::new()));

        // Reference orientation tracking (fixes ZERO_MATCH bug in Case 3 terminator segments)
        let reference_orientations: Arc<RwLock<BTreeMap<u32, bool>>> = Arc::new(RwLock::new(BTreeMap::new()));

        // FFI Grouping Engine - C++ AGC-compatible group assignment
        #[cfg(feature = "cpp_agc")]
        let grouping_engine = Arc::new(Mutex::new(crate::ragc_ffi::GroupingEngine::new(
            config.k as u32,
            NO_RAW_GROUPS,  // Start group IDs at 16 (group 0 reserved for orphan segments)
        )));

        // Priority tracking for interleaved processing (matches C++ AGC)
        let sample_priorities: Arc<RwLock<BTreeMap<String, i32>>> = Arc::new(RwLock::new(BTreeMap::new()));
        let last_sample_name: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None)); // Track last sample for boundary detection
        let next_priority = Arc::new(Mutex::new(i32::MAX)); // Start high, decrease for each sample
        let next_sequence = Arc::new(std::sync::atomic::AtomicU64::new(0)); // Increases for each contig (FASTA order)
        let global_contig_count = Arc::new(AtomicUsize::new(0)); // GLOBAL counter across all samples (C++ AGC: cnt_contigs_in_sample)

        // Batch-local group assignment (matches C++ AGC m_kmers per-batch behavior)
        let batch_samples: Arc<Mutex<HashSet<String>>> = Arc::new(Mutex::new(HashSet::new()));
        let batch_local_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>> = Arc::new(Mutex::new(BTreeMap::new()));
        let batch_local_terminators: Arc<Mutex<BTreeMap<u64, Vec<u64>>>> = Arc::new(Mutex::new(BTreeMap::new()));
        let pending_batch_segments: Arc<Mutex<Vec<PendingSegment>>> = Arc::new(Mutex::new(Vec::new()));
        // Two-tier segment buffering for C++ AGC 4-phase parallel pattern
        let buffered_seg_part: Arc<BufferedSegPart> = Arc::new(BufferedSegPart::new(NO_RAW_GROUPS as usize));
        let map_fallback_minimizers: Arc<Mutex<BTreeMap<u64, Vec<(u64, u64)>>>> = Arc::new(Mutex::new(BTreeMap::new()));

        // Initialize barrier for sample boundary synchronization (matches C++ AGC barrier)
        // All workers must synchronize at sample boundaries to ensure batch flush completes before processing new samples
        let barrier = Arc::new(std::sync::Barrier::new(config.num_threads));

        // Parallel Phase 3 state for atomic work-stealing (matches C++ AGC architecture)
        let parallel_state = Arc::new(ParallelFlushState::new());

        // Per-stream write buffer for parallel Phase 3 (C++ AGC pattern: per-segment mutex)
        // Workers on different streams can buffer writes concurrently without contention
        let write_buffer = Arc::new(ParallelWriteBuffer::new());

        // RAW segment buffers for deferred classification (parallel Phase 1 optimization)
        // PER-WORKER buffers eliminate contention: each worker pushes to its own buffer
        let raw_segment_buffers: Arc<Vec<Mutex<Vec<RawBufferedSegment>>>> = Arc::new(
            (0..config.num_threads)
                .map(|_| Mutex::new(Vec::new()))
                .collect()
        );

        // Spawn worker threads
        let mut workers = Vec::new();
        for worker_id in 0..config.num_threads {
            let queue = Arc::clone(&queue);
            let collection = Arc::clone(&collection);
            let splitters = Arc::clone(&splitters);
            let ref_singletons = Arc::clone(&ref_singletons);
            let ref_duplicates = Arc::clone(&ref_duplicates);
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
            let batch_samples = Arc::clone(&batch_samples);
            let batch_local_groups = Arc::clone(&batch_local_groups);
            let batch_local_terminators = Arc::clone(&batch_local_terminators);
            let pending_batch_segments = Arc::clone(&pending_batch_segments);
            let buffered_seg_part = Arc::clone(&buffered_seg_part);
            let map_fallback_minimizers = Arc::clone(&map_fallback_minimizers);
            let barrier = Arc::clone(&barrier);
            let parallel_state = Arc::clone(&parallel_state);
            let write_buffer = Arc::clone(&write_buffer);
            let raw_segment_buffers = Arc::clone(&raw_segment_buffers);
            let config = config.clone();

            let handle = thread::spawn(move || {
                worker_thread(
                    worker_id,
                    queue,
                    collection,
                    splitters,
                    ref_singletons,
                    ref_duplicates,
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
                    batch_samples,
                    batch_local_groups,
                    batch_local_terminators,
                    pending_batch_segments,
                    buffered_seg_part,
                    map_fallback_minimizers,
                    raw_segment_buffers,
                    barrier,
                    parallel_state,
                    write_buffer,
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
            barrier,
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
            last_sample_name,
            next_priority,
            batch_samples,
            batch_local_groups,
            batch_local_terminators,
            pending_batch_segments,
            buffered_seg_part,
            map_fallback_minimizers,
            next_sequence,
            global_contig_count,
            // Deferred metadata streams (written at end for C++ AGC compatibility)
            deferred_file_type_info,
            deferred_params,
            deferred_splitters,
            deferred_segment_splitters,
            // Dynamic splitter discovery - MUST use the SAME Arcs passed to workers!
            // (empty by default - populated with_full_splitter_data)
            ref_singletons,
            ref_duplicates,
            // Parallel Phase 3 state
            parallel_state,
            // Per-stream write buffer
            write_buffer,
            // Raw segment buffers for deferred classification (per-worker)
            raw_segment_buffers,
        })
    }

    /// Create a new streaming compressor with full splitter data for dynamic discovery
    ///
    /// This is the preferred constructor when using adaptive mode. It accepts:
    /// - `splitters`: Pre-computed splitter k-mers from reference (for initial segmentation)
    /// - `singletons`: All singleton k-mers from reference (for exclusion in find_new_splitters)
    /// - `duplicates`: All duplicate k-mers from reference (for exclusion in find_new_splitters)
    ///
    /// # Arguments
    /// * `output_path` - Path to output AGC archive
    /// * `config` - Compression configuration
    /// * `splitters` - Pre-computed splitter k-mers
    /// * `singletons` - Reference singleton k-mers (sorted Vec for binary search)
    /// * `duplicates` - Reference duplicate k-mers
    pub fn with_full_splitter_data(
        output_path: impl AsRef<Path>,
        config: StreamingQueueConfig,
        splitters: AHashSet<u64>,
        singletons: Vec<u64>,
        duplicates: AHashSet<u64>,
    ) -> Result<Self> {
        // Sort singletons for binary search before creating compressor
        let mut sorted_singletons = singletons;
        sorted_singletons.sort_unstable();

        let verbosity = config.verbosity;
        let ref_singletons = Arc::new(sorted_singletons);
        let ref_duplicates = Arc::new(duplicates);

        if verbosity > 0 {
            eprintln!(
                "  Dynamic splitter discovery enabled: {} ref singletons, {} ref duplicates",
                ref_singletons.len(),
                ref_duplicates.len()
            );
        }

        // Call internal constructor with ref data so workers get the correct Arcs
        Self::with_splitters_internal(output_path, config, splitters, ref_singletons, ref_duplicates)
    }

    /// Create compressor and determine splitters from first contig
    ///
    /// **Note**: This requires at least one contig to be pushed before workers start.
    /// Consider using `with_splitters()` instead if you have a reference genome.
    pub fn new(output_path: impl AsRef<Path>, config: StreamingQueueConfig) -> Result<Self> {
        // Start with empty splitters - will be determined from first push
        Self::with_splitters(output_path, config, AHashSet::new())
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
                let ref_singletons = Arc::clone(&self.ref_singletons);
                let ref_duplicates = Arc::clone(&self.ref_duplicates);
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
                let batch_samples = Arc::clone(&self.batch_samples);
                let batch_local_groups = Arc::clone(&self.batch_local_groups);
                let batch_local_terminators = Arc::clone(&self.batch_local_terminators);
                let pending_batch_segments = Arc::clone(&self.pending_batch_segments);
                let buffered_seg_part = Arc::clone(&self.buffered_seg_part);
                let map_fallback_minimizers = Arc::clone(&self.map_fallback_minimizers);
                let raw_segment_buffers = Arc::clone(&self.raw_segment_buffers);
                let barrier = Arc::clone(&self.barrier);
                let parallel_state = Arc::clone(&self.parallel_state);
                let write_buffer = Arc::clone(&self.write_buffer);
                let config = self.config.clone();

                let handle = thread::spawn(move || {
                    worker_thread(
                        worker_id,
                        queue,
                        collection,
                        splitters,
                        ref_singletons,
                        ref_duplicates,
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
                        batch_samples,
                        batch_local_groups,
                        batch_local_terminators,
                        pending_batch_segments,
                        buffered_seg_part,
                        map_fallback_minimizers,
                        raw_segment_buffers,
                        barrier,
                        parallel_state,
                        write_buffer,
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
        // C++ AGC also decrements priority every 50 contigs WITHIN a sample (max_no_contigs_before_synchronization)
        let sample_priority = {
            let mut priorities = self.sample_priorities.write().unwrap();
            let current_priority = *priorities.entry(sample_name.clone()).or_insert_with(|| {
                // First time seeing this sample - assign new priority
                let mut next_p = self.next_priority.lock().unwrap();
                let priority = *next_p;
                *next_p -= 1; // Decrement for next sample (C++ AGC uses --sample_priority)
                priority
            });

            // Track GLOBAL contig count and insert sync tokens every 50 contigs (pack_cardinality)
            // C++ AGC: if (++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
            // NOTE: Despite the name, C++ AGC's cnt_contigs_in_sample is GLOBAL, not per-sample!
            // FIX 5: Only send PACK_BOUNDARY sync tokens in concatenated mode (single file)
            // In non-concatenated mode (multiple files), only SAMPLE_BOUNDARY sync tokens are sent
            let count = self.global_contig_count.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            let need_sync = self.config.concatenated_genomes && (count + 1) % self.config.pack_size == 0;

            if need_sync {
                // Reached synchronization point (every 50 contigs GLOBALLY)
                // C++ AGC does: cnt_contigs_in_sample = 0; --sample_priority;
                if let Some(priority) = priorities.get_mut(&sample_name) {
                    *priority -= 1;
                }

                // Get the NEW priority (after decrement) for sync tokens
                let new_priority = *priorities.get(&sample_name).unwrap();

                // Drop locks before inserting sync tokens to avoid deadlock
                drop(priorities);

                // Insert sync tokens (matches C++ AGC EmplaceManyNoCost)
                // CRITICAL: Sync tokens must have HIGHER priority than subsequent contigs
                // to ensure they're processed before any contigs with the new_priority.
                if self.config.verbosity > 0 {
                    eprintln!(
                        "PACK_BOUNDARY: Inserting {} sync tokens after {} contigs (global count)",
                        self.config.num_threads, count + 1
                    );
                }

                for _ in 0..self.config.num_threads {
                    let sync_token = ContigTask {
                        sample_name: sample_name.clone(),
                        contig_name: String::from("<SYNC>"),
                        data: Vec::new(),
                        // Use large priority boost to ensure sync tokens are processed BEFORE any contigs
                        // With +1, contigs with same priority but higher cost were being popped first
                        // This caused barrier deadlock when some workers exited before others got sync tokens
                        sample_priority: new_priority + 1_000_000,
                        cost: 0,
                        sequence,
                        is_sync_token: true,
                    };
                    self.queue.push(sync_token, 0)?;
                }

                // Return NEW priority for subsequent contigs
                new_priority
            } else {
                current_priority  // Use priority BEFORE potential decrement (this contig uses current priority)
            }
        };

        // Insert sync tokens at sample boundaries (matches C++ AGC registration tokens)
        // OPTIMIZATION: In multi-file mode, SKIP per-sample sync tokens for better parallelism
        // This batches all samples together - sync only happens at finalization
        // Set RAGC_SYNC_PER_SAMPLE=1 to force per-sample sync (matches old behavior)
        {
            let mut last_sample = self.last_sample_name.lock().unwrap();
            if let Some(ref last) = *last_sample {
                if last != &sample_name {
                    // Sample boundary detected
                    // Only insert sync tokens if forced by env var (for debugging/compatibility)
                    let force_sync = std::env::var("RAGC_SYNC_PER_SAMPLE").map(|v| v == "1").unwrap_or(false);

                    if force_sync {
                        if self.config.verbosity > 0 {
                            eprintln!(
                                "SAMPLE_BOUNDARY: Inserting {} sync tokens (transitioning from {} to {})",
                                self.config.num_threads, last, sample_name
                            );
                        }

                        // Insert num_threads sync tokens (matches C++ AGC EmplaceManyNoCost)
                        // All workers must pop a token and synchronize before processing new sample
                        // CRITICAL: Sync tokens must have MUCH HIGHER priority than any contigs
                        // to ensure they're pulled and processed BEFORE any contigs.
                        // Use large priority boost (+1_000_000) to overcome cost-based tie-breaking
                        // which was causing contigs to be popped before sync tokens at same priority.
                        for _ in 0..self.config.num_threads {
                            let sync_token = ContigTask {
                                sample_name: sample_name.clone(),
                                contig_name: String::from("<SYNC>"),
                                data: Vec::new(), // Empty data for sync token
                                sample_priority: sample_priority + 1_000_000, // Much higher priority than any contigs
                                cost: 0, // No cost for sync tokens
                                sequence,
                                is_sync_token: true,
                            };
                            self.queue.push(sync_token, 0)?; // 0 size for sync tokens
                        }
                    } else if self.config.verbosity > 1 {
                        eprintln!(
                            "SAMPLE_BOUNDARY: SKIPPING sync tokens (multi-file batching: {} -> {})",
                            last, sample_name
                        );
                    }
                }
            }
            // Update last sample name
            *last_sample = Some(sample_name.clone());
        }

        // Create task with priority information
        // NOTE: sequence is used for FASTA ordering (lower = processed first)
        let cost = data.len(); // C++ AGC: auto cost = contig.size()
        let task = ContigTask {
            sample_name: sample_name.clone(),
            contig_name,
            data,
            sample_priority,
            cost,
            sequence,
            is_sync_token: false, // Normal contig task, not a sync token
        };

        // Push to queue (BLOCKS if queue is full!)
        // Queue is now a priority queue - highest priority processed first
        // eprintln!("[RAGC PUSH] sample={} contig={} priority={} cost={} sequence={}",
        //           &task.sample_name, &task.contig_name, task.sample_priority, task.cost, task.sequence);
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

    /// Insert sync tokens to trigger incremental compression of buffered segments.
    /// Call this after pushing a batch of samples to process them incrementally
    /// instead of waiting for finalize().
    pub fn sync_and_flush(&self, sample_name: &str) -> Result<()> {
        // Insert sync tokens for each worker
        let sequence = self.next_sequence.fetch_add(1, std::sync::atomic::Ordering::SeqCst);

        for _ in 0..self.config.num_threads {
            let sync_token = ContigTask {
                sample_name: format!("<SYNC:{}>", sample_name),
                contig_name: String::from("<SYNC>"),
                data: Vec::new(),
                sample_priority: 1_000_000_i32, // High priority = processed after pending contigs
                cost: 0,
                sequence,
                is_sync_token: true,
            };
            self.queue.push(sync_token, 0)?;
        }

        // Wait for sync tokens to be processed (queue empty)
        while self.queue.len() > 0 {
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        Ok(())
    }

    pub fn finalize(self) -> Result<()> {
        if self.config.verbosity > 0 {
            eprintln!("Finalizing compression...");
        }

        // CRITICAL: Insert FINAL sync tokens before closing queue
        // This ensures buffered_seg_part data is processed and flushed
        // (matches C++ AGC line 2236-2244: final sync at end of input)
        if self.config.verbosity > 0 {
            eprintln!("  Inserting {} final sync tokens...", self.config.num_threads);
        }

        // Use sequence 0 and high priority to ensure sync tokens are processed last
        let sequence = 0;

        for _ in 0..self.config.num_threads {
            let sync_token = ContigTask {
                sample_name: String::from("<FINAL>"),
                contig_name: String::from("<SYNC>"),
                data: Vec::new(),
                sample_priority: 1_000_000_i32, // Very high priority = processed after all real contigs
                cost: 0,
                sequence,
                is_sync_token: true,
            };
            self.queue.push(sync_token, 0)?;
        }

        if self.config.verbosity > 0 {
            eprintln!("  Closing queue...");
        }

        // Close queue - no more pushes allowed
        self.queue.close();

        if self.config.verbosity > 0 {
            eprintln!("  Waiting for {} workers to finish...", self.workers.len());
        }

        let wait_start = std::time::Instant::now();
        // Wait for all workers to finish
        for (i, handle) in self.workers.into_iter().enumerate() {
            handle
                .join()
                .expect("Worker thread panicked")
                .with_context(|| format!("Worker {} failed", i))?;
        }

        if self.config.verbosity > 0 {
            eprintln!("FINALIZE_TIMING: Wait for workers took {:?}", wait_start.elapsed());
            eprintln!("All workers finished!");
            eprintln!("Flushing remaining segment packs...");
        }

        // Flush all remaining partial packs using PARALLEL compression
        let flush_start = std::time::Instant::now();
        {
            use crate::segment_compression::compress_segment_configured;
            use rayon::prelude::*;

            let mut groups = self.segment_groups.lock().unwrap();
            let num_groups = groups.len();

            // Phase 1: Flush any groups with pending segments (rare, usually 0-1)
            for (key, buffer) in groups.iter_mut() {
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
            }

            // Phase 2: Collect and PARALLEL compress pending_deltas
            // Each entry: (stream_id, raw_data, compressed_data, raw_size)
            struct PartialPackData {
                stream_id: usize,
                raw_data: Vec<u8>,
                compressed: Vec<u8>,
                raw_size: usize,
                use_compressed: bool,
            }

            let compression_level = self.config.compression_level;
            let verbosity = self.config.verbosity;

            // Extract work items from groups
            let work_items: Vec<_> = groups.iter_mut()
                .filter(|(_, buffer)| !buffer.pending_deltas.is_empty())
                .map(|(_, buffer)| {
                    let use_lz_encoding = buffer.group_id > 0;
                    let mut packed_data = Vec::new();

                    if !use_lz_encoding && !buffer.raw_placeholder_written {
                        packed_data.push(0x7f);
                        packed_data.push(CONTIG_SEPARATOR);
                    }

                    for delta in buffer.pending_deltas.iter() {
                        packed_data.extend_from_slice(delta);
                        packed_data.push(CONTIG_SEPARATOR);
                    }

                    let stream_id = buffer.stream_id as usize;
                    let group_id = buffer.group_id;
                    let delta_count = buffer.pending_deltas.len();

                    buffer.pending_deltas.clear();
                    buffer.pending_delta_ids.clear();

                    (stream_id, packed_data, group_id, delta_count)
                })
                .collect();

            // Parallel compression using rayon
            let compressed_packs: Vec<PartialPackData> = work_items
                .into_par_iter()
                .filter_map(|(stream_id, packed_data, group_id, delta_count)| {
                    if packed_data.is_empty() {
                        return None;
                    }

                    let raw_size = packed_data.len();
                    let mut compressed = match compress_segment_configured(&packed_data, compression_level) {
                        Ok(c) => c,
                        Err(e) => {
                            eprintln!("Error compressing final partial pack for group {}: {}", group_id, e);
                            return None;
                        }
                    };
                    compressed.push(0); // Marker 0 = plain ZSTD

                    let use_compressed = compressed.len() < raw_size;

                    if verbosity > 1 {
                        eprintln!(
                            "  Compressed final partial pack for group {} with {} deltas",
                            group_id, delta_count
                        );
                    }

                    Some(PartialPackData {
                        stream_id,
                        raw_data: packed_data,
                        compressed,
                        raw_size,
                        use_compressed,
                    })
                })
                .collect();

            // Phase 3: Sequential writes to archive (sorted by stream_id for determinism)
            let mut sorted_packs = compressed_packs;
            sorted_packs.sort_by_key(|p| p.stream_id);

            let mut arch = self.archive.lock().unwrap();
            for pack in sorted_packs {
                if pack.use_compressed {
                    arch.add_part(pack.stream_id, &pack.compressed, pack.raw_size as u64)
                        .context("Failed to write compressed final partial pack")?;
                } else {
                    arch.add_part(pack.stream_id, &pack.raw_data, 0)
                        .context("Failed to write uncompressed final partial pack")?;
                }
            }
            drop(arch);

            if self.config.verbosity > 0 {
                eprintln!("Flushed {} segment groups", num_groups);
                eprintln!("FINALIZE_TIMING: Flush took {:?}", flush_start.elapsed());
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
/// Pre-compressed data ready for archive write (no locks needed during compression)
struct PreCompressedPart {
    stream_id: usize,
    data: Vec<u8>,
    metadata: u64,
}

/// Segment registration data for collection (batched for single lock acquisition)
struct SegmentRegistration {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    group_id: u32,
    in_group_id: u32,
    is_rev_comp: bool,
    raw_length: u32,
}

/// Result of parallel compression phase (for deterministic sequential writes)
/// Workers produce these in parallel, then Thread 0 writes them in sorted order
struct FlushPackResult {
    group_id: u32,
    archive_writes: Vec<PreCompressedPart>,
    registrations: Vec<SegmentRegistration>,
    ref_to_store: Option<(u32, Vec<u8>)>,
}

fn flush_pack(
    buffer: &mut SegmentGroupBuffer,
    collection: &Arc<Mutex<CollectionV3>>,
    archive: &Arc<Mutex<Archive>>,
    config: &StreamingQueueConfig,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
) -> Result<()> {
    use crate::segment_compression::{compress_reference_segment, compress_segment_configured};

    // Skip if no segments to write (but still write reference if present)
    if buffer.segments.is_empty() && buffer.ref_written {
        return Ok(());
    }

    let use_lz_encoding = buffer.group_id >= NO_RAW_GROUPS;

    // CRITICAL FIX: Sort ALL segments FIRST by (sample_name, contig_name, seg_part_no)
    // BEFORE picking the reference. This ensures the lexicographically first segment
    // becomes the reference, matching C++ AGC's behavior.
    // (Previous code sorted AFTER picking reference, causing wrong reference selection)
    buffer.segments.sort();

    // ============================================================
    // PHASE 1: Compress everything WITHOUT holding any locks
    // ============================================================

    // Collect all pre-compressed writes for batched archive write
    let mut archive_writes: Vec<PreCompressedPart> = Vec::new();
    // Collect all segment registrations for batched collection update
    let mut registrations: Vec<SegmentRegistration> = Vec::new();
    // Reference data to store in global map (if any)
    let mut ref_to_store: Option<(u32, Vec<u8>)> = None;

    // Write reference segment if not already written (first pack for this group)
    // Extract reference from sorted segments (matching C++ AGC: first segment after sort becomes reference)
    // NOTE: Raw groups (0-15) do NOT have a reference - all segments stored raw
    if use_lz_encoding && !buffer.ref_written && !buffer.segments.is_empty() {
        // Remove first segment (alphabetically first after sorting) to use as reference
        let ref_seg = buffer.segments.remove(0);

        if crate::env_cache::debug_ref_write() {
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

        // Compress reference using adaptive compression (NO LOCK)
        let (mut compressed, marker) = compress_reference_segment(&ref_seg.data)
            .context("Failed to compress reference")?;
        compressed.push(marker);

        // Metadata stores the uncompressed size
        let ref_size = ref_seg.data.len() as u64;

        // CRITICAL: Check if compression helped (matching C++ AGC segment.h lines 179, 204)
        // C++ AGC: if(packed_size + 1u < (uint32_t) data.size())
        // If compression didn't help, write UNCOMPRESSED raw data with metadata=0
        if compressed.len() < ref_seg.data.len() {
            // Compression helped - write compressed data with metadata=original_size
            archive_writes.push(PreCompressedPart {
                stream_id: buffer.ref_stream_id,
                data: compressed,
                metadata: ref_size,
            });
        } else {
            // Compression didn't help - write UNCOMPRESSED data with metadata=0
            archive_writes.push(PreCompressedPart {
                stream_id: buffer.ref_stream_id,
                data: ref_seg.data.clone(),
                metadata: 0,
            });
        }

        // Queue reference registration
        registrations.push(SegmentRegistration {
            sample_name: ref_seg.sample_name.clone(),
            contig_name: ref_seg.contig_name.clone(),
            seg_part_no: ref_seg.seg_part_no,
            group_id: buffer.group_id,
            in_group_id: 0, // Reference is always at position 0
            is_rev_comp: ref_seg.is_rev_comp,
            raw_length: ref_seg.data.len() as u32,
        });

        buffer.ref_written = true;

        // Queue reference for global map storage
        ref_to_store = Some((buffer.group_id, ref_seg.data.clone()));

        buffer.reference_segment = Some(ref_seg.clone()); // Store for LZ encoding

        // Prepare LZ encoder with reference (matching C++ AGC segment.cpp line 43: lz_diff->Prepare(s))
        // This is done ONCE when the reference is written, then reused for all subsequent segments
        if use_lz_encoding {
            let mut lz = LZDiff::new(config.min_match_len as u32);
            lz.prepare(&ref_seg.data);
            buffer.lz_diff = Some(lz);
        }
    }

    // NOTE: Segments are already sorted at the start of flush_pack (line ~1003)
    // This sort was moved earlier to ensure correct reference selection.

    // Pack segments together with delta deduplication (matching C++ AGC segment.cpp lines 66-74)
    // Note: segments do NOT include the reference - it's stored separately
    //
    // CRITICAL FIX: Partial packs must persist across flush_pack calls to ensure pack boundaries
    // align with decompression expectations. Pack N must contain entries for in_group_ids
    // (N*50)+1 to (N+1)*50. Only write a pack when it has exactly 50 entries (or at finalization).
    // Use buffer.pending_deltas and buffer.pending_delta_ids to persist partial packs.

    let mut segment_in_group_ids: Vec<(usize, u32)> = Vec::new(); // (segment_index, in_group_id) for each segment

    // Helper function to compress a complete pack (exactly 50 entries) - NO LOCK
    let compress_pack = |deltas: &[Vec<u8>], needs_raw_placeholder: bool, stream_id: usize, compression_level: i32| -> Result<PreCompressedPart> {
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

        let total_raw_size = packed_data.len();
        let mut compressed = compress_segment_configured(&packed_data, compression_level)
            .context("Failed to compress pack")?;
        compressed.push(0); // Marker 0 = plain ZSTD

        if compressed.len() < total_raw_size {
            Ok(PreCompressedPart {
                stream_id,
                data: compressed,
                metadata: total_raw_size as u64,
            })
        } else {
            Ok(PreCompressedPart {
                stream_id,
                data: packed_data,
                metadata: 0,
            })
        }
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

            // Compare with C++ AGC encode (TEST HARNESS)
            #[cfg(feature = "cpp_agc")]
            if crate::env_cache::test_lz_encoding() {
                if let Some(ref_seg) = &buffer.reference_segment {
                    if let Some(cpp_encoded) = crate::ragc_ffi::lzdiff_v2_encode(
                        &ref_seg.data,
                        &seg.data,
                        config.min_match_len as u32,
                    ) {
                        if ragc_encoded != cpp_encoded {
                            eprintln!("\n========================================");
                            eprintln!("🔥 LZ ENCODING MISMATCH DETECTED!");
                            eprintln!("========================================");
                            eprintln!("Group:          {}", buffer.group_id);
                            eprintln!("Sample:         {}", seg.sample_name);
                            eprintln!("Contig:         {}", seg.contig_name);
                            eprintln!("Segment:        {}", seg.seg_part_no);
                            eprintln!("Reference len:  {}", ref_seg.data.len());
                            eprintln!("Target len:     {}", seg.data.len());
                            eprintln!("RAGC encoded:   {} bytes", ragc_encoded.len());
                            eprintln!("C++ AGC encoded: {} bytes", cpp_encoded.len());
                            eprintln!("Difference:     {} bytes", (ragc_encoded.len() as i64 - cpp_encoded.len() as i64).abs());
                            eprintln!();

                            // Find first difference
                            let mut first_diff_byte = None;
                            for (i, (r, c)) in ragc_encoded.iter().zip(cpp_encoded.iter()).enumerate() {
                                if r != c {
                                    first_diff_byte = Some(i);
                                    break;
                                }
                            }

                            if let Some(i) = first_diff_byte {
                                eprintln!("First difference at byte {}", i);
                                let start = if i > 20 { i - 20 } else { 0 };
                                let end = (i + 30).min(ragc_encoded.len()).min(cpp_encoded.len());

                                eprintln!("\nRAGC output around difference:");
                                let ragc_hex: Vec<_> = ragc_encoded[start..end].iter().map(|b| format!("{:02x}", b)).collect();
                                let ragc_ascii: String = ragc_encoded[start..end].iter().map(|&b| {
                                    if b >= 32 && b < 127 { b as char } else { '.' }
                                }).collect();
                                eprintln!("  Hex:   {}", ragc_hex.join(" "));
                                eprintln!("  ASCII: {}", ragc_ascii);

                                eprintln!("\nC++ AGC output around difference:");
                                let cpp_hex: Vec<_> = cpp_encoded[start..end].iter().map(|b| format!("{:02x}", b)).collect();
                                let cpp_ascii: String = cpp_encoded[start..end].iter().map(|&b| {
                                    if b >= 32 && b < 127 { b as char } else { '.' }
                                }).collect();
                                eprintln!("  Hex:   {}", cpp_hex.join(" "));
                                eprintln!("  ASCII: {}", cpp_ascii);

                                eprintln!("\nByte at position {}:", i);
                                eprintln!("  RAGC:    0x{:02x} ('{}')", ragc_encoded[i],
                                    if ragc_encoded[i] >= 32 && ragc_encoded[i] < 127 { ragc_encoded[i] as char } else { '?' });
                                eprintln!("  C++ AGC: 0x{:02x} ('{}')", cpp_encoded[i],
                                    if cpp_encoded[i] >= 32 && cpp_encoded[i] < 127 { cpp_encoded[i] as char } else { '?' });
                            } else if ragc_encoded.len() != cpp_encoded.len() {
                                eprintln!("Encodings match for first {} bytes, but lengths differ",
                                    ragc_encoded.len().min(cpp_encoded.len()));
                                if ragc_encoded.len() > cpp_encoded.len() {
                                    let extra_start = cpp_encoded.len();
                                    let extra_hex: Vec<_> = ragc_encoded[extra_start..].iter().take(40).map(|b| format!("{:02x}", b)).collect();
                                    let extra_ascii: String = ragc_encoded[extra_start..].iter().take(40).map(|&b| {
                                        if b >= 32 && b < 127 { b as char } else { '.' }
                                    }).collect();
                                    eprintln!("RAGC has {} extra bytes:", ragc_encoded.len() - cpp_encoded.len());
                                    eprintln!("  Hex:   {}", extra_hex.join(" "));
                                    eprintln!("  ASCII: {}", extra_ascii);
                                } else {
                                    let extra_start = ragc_encoded.len();
                                    let extra_hex: Vec<_> = cpp_encoded[extra_start..].iter().take(40).map(|b| format!("{:02x}", b)).collect();
                                    let extra_ascii: String = cpp_encoded[extra_start..].iter().take(40).map(|&b| {
                                        if b >= 32 && b < 127 { b as char } else { '.' }
                                    }).collect();
                                    eprintln!("C++ AGC has {} extra bytes:", cpp_encoded.len() - ragc_encoded.len());
                                    eprintln!("  Hex:   {}", extra_hex.join(" "));
                                    eprintln!("  ASCII: {}", extra_ascii);
                                }
                            }

                            // Show last 10 bytes of each
                            eprintln!("\nLast 10 bytes of each encoding:");
                            let ragc_tail_start = if ragc_encoded.len() > 10 { ragc_encoded.len() - 10 } else { 0 };
                            let ragc_tail_hex: Vec<_> = ragc_encoded[ragc_tail_start..].iter().map(|b| format!("{:02x}", b)).collect();
                            let ragc_tail_ascii: String = ragc_encoded[ragc_tail_start..].iter().map(|&b| {
                                if b >= 32 && b < 127 { b as char } else { '.' }
                            }).collect();
                            eprintln!("RAGC    (bytes {}-{}):", ragc_tail_start, ragc_encoded.len()-1);
                            eprintln!("  Hex:   {}", ragc_tail_hex.join(" "));
                            eprintln!("  ASCII: {}", ragc_tail_ascii);

                            let cpp_tail_start = if cpp_encoded.len() > 10 { cpp_encoded.len() - 10 } else { 0 };
                            let cpp_tail_hex: Vec<_> = cpp_encoded[cpp_tail_start..].iter().map(|b| format!("{:02x}", b)).collect();
                            let cpp_tail_ascii: String = cpp_encoded[cpp_tail_start..].iter().map(|&b| {
                                if b >= 32 && b < 127 { b as char } else { '.' }
                            }).collect();
                            eprintln!("C++ AGC (bytes {}-{}):", cpp_tail_start, cpp_encoded.len()-1);
                            eprintln!("  Hex:   {}", cpp_tail_hex.join(" "));
                            eprintln!("  ASCII: {}", cpp_tail_ascii);

                            eprintln!("\n========================================");
                            eprintln!("Aborting on first LZ encoding mismatch!");
                            eprintln!("========================================\n");

                            panic!("LZ encoding mismatch detected - see details above");
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
                // Compress pack WITHOUT holding any lock
                let needs_placeholder = !use_lz_encoding && !buffer.raw_placeholder_written;
                let pack = compress_pack(&buffer.pending_deltas, needs_placeholder, buffer.stream_id, config.compression_level)?;
                archive_writes.push(pack);
                buffer.raw_placeholder_written = true;

                // Clear for next pack - deduplication starts fresh
                buffer.pending_deltas.clear();
                buffer.pending_delta_ids.clear();
            }
        }
    }

    // DO NOT write partial pack here - leave it in buffer.pending_deltas for next flush_pack call
    // Partial packs are only written in finalize() to ensure pack boundaries align with decompression

    // Queue segment registrations (batched for single lock acquisition)
    for &(seg_idx, in_group_id) in segment_in_group_ids.iter() {
        let seg = &buffer.segments[seg_idx];
        registrations.push(SegmentRegistration {
            sample_name: seg.sample_name.clone(),
            contig_name: seg.contig_name.clone(),
            seg_part_no: seg.seg_part_no,
            group_id: buffer.group_id,
            in_group_id,
            is_rev_comp: seg.is_rev_comp,
            raw_length: seg.data.len() as u32,
        });
    }

    // ============================================================
    // PHASE 2: Batched writes with minimal lock duration
    // ============================================================

    // Write all pre-compressed data to archive (SINGLE lock acquisition)
    if !archive_writes.is_empty() {
        let mut arch = archive.lock().unwrap();
        for part in archive_writes {
            arch.add_part(part.stream_id, &part.data, part.metadata)
                .context("Failed to write to archive")?;
        }
    }

    // Store reference in global map (if any)
    if let Some((group_id, ref_data)) = ref_to_store {
        let mut ref_segs = reference_segments.write().unwrap();
        ref_segs.insert(group_id, ref_data);
    }

    // Register all segments in collection (SINGLE lock acquisition)
    if !registrations.is_empty() {
        let mut coll = collection.lock().unwrap();
        for reg in registrations {
            coll.add_segment_placed(
                &reg.sample_name,
                &reg.contig_name,
                reg.seg_part_no,
                reg.group_id,
                reg.in_group_id,
                reg.is_rev_comp,
                reg.raw_length,
            )
            .context("Failed to register segment")?;
        }
    }

    // Clear segments for next batch (but keep pending_deltas!)
    buffer.segments.clear();

    Ok(())
}

/// Compress-only version of flush_pack for deterministic parallel compression.
/// Workers call this in parallel to produce FlushPackResult, then Thread 0
/// writes all results in sorted group_id order for deterministic archives.
fn flush_pack_compress_only(
    buffer: &mut SegmentGroupBuffer,
    config: &StreamingQueueConfig,
) -> Result<FlushPackResult> {
    use crate::segment_compression::{compress_reference_segment, compress_segment_configured};

    let mut archive_writes: Vec<PreCompressedPart> = Vec::new();
    let mut registrations: Vec<SegmentRegistration> = Vec::new();
    let mut ref_to_store: Option<(u32, Vec<u8>)> = None;

    // Skip if no segments to write (but still write reference if present)
    if buffer.segments.is_empty() && buffer.ref_written {
        return Ok(FlushPackResult {
            group_id: buffer.group_id,
            archive_writes,
            registrations,
            ref_to_store,
        });
    }

    let use_lz_encoding = buffer.group_id >= NO_RAW_GROUPS;

    // Sort segments for deterministic reference selection
    buffer.segments.sort();

    // Write reference segment if not already written
    if use_lz_encoding && !buffer.ref_written && !buffer.segments.is_empty() {
        let ref_seg = buffer.segments.remove(0);

        // Compress reference
        let (mut compressed, marker) = compress_reference_segment(&ref_seg.data)
            .context("Failed to compress reference")?;
        compressed.push(marker);

        let ref_size = ref_seg.data.len() as u64;

        if compressed.len() < ref_seg.data.len() {
            archive_writes.push(PreCompressedPart {
                stream_id: buffer.ref_stream_id,
                data: compressed,
                metadata: ref_size,
            });
        } else {
            archive_writes.push(PreCompressedPart {
                stream_id: buffer.ref_stream_id,
                data: ref_seg.data.clone(),
                metadata: 0,
            });
        }

        registrations.push(SegmentRegistration {
            sample_name: ref_seg.sample_name.clone(),
            contig_name: ref_seg.contig_name.clone(),
            seg_part_no: ref_seg.seg_part_no,
            group_id: buffer.group_id,
            in_group_id: 0,
            is_rev_comp: ref_seg.is_rev_comp,
            raw_length: ref_seg.data.len() as u32,
        });

        buffer.ref_written = true;
        ref_to_store = Some((buffer.group_id, ref_seg.data.clone()));
        buffer.reference_segment = Some(ref_seg.clone());

        if use_lz_encoding {
            let mut lz = LZDiff::new(config.min_match_len as u32);
            lz.prepare(&ref_seg.data);
            buffer.lz_diff = Some(lz);
        }
    }

    // Compress pack helper (same as flush_pack)
    let compress_pack = |deltas: &[Vec<u8>], needs_raw_placeholder: bool, stream_id: usize, compression_level: i32| -> Result<PreCompressedPart> {
        let mut packed_data = Vec::new();

        if needs_raw_placeholder {
            packed_data.push(0x7f);
            packed_data.push(CONTIG_SEPARATOR);
        }

        for delta in deltas.iter() {
            packed_data.extend_from_slice(delta);
            packed_data.push(CONTIG_SEPARATOR);
        }

        let total_raw_size = packed_data.len();
        let mut compressed = compress_segment_configured(&packed_data, compression_level)
            .context("Failed to compress pack")?;
        compressed.push(0);

        if compressed.len() < total_raw_size {
            Ok(PreCompressedPart {
                stream_id,
                data: compressed,
                metadata: total_raw_size as u64,
            })
        } else {
            Ok(PreCompressedPart {
                stream_id,
                data: packed_data,
                metadata: 0,
            })
        }
    };

    let mut segment_in_group_ids: Vec<(usize, u32)> = Vec::new();

    for (seg_idx, seg) in buffer.segments.iter().enumerate() {
        let contig_data = if !use_lz_encoding || buffer.reference_segment.is_none() {
            seg.data.clone()
        } else {
            buffer.lz_diff.as_mut()
                .expect("lz_diff should be prepared")
                .encode(&seg.data)
        };

        if use_lz_encoding && contig_data.is_empty() {
            segment_in_group_ids.push((seg_idx, 0));
            continue;
        }

        if let Some(existing_idx) = buffer.pending_deltas.iter().position(|d| d == &contig_data) {
            let reused_id = buffer.pending_delta_ids[existing_idx];
            segment_in_group_ids.push((seg_idx, reused_id));
        } else {
            buffer.segments_written = buffer.segments_written.max(1);
            let in_group_id = buffer.segments_written;
            buffer.segments_written += 1;
            buffer.pending_delta_ids.push(in_group_id);
            segment_in_group_ids.push((seg_idx, in_group_id));
            buffer.pending_deltas.push(contig_data);

            if buffer.pending_deltas.len() == PACK_CARDINALITY {
                let needs_placeholder = !use_lz_encoding && !buffer.raw_placeholder_written;
                let pack = compress_pack(&buffer.pending_deltas, needs_placeholder, buffer.stream_id, config.compression_level)?;
                archive_writes.push(pack);
                buffer.raw_placeholder_written = true;
                buffer.pending_deltas.clear();
                buffer.pending_delta_ids.clear();
            }
        }
    }

    for &(seg_idx, in_group_id) in segment_in_group_ids.iter() {
        let seg = &buffer.segments[seg_idx];
        registrations.push(SegmentRegistration {
            sample_name: seg.sample_name.clone(),
            contig_name: seg.contig_name.clone(),
            seg_part_no: seg.seg_part_no,
            group_id: buffer.group_id,
            in_group_id,
            is_rev_comp: seg.is_rev_comp,
            raw_length: seg.data.len() as u32,
        });
    }

    buffer.segments.clear();

    Ok(FlushPackResult {
        group_id: buffer.group_id,
        archive_writes,
        registrations,
        ref_to_store,
    })
}

/// Write reference segment immediately when first segment arrives in group
/// (Matches C++ AGC segment.cpp lines 41-48: if (no_seqs == 0) writes reference right away)
/// This ensures LZ encoding works correctly for subsequent segments
fn write_reference_immediately(
    segment: &BufferedSegment,
    buffer: &mut SegmentGroupBuffer,
    collection: &Arc<Mutex<CollectionV3>>,
    archive: &Arc<Mutex<Archive>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    reference_orientations: &Arc<RwLock<BTreeMap<u32, bool>>>,
    config: &StreamingQueueConfig,
) -> Result<()> {
    use crate::segment_compression::compress_reference_segment;

    if crate::env_cache::debug_ref_write() {
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
        let mut ref_segs = reference_segments.write().unwrap();
        ref_segs.insert(buffer.group_id, segment.data.clone());
    }

    // 4c. Store reference orientation for ZERO_MATCH bug fix
    // When a delta segment joins this group later, it MUST use the same orientation
    // as the reference to ensure LZ encoding works correctly
    {
        let mut ref_orients = reference_orientations.write().unwrap();
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
    map_segments_terminators: &Arc<RwLock<BTreeMap<u64, Vec<u64>>>>,
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    config: &StreamingQueueConfig,
) -> (u64, u64, bool) {
    let segment_len = segment_data.len();
    use crate::segment::MISSING_KMER;

    // Look up kmer in terminators map to find connected k-mers
    let connected_kmers = {
        let terminators = map_segments_terminators.read().unwrap();
        match terminators.get(&kmer) {
            Some(vec) => vec.clone(),
            None => {
                // No connections found - create new group with MISSING
                // Match C++ AGC lines 1671-1679: check is_dir_oriented()
                // Debug: log entry to no-connection path
                if crate::env_cache::debug_is_dir() {
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
    if crate::env_cache::debug_is_dir() {
        eprintln!("RAGC_FIND_GROUP_FOUND_CONN: kmer={} kmer_is_dir={} connections={:?}",
            kmer, kmer_is_dir, connected_kmers);
    }

    // Build list of candidate groups
    // Each candidate: (key_front, key_back, needs_rc, ref_segment_size)
    let mut candidates: Vec<(u64, u64, bool, usize)> = Vec::new();

    // OPTIMIZATION: Reduce lock scope - first collect candidate keys, then look up ref sizes
    // This minimizes the time segment_groups.lock() is held

    // Phase 1: Build candidate orderings (no locks needed)
    let mut candidate_keys: Vec<(u64, u64, bool, SegmentGroupKey)> = Vec::new();
    for &cand_kmer in &connected_kmers {
        // Create candidate group key normalized (smaller, larger)
        // C++ AGC lines 1691-1704
        //
        // IMPORTANT: When cand_kmer is MISSING, we need to try BOTH orderings!
        // Groups with MISSING k-mers can be stored as either (MISSING, kmer) or (kmer, MISSING)
        // depending on kmer_is_dir when they were created. We must match the actual stored key.
        let orderings: Vec<(u64, u64, bool)> = if cand_kmer == MISSING_KMER {
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

        for (key_front, key_back, needs_rc) in orderings {
            let cand_key = SegmentGroupKey {
                kmer_front: key_front,
                kmer_back: key_back,
            };
            candidate_keys.push((key_front, key_back, needs_rc, cand_key));
        }
    }

    // Phase 2: Quick check which candidates exist (brief locks)
    // First pass: check global registry (RwLock - can be concurrent)
    let mut existing_candidates: Vec<(u64, u64, bool, Option<u32>)> = Vec::new();
    {
        let seg_map = map_segments.read().unwrap();
        for (key_front, key_back, needs_rc, cand_key) in &candidate_keys {
            if let Some(&group_id) = seg_map.get(cand_key) {
                existing_candidates.push((*key_front, *key_back, *needs_rc, Some(group_id)));
            }
        }
    } // seg_map lock released

    // Second pass: check batch-local buffer for remaining candidates (Mutex - exclusive)
    // Only if we didn't find candidates in global registry
    if existing_candidates.is_empty() {
        let groups = segment_groups.lock().unwrap();
        let mut already_found = std::collections::HashSet::new();
        for (key_front, key_back, needs_rc, cand_key) in &candidate_keys {
            if groups.contains_key(cand_key) {
                // Get ref size from buffer
                let ref_size = if let Some(group_buffer) = groups.get(cand_key) {
                    if let Some(ref_seg) = &group_buffer.reference_segment {
                        ref_seg.data.len()
                    } else {
                        segment_len
                    }
                } else {
                    segment_len
                };

                // Debug trace
                if crate::env_cache::debug_is_dir() {
                    eprintln!("RAGC_FIND_GROUP_CAND_CHECK: cand_key=({},{}) exists_in_groups=true ref_size={}",
                        key_front, key_back, ref_size);
                }

                // Use cand_kmer as key to deduplicate (only one match per connected_kmer)
                let connected = if *key_front == kmer { *key_back } else { *key_front };
                if !already_found.contains(&connected) {
                    candidates.push((*key_front, *key_back, *needs_rc, ref_size));
                    already_found.insert(connected);
                }
            }
        }
    } // groups lock released

    // Phase 3: Get ref sizes for global candidates (brief RwLock)
    if !existing_candidates.is_empty() {
        let ref_segs = reference_segments.read().unwrap();
        let mut already_found = std::collections::HashSet::new();
        for (key_front, key_back, needs_rc, group_id_opt) in existing_candidates {
            let ref_size = if let Some(group_id) = group_id_opt {
                if let Some(ref_data) = ref_segs.get(&group_id) {
                    ref_data.len()
                } else {
                    segment_len
                }
            } else {
                segment_len
            };

            // Debug trace
            if crate::env_cache::debug_is_dir() {
                eprintln!("RAGC_FIND_GROUP_CAND_CHECK: cand_key=({},{}) exists_in_seg_map=true ref_size={}",
                    key_front, key_back, ref_size);
            }

            // Use cand_kmer as key to deduplicate (only one match per connected_kmer)
            let connected = if key_front == kmer { key_back } else { key_front };
            if !already_found.contains(&connected) {
                candidates.push((key_front, key_back, needs_rc, ref_size));
                already_found.insert(connected);
            }
        }
    } // ref_segs lock released

    if candidates.is_empty() {
        // No existing groups found - create new with MISSING
        // Must match C++ AGC is_dir_oriented logic (same as no-connections case above)
        if crate::env_cache::debug_is_dir() {
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
        let seg_map = map_segments.read().unwrap();
        let ref_segs = reference_segments.read().unwrap();

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
    map_fallback_minimizers: &Arc<Mutex<BTreeMap<u64, Vec<(u64, u64)>>>>,
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
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
    let mut cand_seg_counts: BTreeMap<(u64, u64), Vec<u64>> = BTreeMap::new(); // BTreeMap for determinism

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
        let seg_map = map_segments.read().unwrap();
        let ref_segs = reference_segments.read().unwrap();

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
    map_fallback_minimizers: &Arc<Mutex<BTreeMap<u64, Vec<(u64, u64)>>>>,
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

// =============================================================================
// Parallel batch processing functions (C++ AGC 4-phase pattern)
// =============================================================================

/// Phase 2: Prepare batch for parallel processing
/// - Processes NEW segments from buffered_seg_part (assigns group_ids)
/// - Drains segments from buffered_seg_part into SegmentGroupBuffer entries
/// - Extracts buffers that need flushing into ParallelFlushState
/// Returns true if there are buffers to flush, false otherwise
#[allow(clippy::too_many_arguments)]
fn prepare_batch_parallel(
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    buffered_seg_part: &Arc<BufferedSegPart>,
    batch_local_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>>,
    batch_local_terminators: &Arc<Mutex<BTreeMap<u64, Vec<u64>>>>,
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    group_counter: &Arc<AtomicU32>,
    raw_group_counter: &Arc<AtomicU32>,
    archive: &Arc<Mutex<Archive>>,
    collection: &Arc<Mutex<CollectionV3>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    reference_orientations: &Arc<RwLock<BTreeMap<u32, bool>>>,
    #[cfg(feature = "cpp_agc")]
    grouping_engine: &Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,
    parallel_state: &ParallelFlushState,
    config: &StreamingQueueConfig,
) -> Result<bool> {
    use crate::segment::MISSING_KMER;

    // Check if there's anything to process
    let batch_map_len = batch_local_groups.lock().unwrap().len();
    let batch_terms_len = batch_local_terminators.lock().unwrap().len();
    let has_buffered_segments = buffered_seg_part.has_segments();

    if batch_map_len == 0 && batch_terms_len == 0 && !has_buffered_segments {
        return Ok(false);
    }

    if config.verbosity > 0 {
        eprintln!("PREPARE_BATCH_PARALLEL: Processing {} batch-local groups, buffered segments: {}, {} terminator keys",
            batch_map_len, has_buffered_segments, batch_terms_len);
    }

    // Phase 2a: Process NEW segments - assign group_ids deterministically
    // Get current group counter and update after process_new
    let mut next_group_id = group_counter.load(Ordering::SeqCst);
    {
        let mut global_map = map_segments.write().unwrap();
        let new_count = buffered_seg_part.process_new(&mut global_map, &mut next_group_id);
        if config.verbosity > 0 && new_count > 0 {
            eprintln!("PREPARE_BATCH_PARALLEL: Assigned {} new group IDs", new_count);
        }
    }
    // Update the shared counter
    group_counter.store(next_group_id, Ordering::SeqCst);

    // Phase 2b: Sort segments within each group for determinism
    buffered_seg_part.sort_known();

    // Phase 2c: Drain segments from buffered_seg_part into SegmentGroupBuffer entries
    let mut groups_map = segment_groups.lock().unwrap();

    // Build reverse lookup map (group_id -> key) ONCE to avoid O(n²) lookups
    let group_id_to_key: std::collections::HashMap<u32, SegmentGroupKey> = {
        let global_map = map_segments.read().unwrap();
        global_map.iter().map(|(k, &gid)| (gid, k.clone())).collect()
    };

    // Phase 2c-1: Collect all segments with their keys (no locks needed)
    let num_groups = buffered_seg_part.num_groups();
    let mut collected_segments: Vec<(u32, SegmentGroupKey, BufferedSegment)> = Vec::new();
    for group_id in 0..num_groups as u32 {
        while let Some(seg) = buffered_seg_part.get_part(group_id) {
            let key = group_id_to_key.get(&group_id).cloned().unwrap_or_else(|| {
                SegmentGroupKey {
                    kmer_front: MISSING_KMER,
                    kmer_back: MISSING_KMER,
                }
            });
            collected_segments.push((group_id, key, seg));
        }
    }

    // Phase 2c-2: Batch update batch_local_groups (ONE lock acquisition)
    {
        let mut batch_map = batch_local_groups.lock().unwrap();
        for (group_id, key, _) in &collected_segments {
            batch_map.insert(key.clone(), *group_id);
        }
    }

    // Phase 2c-3: Register with FFI engine (ONE lock acquisition)
    #[cfg(feature = "cpp_agc")]
    {
        let mut eng = grouping_engine.lock().unwrap();
        for (group_id, key, _) in &collected_segments {
            if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                eng.register_group(key.kmer_front, key.kmer_back, *group_id);
            }
        }
    }

    // Phase 2c-4: Batch update terminators (ONE lock acquisition)
    {
        let mut term_map = batch_local_terminators.lock().unwrap();
        for (_, key, _) in &collected_segments {
            if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
                term_map.entry(key.kmer_front)
                    .or_insert_with(Vec::new)
                    .push(key.kmer_back);
                if key.kmer_front != key.kmer_back {
                    term_map.entry(key.kmer_back)
                        .or_insert_with(Vec::new)
                        .push(key.kmer_front);
                }
            }
        }
    }

    // Phase 2c-5: Pre-register all streams for new groups (ONE lock acquisition)
    // Build a set of existing group_ids first for O(1) lookup
    let existing_group_ids: std::collections::HashSet<u32> = groups_map
        .values()
        .map(|b| b.group_id)
        .collect();

    // Collect unique new group_ids (O(n) instead of O(n×m))
    let new_group_ids: std::collections::HashSet<u32> = collected_segments.iter()
        .map(|(gid, _, _)| *gid)
        .filter(|gid| !existing_group_ids.contains(gid))
        .collect();

    // Pre-register all streams in one lock acquisition
    let stream_registrations: std::collections::HashMap<u32, (usize, usize)> = if !new_group_ids.is_empty() {
        let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
        let mut arch = archive.lock().unwrap();
        new_group_ids.iter().map(|&group_id| {
            let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
            let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);
            let stream_id = arch.register_stream(&delta_stream_name);
            let ref_stream_id = arch.register_stream(&ref_stream_name);
            (group_id, (stream_id, ref_stream_id))
        }).collect()
    } else {
        std::collections::HashMap::new()
    };

    // Phase 2c-6: Add segments to buffers (groups_map already locked)
    for (group_id, key, seg) in collected_segments {
        let buffer = groups_map.entry(key.clone()).or_insert_with(|| {
            let (stream_id, ref_stream_id) = stream_registrations.get(&group_id)
                .copied()
                .unwrap_or_else(|| {
                    // Fallback: register now (shouldn't happen if logic is correct)
                    let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
                    let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
                    let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);
                    let mut arch = archive.lock().unwrap();
                    let sid = arch.register_stream(&delta_stream_name);
                    let rsid = arch.register_stream(&ref_stream_name);
                    (sid, rsid)
                });
            SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
        });
        buffer.segments.push(seg);
    }

    // Clear the buffered_seg_part after draining
    buffered_seg_part.clear();

    // Extract buffers that need flushing
    let mut extracted: Vec<(SegmentGroupKey, SegmentGroupBuffer)> = Vec::new();
    let keys_to_remove: Vec<SegmentGroupKey> = groups_map
        .iter()
        .filter(|(_, buffer)| !buffer.segments.is_empty() || !buffer.ref_written)
        .map(|(k, _)| k.clone())
        .collect();

    // Sort keys for deterministic processing order
    let mut sorted_keys = keys_to_remove;
    sorted_keys.sort();

    for key in sorted_keys {
        if let Some(buffer) = groups_map.remove(&key) {
            extracted.push((key, buffer));
        }
    }

    let has_work = !extracted.is_empty();

    if config.verbosity > 0 {
        eprintln!("PREPARE_BATCH_PARALLEL: Extracted {} buffers for parallel flush", extracted.len());
    }

    // Populate ParallelFlushState
    parallel_state.prepare(extracted);

    Ok(has_work)
}

/// Phase 4: Cleanup after parallel processing
/// - Re-inserts processed buffers
/// - Updates global maps
/// - Clears batch-local state
fn cleanup_batch_parallel(
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    batch_local_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>>,
    batch_local_terminators: &Arc<Mutex<BTreeMap<u64, Vec<u64>>>>,
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: &Arc<RwLock<BTreeMap<u64, Vec<u64>>>>,
    parallel_state: &ParallelFlushState,
    config: &StreamingQueueConfig,
) {
    // Re-insert processed buffers
    let processed = parallel_state.drain_buffers();
    {
        let mut groups_map = segment_groups.lock().unwrap();
        for (key, buffer) in processed {
            groups_map.insert(key, buffer);
        }
    }

    // Update global registry with batch-local groups
    {
        let batch_map = batch_local_groups.lock().unwrap();
        let mut global_map = map_segments.write().unwrap();
        for (key, group_id) in batch_map.iter() {
            global_map.entry(key.clone()).or_insert(*group_id);
        }
    }

    // Merge batch-local terminators into global terminators
    {
        let batch_terms = batch_local_terminators.lock().unwrap();
        let mut global_terms = map_segments_terminators.write().unwrap();
        for (kmer, connections) in batch_terms.iter() {
            let entry = global_terms.entry(*kmer).or_insert_with(Vec::new);
            entry.extend(connections.iter().cloned());
            entry.sort_unstable();
            entry.dedup();
        }
    }

    // Clear batch-local state
    batch_local_groups.lock().unwrap().clear();
    batch_local_terminators.lock().unwrap().clear();

    if config.verbosity > 0 {
        eprintln!("CLEANUP_BATCH_PARALLEL: Batch cleanup complete");
    }
}

/// Classify raw segments at barrier (Thread 0 only)
/// This eliminates lock contention by doing all classification single-threaded.
/// Raw segments are sorted for determinism, then classified using the same
/// Case 2/3a/3b logic as before, just without contention.
fn classify_raw_segments_at_barrier(
    raw_segment_buffers: &Arc<Vec<Mutex<Vec<RawBufferedSegment>>>>,
    buffered_seg_part: &Arc<BufferedSegPart>,
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: &Arc<RwLock<BTreeMap<u64, Vec<u64>>>>,
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    config: &StreamingQueueConfig,
) {
    use crate::segment::MISSING_KMER;

    // Drain all raw segments from ALL per-worker buffers into one Vec
    let mut raw_segs: Vec<RawBufferedSegment> = Vec::new();
    for buffer in raw_segment_buffers.iter() {
        let mut worker_segs = buffer.lock().unwrap();
        raw_segs.append(&mut *worker_segs);
    }

    if raw_segs.is_empty() {
        return;
    }

    // Sort for determinism: by sample_name, contig_name, original_place
    raw_segs.sort();

    if config.verbosity > 0 {
        eprintln!("CLASSIFY_RAW_BARRIER: Processing {} raw segments (single-threaded)", raw_segs.len());
    }

    // Classify each segment (NO LOCK CONTENTION - single-threaded)
    for raw_seg in raw_segs.drain(..) {
        // Case 2/3a/3b classification (same logic as before)
        let (key_front, key_back, should_reverse) =
            if raw_seg.front_kmer != MISSING_KMER && raw_seg.back_kmer != MISSING_KMER {
                // Case 2: Both k-mers present
                if raw_seg.front_kmer < raw_seg.back_kmer {
                    (raw_seg.front_kmer, raw_seg.back_kmer, false)
                } else {
                    (raw_seg.back_kmer, raw_seg.front_kmer, true)
                }
            } else if raw_seg.front_kmer != MISSING_KMER {
                // Case 3a: Only front k-mer present
                let (kf, kb, sr) = find_group_with_one_kmer(
                    raw_seg.front_kmer,
                    raw_seg.front_kmer_is_dir,
                    &raw_seg.data,
                    &raw_seg.data_rc,
                    map_segments_terminators,
                    map_segments,
                    segment_groups,
                    reference_segments,
                    config,
                );
                (kf, kb, sr)
            } else if raw_seg.back_kmer != MISSING_KMER {
                // Case 3b: Only back k-mer present
                let kmer_is_dir_after_swap = !raw_seg.back_kmer_is_dir;
                let (kf, kb, mut sr) = find_group_with_one_kmer(
                    raw_seg.back_kmer,
                    kmer_is_dir_after_swap,
                    &raw_seg.data_rc,
                    &raw_seg.data,
                    map_segments_terminators,
                    map_segments,
                    segment_groups,
                    reference_segments,
                    config,
                );
                sr = !sr;
                (kf, kb, sr)
            } else {
                // Case 1: Both MISSING - use raw grouping
                (MISSING_KMER, MISSING_KMER, false)
            };

        let key = SegmentGroupKey {
            kmer_front: key_front,
            kmer_back: key_back,
        };

        // Prepare segment data (reverse complement if needed)
        let segment_data = if should_reverse {
            raw_seg.data_rc.clone()
        } else {
            raw_seg.data
        };

        // Check if group exists (NO CONTENTION - we're single-threaded)
        let group_id_opt = {
            let seg_map = map_segments.read().unwrap();
            seg_map.get(&key).copied()
        };

        if let Some(group_id) = group_id_opt {
            // KNOWN: add to per-group buffer
            buffered_seg_part.add_known(group_id, BufferedSegment {
                sample_name: raw_seg.sample_name,
                contig_name: raw_seg.contig_name,
                seg_part_no: raw_seg.original_place,
                data: segment_data,
                is_rev_comp: should_reverse,
                sample_priority: raw_seg.sample_priority,
            });
        } else {
            // NEW: add to s_seg_part
            buffered_seg_part.add_new(NewSegment {
                kmer_front: key.kmer_front,
                kmer_back: key.kmer_back,
                sample_priority: raw_seg.sample_priority,
                sample_name: raw_seg.sample_name,
                contig_name: raw_seg.contig_name,
                seg_part_no: raw_seg.original_place,
                data: segment_data,
                should_reverse,
            });
        }
    }

    if config.verbosity > 0 {
        eprintln!("CLASSIFY_RAW_BARRIER: Classification complete");
    }
}

/// Flush batch-local groups to global state (matches C++ AGC batch boundary)
/// This updates the global map_segments registry with batch-local groups,
/// then clears the batch-local state (like C++ AGC destroying m_kmers at batch end)
fn flush_batch(
    segment_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    pending_batch_segments: &Arc<Mutex<Vec<PendingSegment>>>,
    batch_local_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>>,
    batch_local_terminators: &Arc<Mutex<BTreeMap<u64, Vec<u64>>>>,
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    group_counter: &Arc<AtomicU32>,
    raw_group_counter: &Arc<AtomicU32>, // FIX 17: Round-robin counter for raw groups (0-15)
    map_segments_terminators: &Arc<RwLock<BTreeMap<u64, Vec<u64>>>>,
    archive: &Arc<Mutex<Archive>>,
    collection: &Arc<Mutex<CollectionV3>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    reference_orientations: &Arc<RwLock<BTreeMap<u32, bool>>>,
    #[cfg(feature = "cpp_agc")]
    grouping_engine: &Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,
    config: &StreamingQueueConfig,
) -> Result<()> {
    use crate::segment::MISSING_KMER;

    // Get pending segments and check if anything needs flushing
    let mut pending = pending_batch_segments.lock().unwrap();
    let batch_map_len = batch_local_groups.lock().unwrap().len();
    let batch_terms_len = batch_local_terminators.lock().unwrap().len();

    if batch_map_len == 0 && batch_terms_len == 0 && pending.is_empty() {
        #[cfg(feature = "verbose_debug")]
        if config.verbosity > 0 {
            eprintln!("FLUSH_BATCH: No pending groups to flush");
        }
        return Ok(());
    }

    #[cfg(feature = "verbose_debug")]
    if config.verbosity > 0 {
        eprintln!("FLUSH_BATCH: Processing {} batch-local groups, {} pending segments, {} terminator keys",
            batch_map_len, pending.len(), batch_terms_len);
    }

    // CRITICAL: Sort pending segments by (sample, contig, place) before assigning group_ids
    // This matches C++ AGC's BTreeSet iteration order (agc_compressor.cpp process_new())
    pending.sort();

    if config.verbosity > 1 && !pending.is_empty() {
        eprintln!("FLUSH_BATCH: Sorted {} pending segments for group_id assignment", pending.len());
    }

    // Process sorted pending segments - assign group_ids and write to archive
    let mut groups_map = segment_groups.lock().unwrap();

    for pend in pending.iter() {
        // Assign group_id: Orphan segments (both k-mers MISSING) distributed across raw groups 0-15
        // For segments with k-mers: lookup existing group, or create new group if not found
        // FIX 17: Distribute orphan segments across groups 0-15 (round-robin) to match C++ AGC's
        // distribute_segments(0, 0, no_raw_groups) behavior (agc_compressor.cpp line 986)
        let group_id = if pend.key.kmer_back == MISSING_KMER && pend.key.kmer_front == MISSING_KMER {
            // Round-robin distribution across raw groups 0-15
            raw_group_counter.fetch_add(1, Ordering::SeqCst) % NO_RAW_GROUPS
        } else {
            // Check if this k-mer pair already has a group assigned
            let mut global_map = map_segments.write().unwrap();
            if let Some(&existing_group_id) = global_map.get(&pend.key) {
                // Use existing group
                if crate::env_cache::trace_group() {
                    eprintln!("GROUPING_LOOKUP_HIT: sample={} contig={} place={} front={} back={} found_group={}",
                        pend.sample_name, pend.contig_name, pend.place,
                        pend.key.kmer_front, pend.key.kmer_back, existing_group_id);
                }
                drop(global_map);
                existing_group_id
            } else {
                // Create new group
                let new_group_id = group_counter.fetch_add(1, Ordering::SeqCst);
                if crate::env_cache::trace_group() {
                    eprintln!("GROUPING_LOOKUP_MISS: sample={} contig={} place={} front={} back={} creating_group={} (map has {} entries)",
                        pend.sample_name, pend.contig_name, pend.place,
                        pend.key.kmer_front, pend.key.kmer_back, new_group_id, global_map.len());
                }
                global_map.insert(pend.key.clone(), new_group_id);
                drop(global_map);
                new_group_id
            }
        };

        if config.verbosity > 2 {
            eprintln!("FLUSH_BATCH_ASSIGN: group_id={} front={} back={} sample={} contig={} place={}",
                group_id, pend.key.kmer_front, pend.key.kmer_back,
                pend.sample_name, pend.contig_name, pend.place);
        }

        // Register orphan segments to global map (non-orphans already registered above)
        if pend.key.kmer_back == MISSING_KMER && pend.key.kmer_front == MISSING_KMER {
            let mut global_map = map_segments.write().unwrap();
            global_map.insert(pend.key.clone(), group_id);
        }

        // TRACE: Log when segments from AAA#0 are registered
        if crate::env_cache::trace_group() && pend.sample_name.contains("AAA#0") {
            let global_map = map_segments.read().unwrap();
            eprintln!("TRACE_REGISTER: sample={} contig={} place={} front={} back={} group_id={} (map_segments now has {} entries)",
                pend.sample_name, pend.contig_name, pend.place,
                pend.key.kmer_front, pend.key.kmer_back, group_id, global_map.len());
        }

        // Register to batch-local map
        {
            let mut batch_map = batch_local_groups.lock().unwrap();
            batch_map.insert(pend.key.clone(), group_id);
        }

        // Register with FFI engine
        #[cfg(feature = "cpp_agc")]
        if pend.key.kmer_front != MISSING_KMER && pend.key.kmer_back != MISSING_KMER {
            let mut eng = grouping_engine.lock().unwrap();
            eng.register_group(pend.key.kmer_front, pend.key.kmer_back, group_id);
        }

        // Update batch-local terminators (will be merged to global below)
        // Only for LZ groups (both k-mers non-MISSING)
        if pend.key.kmer_front != MISSING_KMER && pend.key.kmer_back != MISSING_KMER {
            let mut term_map = batch_local_terminators.lock().unwrap();

            term_map.entry(pend.key.kmer_front)
                .or_insert_with(Vec::new)
                .push(pend.key.kmer_back);

            if pend.key.kmer_front != pend.key.kmer_back {
                term_map.entry(pend.key.kmer_back)
                    .or_insert_with(Vec::new)
                    .push(pend.key.kmer_front);
            }
        }

        // Get or create SegmentGroupBuffer for this group
        let buffer = groups_map.entry(pend.key.clone()).or_insert_with(|| {
            // Register streams
            let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
            let delta_stream_name = ragc_common::stream_delta_name(archive_version, group_id);
            let ref_stream_name = ragc_common::stream_ref_name(archive_version, group_id);

            let mut arch = archive.lock().unwrap();
            let stream_id = arch.register_stream(&delta_stream_name);
            let ref_stream_id = arch.register_stream(&ref_stream_name);
            drop(arch);

            SegmentGroupBuffer::new(group_id, stream_id, ref_stream_id)
        });

        // Add to buffer or write as reference
        let is_raw_group = group_id < NO_RAW_GROUPS;  // Groups 0-15 are raw groups (match C++ AGC)
        if !is_raw_group && buffer.reference_segment.is_none() && buffer.segments.is_empty() {
            // First segment in LZ group - write as reference immediately
            // Create BufferedSegment with original orientation (reference sets the group orientation)
            let buffered = BufferedSegment {
                sample_name: pend.sample_name.clone(),
                contig_name: pend.contig_name.clone(),
                seg_part_no: pend.place,
                data: pend.segment_data.clone(),
                is_rev_comp: pend.should_reverse,
                sample_priority: pend.sample_priority,
            };
            if let Err(e) = write_reference_immediately(
                &buffered, buffer, collection, archive, reference_segments, reference_orientations, config
            ) {
                eprintln!("ERROR in flush_batch: Failed to write reference: {}", e);
                buffer.segments.push(buffered);
            }
        } else {
            // Delta segment (joining existing group)
            // FIX 18: Do NOT adjust orientation to match reference - C++ AGC stores each segment
            // with its own computed is_rev_comp based on k-mer comparison (front < back -> false,
            // front >= back -> true). Segments in the same group can have different is_rev_comp.
            let buffered = BufferedSegment {
                sample_name: pend.sample_name.clone(),
                contig_name: pend.contig_name.clone(),
                seg_part_no: pend.place,
                data: pend.segment_data.clone(),
                is_rev_comp: pend.should_reverse,
                sample_priority: pend.sample_priority,
            };
            buffer.segments.push(buffered);
        }

        // FIX 4: Removed mid-batch pack flush to match C++ AGC's batch-level sorting
        // C++ AGC calls sort_known() on ALL segments in batch BEFORE writing ANY
        // Flushing mid-batch would write segments in pack-level sorted order, not batch-level
        // All groups will be flushed at end of batch (after loop) instead
    }

    // Clear pending segments
    pending.clear();
    drop(pending);

    // FIX 4: Flush all group buffers at end of batch (match C++ AGC's sort_known + store_segments)
    // This ensures segments within each group are sorted globally across the entire batch,
    // not just within individual packs. Matches C++ AGC architecture:
    // - C++ AGC: register_segments() calls sort_known() on ALL segments, then store_segments() writes ALL
    // - RAGC: Accumulate all segments for batch, then flush_pack() sorts + writes at end
    for (_key, buffer) in groups_map.iter_mut() {
        if !buffer.segments.is_empty() || !buffer.ref_written {
            flush_pack(buffer, collection, archive, config, reference_segments)
                .context("Failed to flush pack at end of batch")?;
        }
    }

    drop(groups_map);

    // Update global registry with batch-local groups (from existing group processing)
    {
        let batch_map = batch_local_groups.lock().unwrap();
        let mut global_map = map_segments.write().unwrap();
        for (key, group_id) in batch_map.iter() {
            global_map.entry(key.clone()).or_insert(*group_id);
        }
    }

    // CRITICAL: Merge batch-local terminators into global terminators
    // This is where C++ AGC makes terminators visible for find_middle in subsequent samples
    {
        let batch_terms = batch_local_terminators.lock().unwrap();
        let mut global_terms = map_segments_terminators.write().unwrap();
        for (kmer, connections) in batch_terms.iter() {
            let entry = global_terms.entry(*kmer).or_insert_with(Vec::new);
            entry.extend(connections.iter().cloned());
            entry.sort_unstable();
            entry.dedup();
        }
    }

    // Clear batch-local state (like C++ AGC destroying m_kmers)
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
    _key: &SegmentGroupKey,
    _map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    _batch_local_groups: &Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>>,
    _reference_orientations: &Arc<RwLock<BTreeMap<u32, bool>>>,
) -> (Vec<u8>, bool) {
    // FIX 18: Do NOT adjust orientation to match reference - C++ AGC stores each segment
    // with its own computed is_rev_comp based on k-mer comparison. Segments in the same
    // group can have different is_rev_comp values.
    (data.to_vec(), should_reverse)
}

/// Worker thread that pulls from queue and compresses
fn worker_thread(
    worker_id: usize,
    queue: Arc<MemoryBoundedQueue<ContigTask>>,
    collection: Arc<Mutex<CollectionV3>>,
    splitters: Arc<AHashSet<u64>>,
    ref_singletons: Arc<Vec<u64>>,      // For dynamic splitter discovery (sorted)
    ref_duplicates: Arc<AHashSet<u64>>,  // For dynamic splitter discovery
    archive: Arc<Mutex<Archive>>,
    segment_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    group_counter: Arc<AtomicU32>,
    raw_group_counter: Arc<AtomicU32>,
    reference_sample_name: Arc<Mutex<Option<String>>>,
    map_segments: Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: Arc<RwLock<BTreeMap<u64, Vec<u64>>>>,
    reference_segments: Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    reference_orientations: Arc<RwLock<BTreeMap<u32, bool>>>,
    split_offsets: Arc<Mutex<BTreeMap<(String, String, usize), usize>>>,
    #[cfg(feature = "cpp_agc")]
    grouping_engine: Arc<Mutex<crate::ragc_ffi::GroupingEngine>>,
    batch_samples: Arc<Mutex<HashSet<String>>>,
    batch_local_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, u32>>>,
    batch_local_terminators: Arc<Mutex<BTreeMap<u64, Vec<u64>>>>,
    pending_batch_segments: Arc<Mutex<Vec<PendingSegment>>>,
    buffered_seg_part: Arc<BufferedSegPart>, // Per-group buffers for parallel Phase 1
    map_fallback_minimizers: Arc<Mutex<BTreeMap<u64, Vec<(u64, u64)>>>>,
    raw_segment_buffers: Arc<Vec<Mutex<Vec<RawBufferedSegment>>>>, // Per-worker buffers for deferred classification
    barrier: Arc<std::sync::Barrier>, // Synchronization barrier for batch boundaries
    parallel_state: Arc<ParallelFlushState>, // Shared state for parallel Phase 3
    write_buffer: Arc<ParallelWriteBuffer>, // Per-stream buffers for parallel writes
    config: StreamingQueueConfig,
) -> Result<()> {
    let mut processed_count = 0;

    // Create fallback filter from config
    let fallback_filter = FallbackFilter::new(config.fallback_frac);

    // Timing accumulators for performance analysis
    let mut total_queue_wait = std::time::Duration::ZERO;
    let mut total_segment_processing = std::time::Duration::ZERO;
    let mut total_barrier_wait = std::time::Duration::ZERO;
    let mut total_sync_processing = std::time::Duration::ZERO;
    let mut contig_count = 0usize;
    let mut sync_count = 0usize;

    loop {
        // Pull from queue (blocks if empty, returns None when closed)
        let queue_start = std::time::Instant::now();
        let Some(task) = queue.pull() else {
            // Print timing summary on exit
            if config.verbosity > 0 {
                eprintln!("Worker {} TIMING: queue_wait={:?} segment_proc={:?} barrier_wait={:?} sync_proc={:?} contigs={} syncs={}",
                    worker_id, total_queue_wait, total_segment_processing, total_barrier_wait, total_sync_processing, contig_count, sync_count);
            }
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
                &raw_group_counter, // FIX 17: Pass raw_group_counter for round-robin distribution
                &map_segments_terminators,
                &archive,
                &collection,
                &reference_segments,
                &reference_orientations,
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

        let queue_wait = queue_start.elapsed();
        total_queue_wait += queue_wait;

        // Handle sync tokens with barrier synchronization (matches C++ AGC registration stage)
        if task.is_sync_token {
            let sync_start = std::time::Instant::now();
            sync_count += 1;
            if config.verbosity > 0 {
                eprintln!("Worker {} hit sync token for sample {}", worker_id, task.sample_name);
            }

            // =================================================================
            // C++ AGC 4-Phase Parallel Pattern
            // =================================================================

            // Barrier 1: All workers arrive at sample boundary
            let barrier_start = std::time::Instant::now();
            barrier.wait();
            total_barrier_wait += barrier_start.elapsed();

            // Phase 2 (Thread 0 only): Classify raw segments and prepare batch
            if worker_id == 0 {
                if config.verbosity > 0 {
                    eprintln!("Worker 0 preparing batch at sample boundary for {}", task.sample_name);
                }

                let phase2_start = std::time::Instant::now();

                // Step 1: Classify all raw segments (deferred from parallel segment loop)
                // This eliminates lock contention by doing classification single-threaded
                classify_raw_segments_at_barrier(
                    &raw_segment_buffers,
                    &buffered_seg_part,
                    &map_segments,
                    &map_segments_terminators,
                    &segment_groups,
                    &reference_segments,
                    &config,
                );

                let classify_time = phase2_start.elapsed();
                if config.verbosity > 0 {
                    eprintln!("TIMING: Classification took {:?}", classify_time);
                }

                // Step 2: Prepare batch for parallel compression
                let prepare_start = std::time::Instant::now();
                prepare_batch_parallel(
                    &segment_groups,
                    &buffered_seg_part,
                    &batch_local_groups,
                    &batch_local_terminators,
                    &map_segments,
                    &group_counter,
                    &raw_group_counter,
                    &archive,
                    &collection,
                    &reference_segments,
                    &reference_orientations,
                    #[cfg(feature = "cpp_agc")]
                    &grouping_engine,
                    &parallel_state,
                    &config,
                )?;
                let prepare_time = prepare_start.elapsed();
                if config.verbosity > 0 {
                    eprintln!("TIMING: Prepare took {:?}", prepare_time);
                }
            }

            // Barrier 2: All workers see prepared buffers
            let barrier_start = std::time::Instant::now();
            barrier.wait();
            total_barrier_wait += barrier_start.elapsed();

            let compress_start = std::time::Instant::now();
            // Phase 3a (ALL workers): Atomic work-stealing to COMPRESS and BUFFER writes
            // Workers compress segments and buffer archive writes (C++ AGC: AddPartBuffered)
            // Buffering is fast (memory only), flush happens after barrier
            loop {
                let Some(idx) = parallel_state.claim_next_idx() else {
                    break;
                };

                if let Some((key, mut buffer)) = parallel_state.get_buffer_at(idx) {
                    // Compress this buffer
                    if !buffer.segments.is_empty() || !buffer.ref_written {
                        match flush_pack_compress_only(&mut buffer, &config) {
                            Ok(mut result) => {
                                // Buffer archive writes using per-stream mutexes (NO global lock!)
                                // Workers on different streams can buffer concurrently
                                for part in result.archive_writes.drain(..) {
                                    write_buffer.buffer_write(part.stream_id, part.data, part.metadata);
                                }
                                // Store result (now without archive_writes)
                                parallel_state.store_result(idx, result);
                            }
                            Err(e) => {
                                eprintln!("Worker {} error compressing group {}: {}", worker_id, buffer.group_id, e);
                            }
                        }
                    }
                    // Return buffer
                    parallel_state.return_buffer(idx, key, buffer);
                }
            }

            // Barrier 3: All workers done with compression and buffering
            let barrier_start = std::time::Instant::now();
            barrier.wait();
            total_barrier_wait += barrier_start.elapsed();

            if worker_id == 0 && config.verbosity > 0 {
                eprintln!("TIMING: Compression took {:?} (all workers)", compress_start.elapsed());
            }

            // Phase 3b + Phase 4 (Thread 0 only): Flush writes, registrations, and cleanup
            // Combined to reduce barrier overhead (was 2 separate barriers)
            if worker_id == 0 {
                // Phase 3b: Flush buffered writes and process registrations
                let sorted_results = parallel_state.drain_results_sorted();

                // Take all locks once at the start
                let mut arch = archive.lock().unwrap();
                let mut ref_segs = reference_segments.write().unwrap();
                let mut coll = collection.lock().unwrap();

                // Flush all buffered writes from per-stream buffer to archive
                // BTreeMap ensures sorted stream_id order for determinism
                if let Err(e) = write_buffer.flush_to_archive(&mut *arch) {
                    eprintln!("Thread 0 error flushing archive buffers: {}", e);
                }
                // Clear write buffer for next batch
                write_buffer.clear();

                // Process ref_to_store and registrations in sorted group_id order
                for result in sorted_results {
                    // Store reference in global map
                    if let Some((group_id, ref_data)) = result.ref_to_store {
                        ref_segs.insert(group_id, ref_data);
                    }

                    // Register segments in collection
                    for reg in result.registrations {
                        if let Err(e) = coll.add_segment_placed(
                            &reg.sample_name,
                            &reg.contig_name,
                            reg.seg_part_no,
                            reg.group_id,
                            reg.in_group_id,
                            reg.is_rev_comp,
                            reg.raw_length,
                        ) {
                            eprintln!("Thread 0 error registering segment: {}", e);
                        }
                    }
                }
                // Locks released here when guards go out of scope
                drop(arch);
                drop(ref_segs);
                drop(coll);

                // Phase 4: Cleanup
                cleanup_batch_parallel(
                    &segment_groups,
                    &batch_local_groups,
                    &batch_local_terminators,
                    &map_segments,
                    &map_segments_terminators,
                    &parallel_state,
                    &config,
                );

                // Clear batch-local state after flush (start fresh for new sample)
                let mut samples = batch_samples.lock().unwrap();
                samples.clear();
            }

            // Barrier 4: All workers ready for next batch (reduced from 2 barriers)
            let barrier_start = std::time::Instant::now();
            barrier.wait();
            total_barrier_wait += barrier_start.elapsed();

            // Track total sync token processing time
            total_sync_processing += sync_start.elapsed();

            // Sync token processed - continue to next task
            continue;
        }

        // Start timing for segment processing
        let segment_start = std::time::Instant::now();
        contig_count += 1;

        // NOTE: Removed per-contig lock on batch_samples - not needed for deferred classification
        // The batch tracking is handled at the barrier level, not per-contig

        // Split into segments
        // Dynamic splitter discovery for non-reference contigs (matches C++ AGC find_new_splitters)
        // OPTIMIZATION: Cache reference sample name check at thread level to avoid lock per contig
        let is_reference_sample = task.sample_priority >= 1_000_000; // Reference sample has boosted priority

        let segments = if !is_reference_sample && !ref_singletons.is_empty() {
            // Non-reference contig with dynamic discovery enabled
            // Find NEW splitter k-mers unique to this contig (not in reference)
            // Position-based selection ensures only optimally-positioned k-mers become splitters
            let new_splitters = find_new_splitters_for_contig(
                &task.data,
                config.k,
                config.segment_size,
                &ref_singletons,
                &ref_duplicates,
            );

            // Combine base splitters with new splitters
            let mut combined_splitters = (*splitters).clone();
            combined_splitters.extend(new_splitters.iter());

            if config.verbosity > 2 && !new_splitters.is_empty() {
                eprintln!(
                    "DYNAMIC_SPLITTER: {} found {} new splitters for {} (total: {})",
                    task.sample_name, new_splitters.len(), task.contig_name, combined_splitters.len()
                );
            }

            split_at_splitters_with_size(&task.data, &combined_splitters, config.k, config.segment_size)
        } else {
            // Reference contig or dynamic discovery disabled - use base splitters only
            split_at_splitters_with_size(&task.data, &splitters, config.k, config.segment_size)
        };

        if config.verbosity > 2 {
            eprintln!(
                "Worker {} processing {} (split into {} segments)",
                worker_id,
                task.contig_name,
                segments.len()
            );
        }

        // NOTE: split_offsets and local_splits are no longer needed in the parallel loop
        // since classification (including splits) is deferred to the barrier.
        // This eliminates a lock acquisition per contig that was causing contention.

        // =================================================================
        // DEFERRED CLASSIFICATION: Buffer raw segments for parallel Phase 1
        // Classification is deferred to Thread 0 at the barrier to eliminate
        // lock contention from find_group_with_one_kmer and split logic.
        // =================================================================
        // OPTIMIZATION: Collect all segments for this contig locally, then push once
        // This reduces lock acquisitions from O(segments_per_contig) to O(1) per contig
        let contig_segments: Vec<RawBufferedSegment> = segments.iter().enumerate()
            .map(|(original_place, segment)| {
                // Precompute reverse complement (NO LOCKS - can run in parallel)
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

                RawBufferedSegment {
                    data: segment.data.clone(),
                    data_rc: segment_data_rc,
                    front_kmer: segment.front_kmer,
                    back_kmer: segment.back_kmer,
                    front_kmer_is_dir: segment.front_kmer_is_dir,
                    back_kmer_is_dir: segment.back_kmer_is_dir,
                    sample_name: task.sample_name.clone(),
                    contig_name: task.contig_name.clone(),
                    original_place,
                    sample_priority: task.sample_priority,
                }
            })
            .collect();

        // ONE lock acquisition for entire contig (reduces contention significantly)
        // Push to this worker's own buffer (NO CONTENTION - each worker has its own buffer)
        raw_segment_buffers[worker_id].lock().unwrap().extend(contig_segments);

        // End timing for segment processing
        total_segment_processing += segment_start.elapsed();

        // OLD CODE BELOW - REPLACED BY DEFERRED CLASSIFICATION
        // This block is preserved but commented out for reference during the transition.
        // The classification logic has been moved to classify_raw_segments_at_barrier().
        #[cfg(feature = "old_immediate_classification")]
        for (original_place, segment) in std::iter::empty::<(usize, &crate::segment::Segment)>() {
            // Calculate adjusted place based on prior splits in this contig
            // (matches C++ AGC lines 2033-2036: increment seg_part_no twice when split occurs)
            // OPTIMIZATION: Count splits before current position from both prior and local sets
            let prior_count = prior_contig_splits.range(..original_place).count();
            let local_count = local_splits.range(..original_place).count();
            let place = original_place + prior_count + local_count;

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
                    if crate::env_cache::debug_is_dir() {
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
            // For raw segments (both k-mers MISSING), use the same key for all
            // This matches C++ AGC: map_segments[make_pair(~0ull, ~0ull)] = 0
            // All raw segments share the same grouping key and will be assigned to the same group
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

            // NOTE: Split check must happen BEFORE creating BufferedSegment
            // to avoid moving segment_data prematurely
            // PERF: Lock deferred - normal path doesn't need segment_groups lock
            {
                // Phase 1: Check if group already exists
                // (matches C++ AGC: seg_map_mtx.lock() then find at line 1020)
                let key_exists = {
                    let seg_map = map_segments.read().unwrap();
                    seg_map.contains_key(&key)
                };

                // Phase 2: Try to split
                // C++ AGC only attempts splits when key doesn't exist (agc_compressor.cpp:1367)
                // This is the condition: p == map_segments.end() && both k-mers valid && both in terminators
                // Set RAGC_SPLIT_ALL=1 to try splitting even when key exists (experimental)
                // CRITICAL: C++ AGC lines 1374-1378 skip segment splitting when front == back!
                // When front == back, it just sets store_rc based on orientation, does NOT call
                // find_cand_segment_with_missing_middle_splitter. We must do the same.
                let split_allowed = if crate::env_cache::split_all() { true } else { !key_exists };

                // Debug: trace split decision
                if crate::env_cache::debug_split() && task.contig_name.contains("chrVII") && place >= 2 && place <= 5 {
                    eprintln!("RAGC_SPLIT_CHECK: contig={} seg={} key=({},{}) key_exists={} split_allowed={} front_missing={} back_missing={} front==back={}",
                        task.contig_name, place, key_front, key_back, key_exists, split_allowed,
                        key_front == MISSING_KMER, key_back == MISSING_KMER, key_front == key_back);
                }

                if split_allowed && key_front != MISSING_KMER && key_back != MISSING_KMER && key_front != key_back {
                    // CRITICAL: First attempt to find middle splitter
                    // Use ONLY global terminators (not batch-local) to match C++ AGC behavior
                    // C++ AGC only sees terminators from previous batches, not the current one
                    let middle_kmer_opt = {
                        let terminators = map_segments_terminators.read().unwrap();
                        let result = find_middle_splitter(key_front, key_back, &terminators);
                        // Debug: trace middle splitter result
                        if crate::env_cache::debug_split() && task.contig_name.contains("chrVII") && place >= 2 && place <= 5 {
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
                        if crate::env_cache::debug_split() {
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
                            let global_map = map_segments.read().unwrap();
                            (global_map.contains_key(&left_key), global_map.contains_key(&right_key))
                        };

                        // EXPERIMENTAL: Allow split even when groups don't exist
                        // Set RAGC_SPLIT_CREATE_GROUPS=1 to enable creating new groups during split
                        // This is needed for streaming mode where non-reference samples may create
                        // new segment groups that the reference sample didn't have.
                        let allow_create_groups = crate::env_cache::split_create_groups();

                        if !left_exists || !right_exists {
                            // Skip split - one or both target groups don't exist yet
                            // This matches C++ AGC behavior where .at() would throw
                            // UNLESS we're in experimental mode where we allow creating groups
                            if config.verbosity > 1 {
                                eprintln!(
                                    "SPLIT_SKIP_NO_GROUP: left_key=({},{}) exists={} right_key=({},{}) exists={} allow_create={}",
                                    left_key.kmer_front, left_key.kmer_back, left_exists,
                                    right_key.kmer_front, right_key.kmer_back, right_exists, allow_create_groups
                                );
                            }
                            if crate::env_cache::debug_split() {
                                eprintln!(
                                    "RAGC_SPLIT_SKIP_NO_GROUP: left=({},{}) exists={} right=({},{}) exists={} allow_create={}",
                                    left_key.kmer_front, left_key.kmer_back, left_exists,
                                    right_key.kmer_front, right_key.kmer_back, right_exists, allow_create_groups
                                );
                            }
                            if !allow_create_groups {
                                // Don't attempt split - fall through to normal segment processing
                            }
                        }

                        // Proceed with split if groups exist OR if we allow creating groups
                        if (left_exists && right_exists) || allow_create_groups {
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
                            allow_create_groups, // Force split at middle k-mer position if refs are empty
                        );

                        if let Some((left_data, right_data, _mid)) = split_result {
                            // PERF: Acquire lock only for split path (rare case)
                            // Normal segments bypass this entirely for better parallelism
                            let mut groups = segment_groups.lock().unwrap();

                            // FIX 27 v4: Compute separate orientations for left and right parts
                            // C++ AGC lines 1526-1536 and 1540-1550:
                            //   store_rc = (kmer_front.data() >= split_match.first)  -- for left
                            //   store2_rc = (split_match.first >= kmer_back.data())  -- for right
                            //
                            // When should_reverse=true, the segment was RC'd before splitting,
                            // so "left" in the split is from original RIGHT, and "right" is from original LEFT.
                            // We need to swap the k-mer comparisons accordingly.
                            let (left_should_rc, right_should_rc) = if should_reverse {
                                // Segment was RC'd: left is from original right, right is from original left
                                // Swap the k-mer associations
                                let left_should_rc = middle_kmer >= segment.back_kmer;  // use back_kmer for "left"
                                let right_should_rc = segment.front_kmer >= middle_kmer;  // use front_kmer for "right"
                                (left_should_rc, right_should_rc)
                            } else {
                                // Normal: left is from original left, right is from original right
                                let left_should_rc = segment.front_kmer >= middle_kmer;
                                let right_should_rc = middle_kmer >= segment.back_kmer;
                                (left_should_rc, right_should_rc)
                            };

                            // Transform data if needed: current state is `should_reverse`
                            // If target state differs, we RC the data
                            let left_data = if left_should_rc != should_reverse {
                                left_data.iter().rev().map(|&base| {
                                    match base {
                                        0 => 3, 1 => 2, 2 => 1, 3 => 0, _ => base
                                    }
                                }).collect::<Vec<u8>>()
                            } else {
                                left_data
                            };
                            let right_data = if right_should_rc != should_reverse {
                                right_data.iter().rev().map(|&base| {
                                    match base {
                                        0 => 3, 1 => 2, 2 => 1, 3 => 0, _ => base
                                    }
                                }).collect::<Vec<u8>>()
                            } else {
                                right_data
                            };

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
                                        // OPTIMIZATION: Read-check-write pattern to reduce lock contention
                                        // First check with read lock (fast path - most groups already exist)
                                        let group_id = {
                                            let global_map = map_segments.read().unwrap();
                                            if let Some(&existing_id) = global_map.get(&left_key) {
                                                existing_id
                                            } else {
                                                drop(global_map);
                                                // Group doesn't exist - upgrade to write lock
                                                let mut global_map = map_segments.write().unwrap();
                                                // Double-check after acquiring write lock (race condition)
                                                if let Some(&existing_id) = global_map.get(&left_key) {
                                                    existing_id
                                                } else {
                                                    // Create new group ID and register IMMEDIATELY to global map
                                                    let new_id = group_counter.fetch_add(1, Ordering::SeqCst);
                                                    global_map.insert(left_key.clone(), new_id);
                                                    drop(global_map);
                                                    // Also register to batch-local for flush tracking
                                                    let mut batch_map = batch_local_groups.lock().unwrap();
                                                    batch_map.insert(left_key.clone(), new_id);
                                                    new_id
                                                }
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
                                            let mut term_map = map_segments_terminators.write().unwrap();
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
                                    // FIX 27 v4: Use left_should_rc instead of should_reverse
                                    let (fixed_left_data, fixed_left_rc) = fix_orientation_for_group(&left_data, left_should_rc, &left_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let left_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: place, data: fixed_left_data, is_rev_comp: fixed_left_rc, sample_priority: task.sample_priority };
                                    left_buffer.segments.push(left_buffered);
                                    // Flush pack if full (matches C++ AGC write-as-you-go behavior)
                                    if left_buffer.should_flush_pack(config.pack_size) {
                                        flush_pack(left_buffer, &collection, &archive, &config, &reference_segments)
                                            .context("Failed to flush left pack")?;
                                    }
                                }
                                if !is_degenerate_right {
                                    let right_buffer = groups.entry(right_key.clone()).or_insert_with(|| {
                                        // OPTIMIZATION: Read-check-write pattern to reduce lock contention
                                        // First check with read lock (fast path - most groups already exist)
                                        let group_id = {
                                            let global_map = map_segments.read().unwrap();
                                            if let Some(&existing_id) = global_map.get(&right_key) {
                                                existing_id
                                            } else {
                                                drop(global_map);
                                                // Group doesn't exist - upgrade to write lock
                                                let mut global_map = map_segments.write().unwrap();
                                                // Double-check after acquiring write lock (race condition)
                                                if let Some(&existing_id) = global_map.get(&right_key) {
                                                    existing_id
                                                } else {
                                                    // Create new group ID and register IMMEDIATELY to global map
                                                    let new_id = group_counter.fetch_add(1, Ordering::SeqCst);
                                                    global_map.insert(right_key.clone(), new_id);
                                                    drop(global_map);
                                                    // Also register to batch-local for flush tracking
                                                    let mut batch_map = batch_local_groups.lock().unwrap();
                                                    batch_map.insert(right_key.clone(), new_id);
                                                    new_id
                                                }
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
                                            let mut term_map = map_segments_terminators.write().unwrap();
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
                                    // FIX 27 v4: Use right_should_rc instead of should_reverse
                                    let (fixed_right_data, fixed_right_rc) = fix_orientation_for_group(&right_data, right_should_rc, &right_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let right_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: seg_part, data: fixed_right_data, is_rev_comp: fixed_right_rc, sample_priority: task.sample_priority };
                                    right_buffer.segments.push(right_buffered);
                                    // Flush pack if full (matches C++ AGC write-as-you-go behavior)
                                    if right_buffer.should_flush_pack(config.pack_size) {
                                        flush_pack(right_buffer, &collection, &archive, &config, &reference_segments)
                                            .context("Failed to flush right pack")?;
                                    }
                                }
                            } else {
                                // reversed: right first
                                if !is_degenerate_right {
                                    let right_buffer = groups.entry(right_key.clone()).or_insert_with(|| {
                                        // BATCH-LOCAL: Check global first, then batch-local (group must exist from earlier)
                                        let group_id = {
                                            let global_map = map_segments.read().unwrap();
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
                                    // FIX 27 v4: Use right_should_rc instead of should_reverse
                                    let (fixed_right_data, fixed_right_rc) = fix_orientation_for_group(&right_data, right_should_rc, &right_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let right_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: place, data: fixed_right_data, is_rev_comp: fixed_right_rc, sample_priority: task.sample_priority };
                                    right_buffer.segments.push(right_buffered);
                                    // Flush pack if full (matches C++ AGC write-as-you-go behavior)
                                    if right_buffer.should_flush_pack(config.pack_size) {
                                        flush_pack(right_buffer, &collection, &archive, &config, &reference_segments)
                                            .context("Failed to flush right pack")?;
                                    }
                                }
                                if !is_degenerate_left {
                                    let left_buffer = groups.entry(left_key.clone()).or_insert_with(|| {
                                        // BATCH-LOCAL: Check global first, then batch-local (group must exist from earlier)
                                        let group_id = {
                                            let global_map = map_segments.read().unwrap();
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
                                    // FIX 27 v4: Use left_should_rc instead of should_reverse
                                    let (fixed_left_data, fixed_left_rc) = fix_orientation_for_group(&left_data, left_should_rc, &left_key, &map_segments, &batch_local_groups, &reference_orientations);
                                    let left_buffered = BufferedSegment { sample_name: task.sample_name.clone(), contig_name: task.contig_name.clone(), seg_part_no: seg_part, data: fixed_left_data, is_rev_comp: fixed_left_rc, sample_priority: task.sample_priority };
                                    left_buffer.segments.push(left_buffered);
                                    // Flush pack if full (matches C++ AGC write-as-you-go behavior)
                                    if left_buffer.should_flush_pack(config.pack_size) {
                                        flush_pack(left_buffer, &collection, &archive, &config, &reference_segments)
                                            .context("Failed to flush left pack")?;
                                    }
                                }
                            }

                            // Optional: assert lengths vs C++ archive if provided
                            if let Some(assert_path) = crate::env_cache::assert_cpp_archive() {
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
                                                    if crate::env_cache::assert_verbose() {
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
                                // OPTIMIZATION: Track locally for this task AND globally for other workers
                                local_splits.insert(original_place);
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
                        let seg_map = map_segments.read().unwrap();
                        seg_map.contains_key(&key)
                    };

                    // Debug: count how many segments could be eligible for secondary fallback
                    if crate::env_cache::debug_fallback2_enabled() {
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

                        if crate::env_cache::debug_fallback2_enabled() {
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
                                let seg_map = map_segments.read().unwrap();
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

                // FIX 18: Do NOT adjust orientation to match reference - C++ AGC stores each segment
                // with its own computed is_rev_comp based on k-mer comparison. Segments in the same
                // group can have different is_rev_comp values.
                let (final_should_reverse, final_segment_data) = (should_reverse, segment_data);

                if config.verbosity > 2 {
                    eprintln!("DEFER_SEGMENT: front={} back={} sample={} contig={} place={}",
                        key_front, key_back, task.sample_name, task.contig_name, place);
                }

                // PHASE 1 (PARALLEL): Add segment to buffered_seg_part
                // Check if group exists (brief read lock on map_segments)
                let group_id_opt = {
                    let seg_map = map_segments.read().unwrap();
                    seg_map.get(&key).copied()
                };

                if let Some(group_id) = group_id_opt {
                    // KNOWN: add to per-group buffer (per-group lock only - PARALLEL)
                    buffered_seg_part.add_known(group_id, BufferedSegment {
                        sample_name: task.sample_name.clone(),
                        contig_name: task.contig_name.clone(),
                        seg_part_no: place,
                        data: final_segment_data,
                        is_rev_comp: final_should_reverse,
                        sample_priority: task.sample_priority,
                    });
                } else {
                    // NEW: add to s_seg_part (brief global lock on BTreeSet)
                    buffered_seg_part.add_new(NewSegment {
                        kmer_front: key.kmer_front,
                        kmer_back: key.kmer_back,
                        sample_priority: task.sample_priority,
                        sample_name: task.sample_name.clone(),
                        contig_name: task.contig_name.clone(),
                        seg_part_no: place,
                        data: final_segment_data,
                        should_reverse: final_should_reverse,
                    });
                }

                // Segment will be handled in flush_batch at barrier synchronization point
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
    terminators: &BTreeMap<u64, Vec<u64>>,
) -> Option<u64> {
    let front_connections = terminators.get(&front_kmer)?;
    let back_connections = terminators.get(&back_kmer)?;

    #[cfg(feature = "cpp_agc")]
    {
        if let Some(m) = crate::ragc_ffi::find_middle(front_connections, back_connections) {
            return Some(m);
        }
        if crate::env_cache::debug_split_find() {
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
        if crate::env_cache::debug_split_find() {
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
    map_segments: &Arc<RwLock<BTreeMap<SegmentGroupKey, u32>>>,
    map_segments_terminators: &Arc<RwLock<BTreeMap<u64, Vec<u64>>>>,
    reference_segments: &Arc<RwLock<BTreeMap<u32, Vec<u8>>>>,
    config: &StreamingQueueConfig,
    should_reverse: bool,
    force_split_on_empty_refs: bool, // When true, split at middle k-mer position even if FFI says no
) -> Option<(Vec<u8>, Vec<u8>, u64)> {
    if config.verbosity > 1 {
        eprintln!(
            "SPLIT_ATTEMPT: front={} back={} middle={}",
            front_kmer, back_kmer, middle_kmer
        );
    }

    // Debug: trace split attempt
    if crate::env_cache::debug_split() {
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
        let map_segments_locked = map_segments.read().unwrap();
        let ref_segments_locked = reference_segments.read().unwrap();

        // C++ AGC uses map_segments[key] which returns 0 (default) if key doesn't exist.
        // v_segments[0] is a raw group initialized with empty_ctg = { 0x7f } (1 byte).
        // This gives maximum LZ cost (no compression) for non-existent groups.
        // To match C++ AGC behavior, use the actual reference if available,
        // otherwise use empty reference (gives max cost like C++ AGC's v_segments[0]).
        let segment_id = map_segments_locked.get(key).copied();

        if let Some(ref_data) = segment_id.and_then(|id| ref_segments_locked.get(&id)) {
            // Reference exists! Prepare LZDiff on-demand
            if crate::env_cache::debug_split_ref() {
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
            if crate::env_cache::debug_split_ref() {
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
    let maybe_best: Option<(usize, usize)> = None; // (best_pos, seg2_start)
    #[cfg(feature = "cpp_agc")]
    {
        // Inspect availability of left/right references and log keys
        let (left_seg_id_opt, right_seg_id_opt) = {
            let map_segments_locked = map_segments.read().unwrap();
            (map_segments_locked.get(left_key).copied(), map_segments_locked.get(right_key).copied())
        };
        let (left_have_ref, right_have_ref) = {
            let ref_segments_locked = reference_segments.read().unwrap();
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
            let term_map = map_segments_terminators.read().unwrap();
            (term_map.get(&front_kmer).cloned().unwrap_or_default(),
             term_map.get(&back_kmer).cloned().unwrap_or_default())
        };

        // Always attempt FFI decision; if refs are missing, C++ will decide no-split
        let (ref_left_opt, ref_right_opt) = {
            let ref_segments_locked = reference_segments.read().unwrap();
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
            if !has_mid { return None; }

            // FFI found middle k-mer but may have said !should due to empty refs
            if should {
                maybe_best = Some((bp, s2));
            } else if force_split_on_empty_refs && ref_left.is_empty() && ref_right.is_empty() {
                // FALLBACK: FFI can't compute costs because refs are empty, but we want to
                // create new groups. Search for ANY terminator k-mer in the segment that can serve as a split point.
                // This handles the case where the exact middle_kmer from reference has a mutation in this sample.
                if config.verbosity > 1 {
                    eprintln!("SPLIT_FALLBACK: FFI said no but force_split_on_empty_refs=true, searching for any terminator k-mer in segment");
                }

                // Build a set of potential middle k-mers from both neighbor lists
                let mut potential_middles: AHashSet<u64> = AHashSet::new();
                for &kmer in front_neighbors.iter() {
                    if kmer != MISSING_KMER && kmer != front_kmer && kmer != back_kmer {
                        potential_middles.insert(kmer);
                    }
                }
                for &kmer in back_neighbors.iter() {
                    if kmer != MISSING_KMER && kmer != front_kmer && kmer != back_kmer {
                        potential_middles.insert(kmer);
                    }
                }

                if config.verbosity > 1 {
                    eprintln!("SPLIT_FALLBACK: {} potential middle k-mers from terminators", potential_middles.len());
                    for &pm in potential_middles.iter().take(5) {
                        eprintln!("  potential_middle: {:#x}", pm);
                    }
                }

                // Search for ANY terminator k-mer in the segment
                let k = config.k;
                if segment_dir.len() >= k && !potential_middles.is_empty() {
                    let mut found_pos: Option<(usize, u64)> = None; // (pos, kmer)
                    let mut kmer_obj = crate::kmer::Kmer::new(k as u32, crate::kmer::KmerMode::Canonical);
                    for (i, &base) in segment_dir.iter().enumerate() {
                        if base > 3 {
                            kmer_obj.reset();
                        } else {
                            kmer_obj.insert(base as u64);
                            if kmer_obj.is_full() {
                                let kmer_at_pos = kmer_obj.data();
                                let pos = i + 1 - k; // Position of k-mer start
                                // Check if this k-mer is in our set of potential middles
                                if potential_middles.contains(&kmer_at_pos) {
                                    // Ensure we're not at the very beginning or end
                                    if pos > k && pos + k + k < segment_dir.len() {
                                        found_pos = Some((pos, kmer_at_pos));
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    if let Some((pos, found_kmer)) = found_pos {
                        // Split at position just after the found k-mer
                        let split_pos = pos + k;
                        if split_pos > k + 1 && split_pos + k + 1 < segment_dir.len() {
                            if config.verbosity > 1 {
                                eprintln!("SPLIT_FALLBACK_FOUND: terminator kmer={:#x} found at pos={}, splitting at {}", found_kmer, pos, split_pos);
                            }
                            maybe_best = Some((split_pos, split_pos));
                        } else if config.verbosity > 1 {
                            eprintln!("SPLIT_FALLBACK_DEGENERATE: pos={} split_pos={} segment_len={}", pos, split_pos, segment_dir.len());
                        }
                    } else {
                        // FALLBACK 2: Terminators not found - discover a NEW singleton k-mer in the segment
                        // Similar to C++ AGC's find_new_splitters() but simpler: just find any singleton
                        if config.verbosity > 1 {
                            eprintln!("SPLIT_FALLBACK_DISCOVER: trying to find singleton k-mer in segment (len={})", segment_dir.len());
                        }

                        // Collect all k-mers in the middle region of the segment
                        let min_margin = k * 2; // Don't split too close to edges
                        let search_start = min_margin;
                        let search_end = segment_dir.len().saturating_sub(min_margin);

                        if search_end > search_start + k {
                            // Enumerate k-mers and find singletons
                            let mut kmer_positions: Vec<(u64, usize)> = Vec::new();
                            let mut kmer_obj2 = crate::kmer::Kmer::new(k as u32, crate::kmer::KmerMode::Canonical);

                            for (i, &base) in segment_dir[search_start..search_end].iter().enumerate() {
                                if base > 3 {
                                    kmer_obj2.reset();
                                } else {
                                    kmer_obj2.insert(base as u64);
                                    if kmer_obj2.is_full() {
                                        let kmer_val = kmer_obj2.data();
                                        let pos = search_start + i + 1 - k;
                                        kmer_positions.push((kmer_val, pos));
                                    }
                                }
                            }

                            // Sort by k-mer value to find duplicates
                            kmer_positions.sort_by_key(|&(kmer, _)| kmer);

                            // Find first singleton (k-mer that appears exactly once)
                            let mut singleton_pos: Option<usize> = None;
                            let mut i = 0;
                            while i < kmer_positions.len() {
                                let (kmer, pos) = kmer_positions[i];
                                let mut j = i + 1;
                                while j < kmer_positions.len() && kmer_positions[j].0 == kmer {
                                    j += 1;
                                }
                                // If exactly one occurrence, it's a singleton
                                if j == i + 1 {
                                    singleton_pos = Some(pos);
                                    if config.verbosity > 1 {
                                        eprintln!("SPLIT_FALLBACK_SINGLETON: found singleton kmer={:#x} at pos={}", kmer, pos);
                                    }
                                    break;
                                }
                                i = j;
                            }

                            if let Some(pos) = singleton_pos {
                                let split_pos = pos + k;
                                if config.verbosity > 1 {
                                    eprintln!("SPLIT_FALLBACK_SINGLETON_SPLIT: splitting at {}", split_pos);
                                }
                                maybe_best = Some((split_pos, split_pos));
                            } else if config.verbosity > 1 {
                                eprintln!("SPLIT_FALLBACK_NO_SINGLETON: no singleton k-mers found in middle region");
                            }
                        } else if config.verbosity > 1 {
                            eprintln!("SPLIT_FALLBACK_TOO_SHORT: segment too short for singleton search");
                        }
                    }
                }
            } else {
                return None;
            }
        } else if config.verbosity > 1 {
            eprintln!("FFI_DECIDE: unavailable (decide_split returned None)");
        }
    }

    // If FFI provided best position, use it; otherwise compute costs in Rust
    let mut v_costs1 = if maybe_best.is_none() {
        if let Some(lz_left) = prepare_on_demand(left_key, "left") {
        #[cfg(feature = "cpp_agc")]
        {
            // Unused path when FFI returns best split; kept for completeness
            let ref_left = {
                let map_segments_locked = map_segments.read().unwrap();
                let ref_segments_locked = reference_segments.read().unwrap();
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
        if crate::env_cache::debug_split() {
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

    let v_costs2 = if maybe_best.is_none() {
        if let Some(lz_right) = prepare_on_demand(right_key, "right") {
        #[cfg(feature = "cpp_agc")]
        {
            let ref_right = {
                let map_segments_locked = map_segments.read().unwrap();
                let ref_segments_locked = reference_segments.read().unwrap();
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
    if crate::env_cache::debug_split_map() && maybe_best.is_none() {
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

    // Apply degenerate position rules ALWAYS to prevent tiny segments at boundaries.
    // Even when FFI or fallback paths provide best_pos, we must enforce this constraint
    // to match C++ AGC behavior (agc_compressor.cpp:1685-1688).
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
            "BOUNDARY_CLAMP: original_best_pos={} clamped_to={} (len={}, k+1={}) source={}",
            original_best_pos, best_pos, v_costs1.len(), k + 1,
            if maybe_best.is_some() { "FFI/fallback" } else { "cost_calc" }
        );
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
