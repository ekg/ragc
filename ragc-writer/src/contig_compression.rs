//! Contig compression with inline segmentation
//!
//! This module implements C++ AGC's compress_contig and add_segment logic
//! for inline segmentation during compression.
//!
//! See docs/INLINE_SEGMENTATION_PATTERN.md for detailed C++ AGC analysis.

use crate::kmer::{Kmer, KmerMode};
use crate::segment_buffer::BufferedSegments;
use std::sync::{Arc, Mutex};

/// Contig type (vector of nucleotides in numeric encoding)
pub type Contig = Vec<u8>;

/// Missing k-mer marker (matches C++ AGC's kmer_t(-1) = u64::MAX)
pub const MISSING_KMER: u64 = u64::MAX;

/// Segment part for inline buffering
#[derive(Debug, Clone)]
pub struct SegmentPart {
    pub kmer1: u64,
    pub kmer2: u64,
    pub sample_name: String,
    pub contig_name: String,
    pub seg_data: Vec<u8>,
    pub is_rev_comp: bool,
    pub seg_part_no: u32,
}

/// Shared state for compression (passed from worker threads)
pub struct CompressionContext {
    pub splitters: Arc<Mutex<std::collections::HashSet<u64>>>,
    pub bloom_splitters: Arc<Mutex<crate::bloom_filter::BloomFilter>>,
    pub buffered_segments: Arc<Mutex<BufferedSegments>>,
    pub kmer_length: usize,
    pub adaptive_mode: bool,
    pub map_segments: Arc<Mutex<std::collections::HashMap<(u64, u64), u32>>>,
    pub map_segments_terminators: Arc<Mutex<std::collections::HashMap<u64, Vec<u64>>>>,
    pub concatenated_genomes: bool,
}

/// Compress a contig with inline segmentation
///
/// Matches C++ AGC's compress_contig (agc_compressor.cpp:2000-2054)
///
/// # Arguments
/// * `sample_name` - Sample identifier
/// * `contig_name` - Contig identifier
/// * `contig` - Nucleotide sequence (numeric encoding: A=0, C=1, G=2, T=3)
/// * `ctx` - Compression context with splitters and buffering
///
/// # Returns
/// * `true` - Contig segmented successfully
/// * `false` - Contig failed to segment (adaptive mode: needs new splitters)
///
/// # Algorithm
/// 1. Scan contig with k-mer sliding window
/// 2. Check each k-mer against splitters (bloom filter + hash set)
/// 3. When splitter found: extract segment, call add_segment()
/// 4. If no splitters found and adaptive mode: return false
/// 5. Handle final segment (after last splitter to end)
pub fn compress_contig(
    sample_name: &str,
    contig_name: &str,
    contig: &Contig,
    ctx: &CompressionContext,
) -> bool {
    let kmer_length = ctx.kmer_length as u32;

    // Create k-mer scanner for canonical k-mers
    let mut kmer = Kmer::new(kmer_length, KmerMode::Canonical);

    // Track segment boundaries
    let mut split_pos: usize = 0; // Start of current segment
    let mut split_kmer = Kmer::new(kmer_length, KmerMode::Canonical); // K-mer at split position
    let mut seg_part_no: u32 = 0;

    // Scan contig for splitters
    for (pos, &base) in contig.iter().enumerate() {
        kmer.insert(base as u64);

        if !kmer.is_full() {
            continue;
        }

        // Get canonical k-mer value
        let kmer_val = kmer.data_canonical();

        // Check if this k-mer is a splitter
        // C++ AGC: bloom_splitters.check(d) && hs_splitters.check(d)
        if is_splitter(kmer_val, ctx) {
            // Found splitter - extract segment from split_pos to current position
            let seg_start = split_pos;
            let seg_end = pos + 1; // Inclusive of current position
            let segment = contig[seg_start..seg_end].to_vec();

            // Add segment with terminal k-mers (split_kmer, current kmer)
            add_segment(
                sample_name,
                contig_name,
                seg_part_no,
                segment,
                &split_kmer,
                &kmer,
                ctx,
            );

            // Update for next segment
            seg_part_no += 1;
            split_pos = pos + 1 - kmer_length as usize; // Start of next segment overlaps by kmer_length
            split_kmer = kmer.clone();
        }
    }

    // Check if contig failed to segment (adaptive mode)
    if ctx.adaptive_mode && seg_part_no == 0 {
        // No splitters found - this is a "hard" contig
        // split_kmer is still empty (never set), which C++ AGC checks as:
        // if (adaptive_compression && split_kmer == CKmer(...))
        return false; // Signal failure - needs new splitters
    }

    // Add final segment (from last splitter to end of contig)
    if split_pos < contig.len() {
        let segment = contig[split_pos..].to_vec();

        add_segment(
            sample_name,
            contig_name,
            seg_part_no,
            segment,
            &split_kmer,
            &Kmer::new(kmer_length, KmerMode::Canonical), // Empty k-mer for end
            ctx,
        );
    }

    true // Success
}

/// Add segment to buffered storage with terminal k-mer detection
///
/// Matches C++ AGC's add_segment (agc_compressor.cpp:1275-1507)
///
/// # Terminal K-mer Cases
/// 1. **No terminators** (k1 = k2 = ~0): Whole-contig segment
/// 2. **Both terminators** (k1, k2 both valid): Normal segment
/// 3. **Front-only** (k1 valid, k2 = ~0): First segment
/// 4. **Back-only** (k1 = ~0, k2 valid): Last segment
///
/// # Split Detection
/// If segment is NEW and both k1 and k2 are in map_segments_terminators:
/// - Try to find shared middle splitter k_mid
/// - Split segment into two overlapping parts
///
/// See INLINE_SEGMENTATION_PATTERN.md Section 2 for detailed algorithm.
fn add_segment(
    sample_name: &str,
    contig_name: &str,
    seg_part_no: u32,
    mut segment: Vec<u8>,
    kmer_front: &Kmer,
    kmer_back: &Kmer,
    ctx: &CompressionContext,
) {
    // Determine terminal k-mers (k1, k2) based on which are full
    let k1 = if kmer_front.is_full() {
        kmer_front.data_canonical()
    } else {
        MISSING_KMER
    };

    let k2 = if kmer_back.is_full() {
        kmer_back.data_canonical()
    } else {
        MISSING_KMER
    };

    // Normalize key to canonical order (C++ AGC uses minmax)
    let pk = if k1 <= k2 { (k1, k2) } else { (k2, k1) };

    // Check if this segment key is already registered
    let map_segments = ctx.map_segments.lock().unwrap();
    let group_id = map_segments.get(&pk).copied();
    drop(map_segments);

    // Try split detection if NEW segment
    let mut segment2: Option<Vec<u8>> = None;
    let mut pk2: Option<(u64, u64)> = None;

    if group_id.is_none() && !ctx.concatenated_genomes {
        // NEW segment - check for split eligibility
        if k1 != MISSING_KMER && k2 != MISSING_KMER {
            // Both terminators exist - check if they're in terminators map
            let terminators = ctx.map_segments_terminators.lock().unwrap();
            let k1_in_map = terminators.contains_key(&k1);
            let k2_in_map = terminators.contains_key(&k2);
            drop(terminators);

            if k1_in_map && k2_in_map {
                // Try to find shared middle splitter
                if let Some((k_middle, left_size, right_size)) =
                    find_cand_segment_with_missing_middle_splitter(kmer_front, kmer_back, ctx)
                {
                    // Determine split outcome
                    if left_size == 0 || right_size == 0 {
                        // No actual split - entire segment goes to one group
                        // (This shouldn't happen but C++ AGC handles it)
                    } else {
                        // Split into 2 segments with overlap
                        let kmer_length = ctx.kmer_length;
                        let seg2_start_pos = left_size - kmer_length / 2;

                        // segment2: [seg2_start_pos, end)
                        segment2 = Some(segment[seg2_start_pos..].to_vec());

                        // segment: [0, seg2_start_pos + kmer_length)
                        segment.resize(seg2_start_pos + kmer_length, 0);

                        // Segment 1 key: (k1, k_middle)
                        pk2 = Some(if k_middle <= k2 {
                            (k_middle, k2)
                        } else {
                            (k2, k_middle)
                        });

                        // Note: pk remains (k1, k_middle) or (k_middle, k1) depending on order
                    }
                }
            }
        }
    }

    // TODO: ZSTD compress segments
    // For now, just store uncompressed
    let compressed_seg = segment.clone();

    // Add first segment (or whole segment if no split)
    if group_id.is_some() {
        // KNOWN segment
        let buffered = ctx.buffered_segments.lock().unwrap();
        buffered.add_known(
            group_id.unwrap(),
            MISSING_KMER,
            MISSING_KMER,
            sample_name.to_string(),
            contig_name.to_string(),
            compressed_seg,
            false, // TODO: Determine is_rev_comp from k-mer comparison
            seg_part_no,
        );
    } else {
        // NEW segment
        let buffered = ctx.buffered_segments.lock().unwrap();
        buffered.add_new(
            pk.0,
            pk.1,
            sample_name.to_string(),
            contig_name.to_string(),
            compressed_seg,
            false, // TODO: Determine is_rev_comp
            seg_part_no,
        );
    }

    // Add second segment if split occurred
    if let (Some(seg2), Some(pk2_key)) = (segment2, pk2) {
        // TODO: ZSTD compress segment2
        let compressed_seg2 = seg2.clone();

        // Check if second segment key is known
        let map_segments = ctx.map_segments.lock().unwrap();
        let group_id2 = map_segments.get(&pk2_key).copied();
        drop(map_segments);

        let buffered = ctx.buffered_segments.lock().unwrap();
        if let Some(gid2) = group_id2 {
            // KNOWN segment
            buffered.add_known(
                gid2,
                MISSING_KMER,
                MISSING_KMER,
                sample_name.to_string(),
                contig_name.to_string(),
                compressed_seg2,
                false, // TODO: is_rev_comp
                seg_part_no + 1,
            );
        } else {
            // NEW segment
            buffered.add_new(
                pk2_key.0,
                pk2_key.1,
                sample_name.to_string(),
                contig_name.to_string(),
                compressed_seg2,
                false, // TODO: is_rev_comp
                seg_part_no + 1,
            );
        }
    }
}

/// Check if k-mer is a splitter
///
/// Matches C++ AGC's two-stage check:
/// 1. Bloom filter (fast, probabilistic)
/// 2. Hash set (exact)
fn is_splitter(kmer: u64, ctx: &CompressionContext) -> bool {
    // Fast check: bloom filter (may have false positives)
    let bloom = ctx.bloom_splitters.lock().unwrap();
    if !bloom.check(kmer) {
        return false; // Definitely not a splitter
    }
    drop(bloom); // Release lock before next check

    // Exact check: hash set (no false positives)
    let splitters = ctx.splitters.lock().unwrap();
    splitters.contains(&kmer)
}

/// Find candidate segment with missing middle splitter (for split detection)
///
/// Matches C++ AGC's find_cand_segment_with_missing_middle_splitter
/// (agc_compressor.cpp:1510-1597)
///
/// # Algorithm
/// 1. Find intersection of terminators for k1 and k2
/// 2. For each shared k_mid in intersection:
///    - Check if (k1, k_mid) and (k_mid, k2) both exist in map_segments
///    - Find k_mid position in segment (via find_middle_splitter)
/// 3. Return first valid split: (middle_kmer, left_size, right_size)
///
/// # Returns
/// - Some((k_middle, left_size, right_size)) if split found
/// - None if no valid split
fn find_cand_segment_with_missing_middle_splitter(
    kmer_front: &Kmer,
    kmer_back: &Kmer,
    ctx: &CompressionContext,
) -> Option<(u64, usize, usize)> {
    let k1 = kmer_front.data_canonical();
    let k2 = kmer_back.data_canonical();

    // Get terminator lists for k1 and k2
    let terminators = ctx.map_segments_terminators.lock().unwrap();
    let k1_terminators = terminators.get(&k1)?;
    let k2_terminators = terminators.get(&k2)?;

    // Find shared terminators (intersection)
    // Both lists are sorted, so we can use two-pointer technique
    let mut shared_splitters = Vec::new();
    let mut i = 0;
    let mut j = 0;

    while i < k1_terminators.len() && j < k2_terminators.len() {
        if k1_terminators[i] == k2_terminators[j] {
            shared_splitters.push(k1_terminators[i]);
            i += 1;
            j += 1;
        } else if k1_terminators[i] < k2_terminators[j] {
            i += 1;
        } else {
            j += 1;
        }
    }

    if shared_splitters.is_empty() {
        return None;
    }

    drop(terminators);

    // Try each shared splitter
    let map_segments = ctx.map_segments.lock().unwrap();

    for &k_middle in &shared_splitters {
        // Create keys (k1, k_middle) and (k_middle, k2) in canonical order
        let pk_left = if k1 <= k_middle {
            (k1, k_middle)
        } else {
            (k_middle, k1)
        };

        let pk_right = if k_middle <= k2 {
            (k_middle, k2)
        } else {
            (k2, k_middle)
        };

        // Check if both groups exist
        if map_segments.contains_key(&pk_left) && map_segments.contains_key(&pk_right) {
            drop(map_segments);

            // TODO: Find actual position of k_middle in segment using find_middle_splitter
            // This requires:
            // - Decompressing both segments from v_segments[pk_left] and v_segments[pk_right]
            // - Finding the position where k_middle occurs
            // - Computing optimal split position via cost analysis
            //
            // For now, return a simplified result with approximate position
            // This is a placeholder - real implementation needed for correctness

            // Placeholder: assume k_middle is in the middle
            // TODO: Implement find_middle_splitter properly
            let left_size = ctx.kmer_length; // Placeholder
            let right_size = ctx.kmer_length; // Placeholder

            return Some((k_middle, left_size, right_size));
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_missing_kmer_constant() {
        // MISSING_KMER must be u64::MAX to match C++ AGC's kmer_t(-1)
        assert_eq!(MISSING_KMER, u64::MAX);
    }

    #[test]
    fn test_segment_part_creation() {
        let part = SegmentPart {
            kmer1: 12345,
            kmer2: 67890,
            sample_name: "sample1".to_string(),
            contig_name: "chr1".to_string(),
            seg_data: vec![0, 1, 2, 3],
            is_rev_comp: false,
            seg_part_no: 0,
        };

        assert_eq!(part.kmer1, 12345);
        assert_eq!(part.kmer2, 67890);
        assert_eq!(part.sample_name, "sample1");
        assert_eq!(part.seg_part_no, 0);
    }

    #[test]
    fn test_compress_contig_no_splitters() {
        use std::collections::HashSet;

        let ctx = CompressionContext {
            splitters: Arc::new(Mutex::new(HashSet::new())),
            bloom_splitters: Arc::new(Mutex::new(crate::bloom_filter::BloomFilter::new(1024))),
            buffered_segments: Arc::new(Mutex::new(BufferedSegments::new(0))),
            kmer_length: 21,
            adaptive_mode: false,
            map_segments: Arc::new(Mutex::new(HashMap::new())),
            map_segments_terminators: Arc::new(Mutex::new(HashMap::new())),
            concatenated_genomes: false,
        };

        // Simple contig with no splitters
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT

        let result = compress_contig("sample1", "chr1", &contig, &ctx);

        // Should succeed (one whole-contig segment)
        assert!(result);

        // Should have 1 NEW segment
        let buffered = ctx.buffered_segments.lock().unwrap();
        assert_eq!(buffered.get_num_new(), 1);
    }

    #[test]
    fn test_compress_contig_adaptive_failure() {
        use std::collections::HashSet;

        let ctx = CompressionContext {
            splitters: Arc::new(Mutex::new(HashSet::new())),
            bloom_splitters: Arc::new(Mutex::new(crate::bloom_filter::BloomFilter::new(1024))),
            buffered_segments: Arc::new(Mutex::new(BufferedSegments::new(0))),
            kmer_length: 21,
            adaptive_mode: true, // Adaptive mode enabled
            map_segments: Arc::new(Mutex::new(HashMap::new())),
            map_segments_terminators: Arc::new(Mutex::new(HashMap::new())),
            concatenated_genomes: false,
        };

        // Short contig with no splitters
        let contig = vec![0, 1, 2, 3, 0, 1, 2, 3];

        let result = compress_contig("sample1", "chr1", &contig, &ctx);

        // Should FAIL in adaptive mode (no splitters found)
        assert!(!result);
    }
}
