// AGC Decompressor
// Extracts genomes from AGC archives

use crate::{
    genome_io::{GenomeWriter, CNV_NUM},
    lz_diff::LZDiff,
    segment_compression::decompress_segment_with_marker,
};

use anyhow::{anyhow, Context, Result};
use ragc_common::{
    stream_delta_name, stream_ref_name, Archive, CollectionV3, Contig, SegmentDesc, AGC_FILE_MAJOR,
    AGC_FILE_MINOR, CONTIG_SEPARATOR,
};
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

/// Configuration for the decompressor
#[derive(Debug, Clone)]
pub struct DecompressorConfig {
    pub verbosity: u32,
}

impl Default for DecompressorConfig {
    fn default() -> Self {
        DecompressorConfig { verbosity: 1 }
    }
}

/// AGC Decompressor
///
/// Provides thread-safe read access to AGC archives.
/// For concurrent access, use `clone_for_thread()` to create independent readers.
///
/// # Thread Safety
/// - `Decompressor` is NOT `Sync` - cannot be shared between threads
/// - Use `clone_for_thread()` to create independent readers (cheap - shares archive data)
/// - Each clone can be used from a separate thread independently
///
/// # Example
/// ```no_run
/// use ragc_core::{Decompressor, DecompressorConfig};
/// use std::thread;
///
/// let dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
/// let samples = dec.list_samples();
///
/// let handles: Vec<_> = samples.into_iter().map(|sample_name| {
///     let mut thread_dec = dec.clone_for_thread().unwrap();
///     thread::spawn(move || {
///         thread_dec.get_sample(&sample_name)
///     })
/// }).collect();
/// # Ok::<(), anyhow::Error>(())
/// ```
pub struct Decompressor {
    config: DecompressorConfig,
    archive: Archive,
    collection: CollectionV3,

    // Cached segment data (group_id -> reference segment)
    segment_cache: HashMap<u32, Contig>,

    // Archive parameters
    _segment_size: u32,
    pub kmer_length: u32,
    pub min_match_len: u32,

    // Archive path for cloning
    archive_path: String,
}

impl Decompressor {
    /// Open an existing archive for decompression
    pub fn open(archive_path: &str, config: DecompressorConfig) -> Result<Self> {
        let mut archive = Archive::new_reader();
        archive
            .open(archive_path)
            .context("Failed to open archive for reading")?;

        let mut collection = CollectionV3::new();

        // Load segment_size, kmer_length, and min_match_len from params stream
        let (segment_size, kmer_length, min_match_len) = Self::load_params(&mut archive)?;

        if config.verbosity > 1 {
            eprintln!("Loaded params: segment_size={segment_size}, kmer_length={kmer_length}");
        }

        // Configure collection
        collection.set_config(segment_size, kmer_length, None);

        // Load collection metadata
        collection.prepare_for_decompression(&archive)?;
        collection.load_batch_sample_names(&mut archive)?;

        if config.verbosity > 0 {
            eprintln!(
                "Loaded archive with {} samples",
                collection.get_no_samples()
            );
            let samples = collection.get_samples_list(false);
            eprintln!("Sample names: {samples:?}");
        }

        Ok(Decompressor {
            config,
            archive,
            collection,
            segment_cache: HashMap::new(),
            _segment_size: segment_size,
            kmer_length,
            min_match_len,
            archive_path: archive_path.to_string(),
        })
    }

    /// Load archive parameters from the params stream
    ///
    /// Supports multiple params formats for C++ AGC compatibility:
    /// - 12 bytes: kmer_length, min_match_len, pack_cardinality (older C++ AGC)
    /// - 16 bytes: kmer_length, min_match_len, pack_cardinality, segment_size (newer C++ AGC)
    /// - 20 bytes: kmer_length, min_match_len, pack_cardinality, segment_size, no_raw_groups (ragc)
    fn load_params(archive: &mut Archive) -> Result<(u32, u32, u32)> {
        // Get params stream
        let stream_id = archive
            .get_stream_id("params")
            .ok_or_else(|| anyhow!("params stream not found in archive"))?;

        // Check that there is exactly one part
        let num_parts = archive.get_num_parts(stream_id);
        if num_parts != 1 {
            anyhow::bail!("Expected 1 part in params stream, found {num_parts}");
        }

        // Read the params part
        let (data, _metadata) = archive.get_part_by_id(stream_id, 0)?;

        // Parse based on length
        if data.len() < 12 {
            anyhow::bail!(
                "params stream too short: {} bytes (expected at least 12)",
                data.len()
            );
        }

        let kmer_length = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
        let min_match_len = u32::from_le_bytes([data[4], data[5], data[6], data[7]]);
        let _pack_cardinality = u32::from_le_bytes([data[8], data[9], data[10], data[11]]);

        // segment_size is optional (not present in older C++ AGC archives)
        let segment_size = if data.len() >= 16 {
            u32::from_le_bytes([data[12], data[13], data[14], data[15]])
        } else {
            // Default segment_size for older archives (matches C++ AGC default)
            60000
        };

        // no_raw_groups is optional (only in ragc archives) - ignored if present

        Ok((segment_size, kmer_length, min_match_len))
    }

    /// Get list of samples in the archive
    pub fn list_samples(&self) -> Vec<String> {
        self.collection.get_samples_list(false)
    }

    /// Get compression statistics for all streams
    /// Returns: Vec<(stream_name, raw_size, packed_size, num_parts)>
    pub fn get_compression_stats(&self) -> Vec<(String, u64, u64, usize)> {
        let mut stats = Vec::new();
        for stream_id in 0..self.archive.get_num_streams() {
            if let Some(name) = self.archive.get_stream_name(stream_id) {
                let raw_size = self.archive.get_raw_size(stream_id);
                let packed_size = self.archive.get_packed_size(stream_id);
                let num_parts = self.archive.get_num_parts(stream_id);
                stats.push((name.to_string(), raw_size, packed_size, num_parts));
            }
        }
        stats
    }

    /// Get list of contigs for a specific sample
    pub fn list_contigs(&mut self, sample_name: &str) -> Result<Vec<String>> {
        if self.config.verbosity > 1 {
            eprintln!("list_contigs called for: {sample_name}");
            eprintln!(
                "get_no_contigs before load: {:?}",
                self.collection.get_no_contigs(sample_name)
            );
        }

        // Load ALL contig batches if sample not found (samples may be in any batch)
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            if self.config.verbosity > 1 {
                eprintln!("Loading all contig batches for sample: {sample_name}");
            }
            let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
            for batch_id in 0..num_batches {
                self.collection
                    .load_contig_batch(&mut self.archive, batch_id)?;
            }

            if self.config.verbosity > 1 {
                eprintln!(
                    "After load: get_no_contigs = {:?}",
                    self.collection.get_no_contigs(sample_name)
                );
            }
        } else if self.config.verbosity > 1 {
            eprintln!(
                "Contig batch already loaded, get_no_contigs = {:?}",
                self.collection.get_no_contigs(sample_name)
            );
        }

        let result = self.collection.get_contig_list(sample_name);
        if self.config.verbosity > 1 {
            eprintln!("get_contig_list result: {result:?}");
        }

        result.ok_or_else(|| anyhow!("Sample not found: {sample_name}"))
    }

    /// Get the length of a contig without decompressing it
    ///
    /// This is O(1) once segment metadata is loaded, as it computes the length
    /// from the `raw_length` fields in segment descriptors.
    ///
    /// # Arguments
    /// * `sample_name` - The sample containing the contig
    /// * `contig_name` - The contig name
    ///
    /// # Returns
    /// The total length of the contig in bases
    pub fn get_contig_length(&mut self, sample_name: &str, contig_name: &str) -> Result<usize> {
        // Load contig batches if needed
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
            for batch_id in 0..num_batches {
                self.collection
                    .load_contig_batch(&mut self.archive, batch_id)?;
            }
        }

        // Get segment descriptors
        let segments = self
            .collection
            .get_contig_desc(sample_name, contig_name)
            .ok_or_else(|| anyhow!("Contig not found: {sample_name}/{contig_name}"))?;

        // Compute length: first segment fully counted, subsequent segments have k-mer overlap
        let kmer_len = self.kmer_length as usize;
        let mut total_length = 0usize;

        for (i, segment) in segments.iter().enumerate() {
            if i == 0 {
                total_length += segment.raw_length as usize;
            } else {
                // Subsequent segments overlap by kmer_length bases
                total_length += segment.raw_length as usize - kmer_len;
            }
        }

        Ok(total_length)
    }

    /// Extract a subsequence from a contig
    ///
    /// This is more efficient than `get_contig()` for small ranges, as it only
    /// decompresses the segments that overlap the requested range.
    ///
    /// # Arguments
    /// * `sample_name` - The sample containing the contig
    /// * `contig_name` - The contig name
    /// * `start` - 0-based start position (inclusive)
    /// * `end` - 0-based end position (exclusive)
    ///
    /// # Returns
    /// The subsequence as `Vec<u8>` in numeric encoding (0=A, 1=C, 2=G, 3=T)
    ///
    /// # Example
    /// ```no_run
    /// use ragc_core::{Decompressor, DecompressorConfig};
    ///
    /// let mut dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
    ///
    /// // Extract bases 1000-2000 from chr1
    /// let seq = dec.get_contig_range("sample", "chr1", 1000, 2000)?;
    /// assert_eq!(seq.len(), 1000);
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn get_contig_range(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Contig> {
        if start >= end {
            return Ok(Vec::new());
        }

        // Load contig batches if needed
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
            for batch_id in 0..num_batches {
                self.collection
                    .load_contig_batch(&mut self.archive, batch_id)?;
            }
        }

        // Get segment descriptors
        let segments = self
            .collection
            .get_contig_desc(sample_name, contig_name)
            .ok_or_else(|| anyhow!("Contig not found: {sample_name}/{contig_name}"))?;

        let kmer_len = self.kmer_length as usize;

        // Calculate segment positions in the contig
        // Each entry is (segment_start_in_contig, segment_end_in_contig, segment_idx)
        let mut segment_ranges: Vec<(usize, usize, usize)> = Vec::with_capacity(segments.len());
        let mut contig_pos = 0usize;

        for (i, segment) in segments.iter().enumerate() {
            let seg_len = segment.raw_length as usize;
            let contribution = if i == 0 { seg_len } else { seg_len - kmer_len };
            let seg_start = contig_pos;
            let seg_end = contig_pos + contribution;
            segment_ranges.push((seg_start, seg_end, i));
            contig_pos = seg_end;
        }

        let contig_len = contig_pos;

        // Clamp end to contig length
        let end = end.min(contig_len);
        if start >= end {
            return Ok(Vec::new());
        }

        // Find segments that overlap [start, end)
        let mut result = Vec::with_capacity(end - start);

        for (seg_start, seg_end, seg_idx) in segment_ranges {
            // Skip segments entirely before the range
            if seg_end <= start {
                continue;
            }
            // Stop if we've passed the range
            if seg_start >= end {
                break;
            }

            // This segment overlaps the range - decompress it
            let segment_desc = &segments[seg_idx];
            let mut segment_data = self.get_segment(segment_desc)?;

            // Apply reverse complement if needed
            if segment_desc.is_rev_comp {
                segment_data = Self::reverse_complement_segment(&segment_data);
            }

            // Determine which part of this segment's contribution we need
            // For segment 0: contribution starts at segment byte 0
            // For segment 1+: contribution starts at segment byte kmer_len
            let contribution_start_in_segment = if seg_idx == 0 { 0 } else { kmer_len };

            // Map [start, end) back to positions within segment_data
            // seg_start..seg_end maps to segment_data[contribution_start_in_segment..]
            //
            // If requested range starts at position P in contig, and segment contributes
            // from seg_start, then the offset within the contribution is P - seg_start.
            // The actual byte in segment_data is contribution_start_in_segment + (P - seg_start).

            let range_start_in_contribution = start.saturating_sub(seg_start);
            let range_end_in_contribution = (end - seg_start).min(seg_end - seg_start);

            let data_start = contribution_start_in_segment + range_start_in_contribution;
            let data_end = contribution_start_in_segment + range_end_in_contribution;

            if data_start < data_end && data_end <= segment_data.len() {
                result.extend_from_slice(&segment_data[data_start..data_end]);
            }
        }

        Ok(result)
    }

    /// Extract a specific contig from a sample
    pub fn get_contig(&mut self, sample_name: &str, contig_name: &str) -> Result<Contig> {
        // Load ALL contig batches if sample not found (samples may be in any batch)
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
            for batch_id in 0..num_batches {
                self.collection
                    .load_contig_batch(&mut self.archive, batch_id)?;
            }
        }

        // Get contig segment descriptors
        let segments = self
            .collection
            .get_contig_desc(sample_name, contig_name)
            .ok_or_else(|| anyhow!("Contig not found: {sample_name}/{contig_name}"))?;

        if self.config.verbosity > 1 {
            eprintln!(
                "Extracting {}/{} ({} segments)",
                sample_name,
                contig_name,
                segments.len()
            );
            for (i, seg) in segments.iter().enumerate() {
                eprintln!(
                    "  Segment[{}]: group_id={}, in_group_id={}, is_rev_comp={}, raw_length={}",
                    i, seg.group_id, seg.in_group_id, seg.is_rev_comp, seg.raw_length
                );
            }
        }

        // Reconstruct contig from segments
        self.reconstruct_contig(&segments)
    }

    /// Public: get segment descriptors for a contig
    pub fn get_contig_segments_desc(
        &mut self,
        sample_name: &str,
        contig_name: &str,
    ) -> Result<Vec<ragc_common::SegmentDesc>> {
        // Ensure batches loaded
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
            for batch_id in 0..num_batches {
                self.collection
                    .load_contig_batch(&mut self.archive, batch_id)?;
            }
        }
        self.collection
            .get_contig_desc(sample_name, contig_name)
            .ok_or_else(|| anyhow!("Contig not found: {}/{}", sample_name, contig_name))
    }

    /// Public: get raw segment data for a specific segment descriptor
    pub fn get_segment_data_by_desc(&mut self, desc: &ragc_common::SegmentDesc) -> Result<Contig> {
        self.get_segment(desc)
    }

    /// Public: get the reference segment for a group (loads and caches if needed)
    pub fn get_reference_segment(&mut self, group_id: u32) -> Result<Contig> {
        if let Some(ref_data) = self.segment_cache.get(&group_id) {
            return Ok(ref_data.clone());
        }

        let archive_version = ragc_common::AGC_FILE_MAJOR * 1000 + ragc_common::AGC_FILE_MINOR;
        let ref_stream_name = stream_ref_name(archive_version, group_id);
        let stream_id = self
            .archive
            .get_stream_id(&ref_stream_name)
            .ok_or_else(|| anyhow!("Reference stream not found: {}", ref_stream_name))?;

        let (mut data, metadata) = self.archive.get_part_by_id(stream_id, 0)?;
        // Decompress if needed; metadata holds original length for packed format
        let decompressed = if data.is_empty() {
            Vec::new()
        } else if data.last() == Some(&0) {
            // Plain ZSTD stream with marker 0
            data.pop();
            decompress_segment_with_marker(&data, 0)?
        } else {
            // Tuple-packed with marker 1
            let marker = data.pop().unwrap();
            decompress_segment_with_marker(&data, marker)?
        };

        // Unpack 2-bit encoded reference to 1-byte bases if needed
        // decompress_segment_with_marker returns bytes in the stored format for references
        // Our helper already returns decompressed raw bytes for references
        let reference = decompressed;
        self.segment_cache.insert(group_id, reference.clone());
        Ok(reference)
    }

    /// Extract all contigs from a sample
    pub fn get_sample(&mut self, sample_name: &str) -> Result<Vec<(String, Contig)>> {
        // Load ALL contig batches if needed (FIX: was only loading batch 0)
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            // Get number of contig batches from the archive
            let num_batches = self.collection.get_no_contig_batches(&self.archive)?;

            // Load ALL contig batches, not just batch 0
            for batch_id in 0..num_batches {
                self.collection
                    .load_contig_batch(&mut self.archive, batch_id)?;
            }
        }

        let sample_desc = self
            .collection
            .get_sample_desc(sample_name)
            .ok_or_else(|| anyhow!("Sample not found: {sample_name}"))?;

        let mut contigs = Vec::new();

        for (contig_name, segments) in sample_desc {
            if self.config.verbosity > 1 {
                eprintln!(
                    "Extracting {}/{} ({} segments)",
                    sample_name,
                    contig_name,
                    segments.len()
                );
            }

            let contig_data = self.reconstruct_contig(&segments)?;
            contigs.push((contig_name, contig_data));
        }

        Ok(contigs)
    }

    /// Unpack 2-bit encoded data (C++ AGC format) to 1-byte-per-base format
    /// C++ AGC stores 4 bases per byte: bits [7:6][5:4][3:2][1:0] = bases [0][1][2][3]
    /// Each 2-bit value: 00=A=0, 01=C=1, 10=G=2, 11=T=3
    fn unpack_2bit(packed: &[u8], expected_length: usize) -> Contig {
        let mut unpacked = Vec::with_capacity(expected_length);

        for &byte in packed {
            // Extract 4 bases from each byte (most significant bits first)
            unpacked.push((byte >> 6) & 0x03);
            unpacked.push((byte >> 4) & 0x03);
            unpacked.push((byte >> 2) & 0x03);
            unpacked.push(byte & 0x03);
        }

        // Truncate to expected length (last byte might not use all 4 slots)
        unpacked.truncate(expected_length);
        unpacked
    }

    /// Apply reverse complement to a segment (reverse order and complement each base)
    ///
    /// C++ AGC compatibility: Only complement bases 0-3 (A, C, G, T).
    /// IUPAC ambiguity codes (4+) are reversed but NOT complemented.
    /// This matches C++ AGC's behavior: `(*p < 4) ? 3 - *p : *p`
    fn reverse_complement_segment(segment: &[u8]) -> Contig {
        segment
            .iter()
            .rev()
            .map(|&base| {
                if base < 4 {
                    3 - base // Complement: A(0)<->T(3), C(1)<->G(2)
                } else {
                    base // Keep IUPAC codes unchanged (just reverse position)
                }
            })
            .collect()
    }

    /// Reconstruct a contig from its segment descriptors
    fn reconstruct_contig(&mut self, segments: &[SegmentDesc]) -> Result<Contig> {
        let mut contig = Contig::new();
        #[cfg(feature = "verbose_debug")]
        let should_debug = crate::env_cache::debug_reconstruct() && segments.len() > 1;
        #[cfg(not(feature = "verbose_debug"))]
        let should_debug = false;

        if should_debug {
            eprintln!(
                "DEBUG_RECONSTRUCT: {} segments, k={}",
                segments.len(),
                self.kmer_length
            );
        }

        for (i, segment_desc) in segments.iter().enumerate() {
            let mut segment_data = self.get_segment(segment_desc)?;

            // Apply reverse complement if the segment was stored in reverse orientation
            if segment_desc.is_rev_comp {
                segment_data = Self::reverse_complement_segment(&segment_data);
            }

            if i == 0 {
                // First segment: add completely
                if should_debug {
                    let last_k: Vec<u8> = segment_data
                        .iter()
                        .rev()
                        .take(self.kmer_length as usize)
                        .rev()
                        .copied()
                        .collect();
                    eprintln!(
                        "  Seg {}: len={} (added fully), last_{}_bytes={:?}, contig_len={}",
                        i,
                        segment_data.len(),
                        self.kmer_length,
                        last_k,
                        contig.len() + segment_data.len()
                    );
                }
                contig.extend_from_slice(&segment_data);
            } else {
                // Subsequent segments: skip first k bases (full k-mer overlap)
                // Each segment includes the FULL k-mer at the start
                let overlap = self.kmer_length as usize;
                if segment_data.len() < overlap {
                    eprintln!(
                        "ERROR: Segment {} too short! Length={}, overlap={}",
                        i,
                        segment_data.len(),
                        overlap
                    );
                    eprintln!(
                        "  Segment desc: group_id={}, in_group_id={}, raw_length={}",
                        segment_desc.group_id, segment_desc.in_group_id, segment_desc.raw_length
                    );
                    anyhow::bail!("Corrupted archive: segment too short (segment {}, got {} bytes, need at least {} bytes)", i, segment_data.len(), overlap);
                }
                let added = segment_data.len() - overlap;
                let first_k: Vec<u8> = segment_data.iter().take(overlap).copied().collect();
                if should_debug {
                    eprintln!(
                        "  Seg {}: len={} (skip {}), first_{}_bytes={:?}, added={}, contig_len={}",
                        i,
                        segment_data.len(),
                        overlap,
                        overlap,
                        first_k,
                        added,
                        contig.len() + added
                    );
                }
                contig.extend_from_slice(&segment_data[overlap..]);
            }
        }

        if should_debug {
            eprintln!("  Final contig length: {}", contig.len());
        }

        Ok(contig)
    }

    /// Unpack contigs from packed data by splitting on CONTIG_SEPARATOR (0xFF)
    fn unpack_contig(packed_data: &[u8], position_in_pack: usize) -> Result<Contig> {
        let mut contig_start = 0;
        let mut current_position = 0;

        for (i, &byte) in packed_data.iter().enumerate() {
            if byte == CONTIG_SEPARATOR {
                if current_position == position_in_pack {
                    // Found the contig we're looking for
                    return Ok(packed_data[contig_start..i].to_vec());
                }
                current_position += 1;
                contig_start = i + 1;
            }
        }

        // Handle last contig (if no trailing separator or position is last)
        if current_position == position_in_pack {
            return Ok(packed_data[contig_start..].to_vec());
        }

        anyhow::bail!(
            "Position {} not found in packed data (only {} contigs found)",
            position_in_pack,
            current_position + 1
        )
    }

    /// Get a single segment (handles reference and LZ diff decoding)
    /// Supports packed-contig mode where multiple contigs are stored in one part
    ///
    /// Two-stream architecture:
    /// - Raw groups (0-15): All segments in delta stream
    /// - LZ groups (16+): Reference in ref stream (part 0), LZ-encoded segments in delta stream
    fn get_segment(&mut self, desc: &SegmentDesc) -> Result<Contig> {
        const PACK_CARDINALITY: usize = 50; // C++ default
        const NO_RAW_GROUPS: u32 = 16;

        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;

        if self.config.verbosity > 1 {
            eprintln!(
                "get_segment: group_id={}, in_group_id={}, raw_length={}",
                desc.group_id, desc.in_group_id, desc.raw_length
            );
        }

        // Handle LZ groups (>= 16) differently from raw groups (< 16)
        if desc.group_id >= NO_RAW_GROUPS {
            // LZ group: two-stream architecture

            // First, ensure we have the reference loaded
            if !self.segment_cache.contains_key(&desc.group_id) {
                // Load reference from ref stream (part 0)
                let ref_stream_name = stream_ref_name(archive_version, desc.group_id);
                let ref_stream_id = self
                    .archive
                    .get_stream_id(&ref_stream_name)
                    .ok_or_else(|| anyhow!("Reference stream not found: {ref_stream_name}"))?;

                let (mut ref_data, ref_metadata) = self.archive.get_part_by_id(ref_stream_id, 0)?;

                // Decompress if needed
                let mut decompressed_ref = if ref_metadata == 0 {
                    ref_data
                } else {
                    if ref_data.is_empty() {
                        anyhow::bail!("Empty compressed reference data");
                    }
                    let marker = ref_data.pop().unwrap();
                    decompress_segment_with_marker(&ref_data, marker)?
                };

                // Check if data is in 2-bit packed format (C++ AGC compatibility)
                // If decompressed length is ~1/4 of raw_length, it's packed
                // Allow small tolerance for rounding (last byte can encode 1-4 bases)
                //
                // CRITICAL: Use ref_metadata (original reference size) NOT desc.raw_length (requested segment size)
                // When ref_metadata != 0, it contains the original uncompressed size of the REFERENCE.
                // desc.raw_length is the size of the segment being requested, which may differ from the reference.
                let expected_ref_len = if ref_metadata != 0 {
                    ref_metadata as usize
                } else {
                    // If ref_metadata == 0, the reference was stored uncompressed, so use actual size
                    decompressed_ref.len()
                };

                let is_packed = decompressed_ref.len() * 4 >= expected_ref_len
                    && decompressed_ref.len() * 4 < expected_ref_len + 8;

                if self.config.verbosity > 2 {
                    eprintln!(
                        "  DEBUG: decompressed_len={}, expected_ref_len={}, decompressed*4={}, is_packed={}",
                        decompressed_ref.len(),
                        expected_ref_len,
                        decompressed_ref.len() * 4,
                        is_packed
                    );
                }

                if is_packed {
                    // Unpack from 2-bit format to 1-byte-per-base format
                    decompressed_ref = Self::unpack_2bit(&decompressed_ref, expected_ref_len);

                    if self.config.verbosity > 1 {
                        eprintln!(
                            "Loaded & unpacked reference for group {}: length={} (was {} bytes packed)",
                            desc.group_id,
                            decompressed_ref.len(),
                            decompressed_ref.len() / 4
                        );
                    }
                } else if self.config.verbosity > 1 {
                    eprintln!(
                        "Loaded reference for group {}: length={} (already unpacked)",
                        desc.group_id,
                        decompressed_ref.len()
                    );
                }

                // Cache the reference
                self.segment_cache.insert(desc.group_id, decompressed_ref);
            }

            // If this IS the reference (in_group_id == 0), return it directly
            if desc.in_group_id == 0 {
                return Ok(self.segment_cache.get(&desc.group_id).unwrap().clone());
            }

            // Otherwise, load LZ-encoded segment from delta stream
            // For delta stream: in_group_id 1 is at position 0, in_group_id 2 is at position 1, etc.
            let delta_position = (desc.in_group_id - 1) as usize;
            let pack_id = delta_position / PACK_CARDINALITY;
            let position_in_pack = delta_position % PACK_CARDINALITY;

            if self.config.verbosity > 1 {
                eprintln!(
                    "  LZ group: delta_position={delta_position}, pack_id={pack_id}, position_in_pack={position_in_pack}"
                );
            }

            let delta_stream_name = stream_delta_name(archive_version, desc.group_id);
            let delta_stream_id =
                self.archive
                    .get_stream_id(&delta_stream_name)
                    .ok_or_else(|| {
                        eprintln!("ERROR: Delta stream not found: {}", delta_stream_name);
                        eprintln!(
                            "  Requested by segment: group_id={}, in_group_id={}, raw_length={}",
                            desc.group_id, desc.in_group_id, desc.raw_length
                        );
                        anyhow!("Delta stream not found: {delta_stream_name}")
                    })?;

            // Check how many parts this stream has
            let num_parts = self.archive.get_num_parts(delta_stream_id);

            if pack_id >= num_parts {
                anyhow::bail!("Pack ID {} out of range (stream has only {} parts). in_group_id={}, delta_position={}",
                    pack_id, num_parts, desc.in_group_id, delta_position);
            }

            let (mut delta_data, delta_metadata) =
                self.archive.get_part_by_id(delta_stream_id, pack_id)?;

            // Decompress if needed
            let decompressed_pack = if delta_metadata == 0 {
                delta_data
            } else {
                if delta_data.is_empty() {
                    anyhow::bail!("Empty compressed delta data");
                }
                let marker = delta_data.pop().unwrap();
                decompress_segment_with_marker(&delta_data, marker)?
            };

            // Unpack the specific LZ-encoded segment
            let lz_encoded = Self::unpack_contig(&decompressed_pack, position_in_pack)?;

            if self.config.verbosity > 1 {
                eprintln!("Unpacked LZ-encoded segment: length={}", lz_encoded.len());
            }

            // Decode using reference
            let reference = self
                .segment_cache
                .get(&desc.group_id)
                .ok_or_else(|| anyhow!("Reference not loaded for group {}", desc.group_id))?;

            let mut lz_diff = LZDiff::new(self.min_match_len);
            lz_diff.prepare(reference);

            // Handle empty encoding (segment equals reference)
            let decoded = if lz_encoded.is_empty() {
                reference.clone()
            } else {
                lz_diff.decode(&lz_encoded)
            };

            // Debug: check for corruption pattern (AAAA = zeros)
            if crate::env_cache::debug_decode() {
                // Critical: check if decoded length mismatches raw_length
                if decoded.len() as u32 != desc.raw_length {
                    eprintln!("RAGC_DEBUG_DECODE_MISMATCH: group={} in_group={} raw_len={} decoded_len={} ref_len={} lz_encoded_len={}",
                        desc.group_id, desc.in_group_id, desc.raw_length, decoded.len(), reference.len(), lz_encoded.len());
                    eprintln!(
                        "  LZ encoded bytes (first 100): {:?}",
                        &lz_encoded[..lz_encoded.len().min(100)]
                    );

                    // Show first few bytes of decoded and expected
                    let dec_start: Vec<u8> = decoded.iter().take(30).copied().collect();
                    let ref_start: Vec<u8> = reference.iter().take(30).copied().collect();
                    eprintln!("  Decoded[0:30]={:?}", dec_start);
                    eprintln!("  Reference[0:30]={:?}", ref_start);
                }
            }

            if self.config.verbosity > 1 {
                eprintln!("Decoded segment: length={}", decoded.len());
            }

            Ok(decoded)
        } else {
            // Raw group (0-15): NO reference stream, ALL data in delta stream
            // This matches C++ AGC's segment.cpp get_raw() function:
            // - part_id = id_seq / contigs_in_pack
            // - seq_in_part_id = id_seq % contigs_in_pack
            // Data is stored raw (not LZ-encoded), separated by contig_separator (0x7f)

            let pack_id = desc.in_group_id as usize / PACK_CARDINALITY;
            let position_in_pack = desc.in_group_id as usize % PACK_CARDINALITY;

            if self.config.verbosity > 1 {
                eprintln!(
                    "  Raw group {}: in_group_id={}, pack_id={}, position_in_pack={}",
                    desc.group_id, desc.in_group_id, pack_id, position_in_pack
                );
            }

            let stream_name = stream_delta_name(archive_version, desc.group_id);
            let stream_id = self
                .archive
                .get_stream_id(&stream_name)
                .ok_or_else(|| anyhow!("Delta stream not found: {stream_name}"))?;

            let (mut data, metadata) = self.archive.get_part_by_id(stream_id, pack_id)?;

            // Decompress if needed
            let decompressed_pack = if metadata == 0 {
                data
            } else {
                if data.is_empty() {
                    anyhow::bail!("Empty compressed data for raw group");
                }
                let marker = data.pop().unwrap();
                decompress_segment_with_marker(&data, marker)?
            };

            // Unpack the raw segment (segments separated by 0x7f)
            let contig_data = Self::unpack_contig(&decompressed_pack, position_in_pack)?;

            if self.config.verbosity > 1 {
                eprintln!("Unpacked raw segment: length={}", contig_data.len());
            }

            Ok(contig_data)
        }
    }

    /// Get all sample names that match a given prefix
    ///
    /// This is useful for extracting all samples from a specific genome or haplotype.
    ///
    /// # Example
    /// ```no_run
    /// use ragc_core::{Decompressor, DecompressorConfig};
    ///
    /// let dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
    ///
    /// // Get all samples starting with "AAA"
    /// let aaa_samples = dec.list_samples_with_prefix("AAA");
    /// // Returns: ["AAA#0", "AAA#1", ...] if they exist
    ///
    /// // Get all samples with haplotype 0
    /// let hap0_samples = dec.list_samples_with_prefix("AAA#0");
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn list_samples_with_prefix(&self, prefix: &str) -> Vec<String> {
        self.list_samples()
            .into_iter()
            .filter(|s| s.starts_with(prefix))
            .collect()
    }

    /// Extract multiple samples matching a prefix
    ///
    /// Returns a HashMap mapping sample names to their contigs.
    /// Each contig is represented as (contig_name, sequence).
    ///
    /// # Example
    /// ```no_run
    /// use ragc_core::{Decompressor, DecompressorConfig};
    ///
    /// let mut dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
    ///
    /// // Extract all AAA samples
    /// let aaa_samples = dec.get_samples_by_prefix("AAA")?;
    ///
    /// for (sample_name, contigs) in aaa_samples {
    ///     println!("Sample {}: {} contigs", sample_name, contigs.len());
    ///     for (contig_name, sequence) in contigs {
    ///         println!("  {}: {} bp", contig_name, sequence.len());
    ///     }
    /// }
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn get_samples_by_prefix(
        &mut self,
        prefix: &str,
    ) -> Result<HashMap<String, Vec<(String, Contig)>>> {
        let sample_names = self.list_samples_with_prefix(prefix);
        let mut results = HashMap::new();

        for sample_name in sample_names {
            let contigs = self.get_sample(&sample_name)?;
            results.insert(sample_name, contigs);
        }

        Ok(results)
    }

    /// Clone this decompressor for use in another thread
    ///
    /// This creates an independent reader that can be used from a different thread.
    /// The clone is lightweight as it shares the archive data structure.
    ///
    /// # Thread Safety
    /// Each clone is fully independent and can be used from a separate thread.
    /// The segment cache is NOT shared between clones.
    ///
    /// # Example
    /// ```no_run
    /// use ragc_core::{Decompressor, DecompressorConfig};
    /// use std::thread;
    ///
    /// let dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
    /// let samples = dec.list_samples();
    ///
    /// // Spawn threads to extract samples in parallel
    /// let handles: Vec<_> = samples.into_iter().map(|sample_name| {
    ///     let mut thread_dec = dec.clone_for_thread().unwrap();
    ///     thread::spawn(move || {
    ///         thread_dec.get_sample(&sample_name)
    ///     })
    /// }).collect();
    ///
    /// // Collect results
    /// for handle in handles {
    ///     let result = handle.join().unwrap()?;
    ///     // Process result...
    /// }
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn clone_for_thread(&self) -> Result<Self> {
        // Re-open the archive (cheap - just re-reads metadata)
        // Archive data is memory-mapped or cached by the OS
        Self::open(&self.archive_path, self.config.clone())
    }

    /// Write a sample to a FASTA file
    pub fn write_sample_fasta(&mut self, sample_name: &str, output_path: &Path) -> Result<()> {
        let contigs = self.get_sample(sample_name)?;

        let mut writer = GenomeWriter::<File>::create(output_path)?;

        for (contig_name, contig_data) in contigs {
            // Convert numeric encoding back to ASCII using CNV_NUM lookup table
            // CNV_NUM[0..16] = [A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V, U]
            // This preserves all IUPAC ambiguity codes, not just ACGT
            let mut ascii_data = Vec::with_capacity(contig_data.len());
            for &base in &contig_data {
                // Use CNV_NUM[0..16] for proper IUPAC code decoding
                // Values 0-15 map to valid bases, anything else (shouldn't happen) maps to N
                let ascii_base = if base < 16 {
                    CNV_NUM[base as usize]
                } else {
                    b'N'
                };
                ascii_data.push(ascii_base);
            }

            writer.save_contig_directly(&contig_name, &ascii_data, 0)?;
        }

        Ok(())
    }

    /// Get archive inspection data for debugging/comparison
    /// Returns: (group_id, num_segments, num_reference_segments, num_delta_segments)
    pub fn get_group_statistics(&mut self) -> Result<Vec<(u32, usize, usize, usize)>> {
        use std::collections::HashMap;

        // Load ALL contig batches to get complete segment information
        let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
        for batch_id in 0..num_batches {
            self.collection
                .load_contig_batch(&mut self.archive, batch_id)?;
        }

        let sample_names = self.list_samples();

        // Now collect statistics about groups
        let mut group_stats: HashMap<u32, (usize, usize, usize)> = HashMap::new();

        for sample_name in &sample_names {
            let contig_list = self
                .collection
                .get_contig_list(sample_name)
                .ok_or_else(|| anyhow!("Failed to get contig list for {}", sample_name))?;

            for contig_name in &contig_list {
                if let Some(segments) = self.collection.get_contig_desc(sample_name, contig_name) {
                    for seg in segments {
                        let entry = group_stats.entry(seg.group_id).or_insert((0, 0, 0));
                        entry.0 += 1; // total segments
                                      // Raw groups (0-15) have NO references - all segments are raw/delta
                                      // LZ groups (16+) have references at in_group_id==0
                        if seg.group_id >= 16 && seg.in_group_id == 0 {
                            entry.1 += 1; // reference segments (LZ groups only)
                        } else {
                            entry.2 += 1; // delta/raw segments
                        }
                    }
                }
            }
        }

        let mut stats: Vec<_> = group_stats
            .into_iter()
            .map(|(gid, (total, refs, deltas))| (gid, total, refs, deltas))
            .collect();
        stats.sort_by_key(|s| s.0);

        Ok(stats)
    }

    /// Get detailed segment information for all contigs
    /// Returns: Vec<(sample_name, contig_name, Vec<SegmentDesc>)>
    pub fn get_all_segments(&mut self) -> Result<Vec<(String, String, Vec<SegmentDesc>)>> {
        // Load ALL contig batches
        let num_batches = self.collection.get_no_contig_batches(&self.archive)?;
        for batch_id in 0..num_batches {
            self.collection
                .load_contig_batch(&mut self.archive, batch_id)?;
        }

        let sample_names = self.list_samples();
        let mut all_segments = Vec::new();

        for sample_name in &sample_names {
            let contig_list = self
                .collection
                .get_contig_list(sample_name)
                .ok_or_else(|| anyhow!("Failed to get contig list for {}", sample_name))?;

            for contig_name in &contig_list {
                if let Some(segments) = self.collection.get_contig_desc(sample_name, contig_name) {
                    all_segments.push((
                        sample_name.clone(),
                        contig_name.clone(),
                        segments.to_vec(),
                    ));
                }
            }
        }

        Ok(all_segments)
    }

    /// Close the archive
    pub fn close(mut self) -> Result<()> {
        self.archive.close()
    }
}
