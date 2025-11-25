// AGC Decompressor
// Extracts genomes from AGC archives

use crate::{
    genome_io::GenomeWriter, kmer::reverse_complement, lz_diff::LZDiff,
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

    // Track which contig batches have been loaded
    loaded_contig_batches: std::collections::HashSet<usize>,

    // Archive parameters
    _segment_size: u32,
    kmer_length: u32,
    min_match_len: u32,

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
            loaded_contig_batches: std::collections::HashSet::new(),
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

    /// Load all contig batches if not already loaded
    fn ensure_contig_batches_loaded(&mut self) -> Result<()> {
        let num_batches = self.collection.get_contig_stream_num_batches(&self.archive)?;

        for batch_id in 0..num_batches {
            if !self.loaded_contig_batches.contains(&batch_id) {
                self.collection.load_contig_batch(&mut self.archive, batch_id)?;
                self.loaded_contig_batches.insert(batch_id);
            }
        }
        Ok(())
    }

    /// Get list of contigs for a specific sample
    pub fn list_contigs(&mut self, sample_name: &str) -> Result<Vec<String>> {
        self.ensure_contig_batches_loaded()?;
        self.collection.get_contig_list(sample_name)
            .ok_or_else(|| anyhow!("Sample not found: {sample_name}"))
    }

    /// Extract a specific contig from a sample
    pub fn get_contig(&mut self, sample_name: &str, contig_name: &str) -> Result<Contig> {
        self.ensure_contig_batches_loaded()?;

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

    /// Extract all contigs from a sample
    pub fn get_sample(&mut self, sample_name: &str) -> Result<Vec<(String, Contig)>> {
        self.ensure_contig_batches_loaded()?;

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
    fn reverse_complement_segment(segment: &[u8]) -> Contig {
        segment
            .iter()
            .rev()
            .map(|&base| reverse_complement(base as u64) as u8)
            .collect()
    }

    /// Reconstruct a contig from its segment descriptors
    fn reconstruct_contig(&mut self, segments: &[SegmentDesc]) -> Result<Contig> {
        let mut contig = Contig::new();

        for (i, segment_desc) in segments.iter().enumerate() {
            let mut segment_data = self.get_segment(segment_desc)?;

            // Apply reverse complement if needed
            if segment_desc.is_rev_comp {
                segment_data = Self::reverse_complement_segment(&segment_data);
            }

            if i == 0 {
                // First segment: add completely
                contig.extend_from_slice(&segment_data);
            } else {
                // Subsequent segments: skip first kmer_length bases (overlap with previous segment)
                if segment_data.len() < self.kmer_length as usize {
                    eprintln!(
                        "ERROR: Segment {} too short! Length={}, kmer_length={}",
                        i,
                        segment_data.len(),
                        self.kmer_length
                    );
                    eprintln!(
                        "  Segment desc: group_id={}, in_group_id={}, raw_length={}",
                        segment_desc.group_id, segment_desc.in_group_id, segment_desc.raw_length
                    );
                    anyhow::bail!("Corrupted archive: segment too short (segment {}, got {} bytes, need at least {} bytes)", i, segment_data.len(), self.kmer_length);
                }
                contig.extend_from_slice(&segment_data[self.kmer_length as usize..]);
            }
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
                let is_packed = decompressed_ref.len() * 4 >= desc.raw_length as usize
                    && decompressed_ref.len() * 4 < desc.raw_length as usize + 8;

                if self.config.verbosity > 2 {
                    eprintln!(
                        "  DEBUG: decompressed_len={}, raw_length={}, decompressed*4={}, is_packed={}",
                        decompressed_ref.len(),
                        desc.raw_length,
                        decompressed_ref.len() * 4,
                        is_packed
                    );
                }

                if is_packed {
                    // Unpack from 2-bit format to 1-byte-per-base format
                    decompressed_ref =
                        Self::unpack_2bit(&decompressed_ref, desc.raw_length as usize);

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

            if self.config.verbosity > 1 {
                eprintln!("Decoded segment: length={}", decoded.len());
            }

            Ok(decoded)
        } else {
            // Raw-only group (0-15): all segments in delta stream
            let pack_id = desc.in_group_id as usize / PACK_CARDINALITY;
            let position_in_pack = desc.in_group_id as usize % PACK_CARDINALITY;

            if self.config.verbosity > 1 {
                eprintln!("  Raw group: pack_id={pack_id}, position_in_pack={position_in_pack}");
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
                    anyhow::bail!("Empty compressed data");
                }
                let marker = data.pop().unwrap();
                decompress_segment_with_marker(&data, marker)?
            };

            // Unpack the raw segment
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
            // Convert numeric encoding back to ASCII (0→A, 1→C, 2→G, 3→T)
            let mut ascii_data = Vec::with_capacity(contig_data.len());
            for &base in &contig_data {
                let ascii_base = match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => b'N',
                };
                ascii_data.push(ascii_base);
            }

            writer.save_contig_directly(&contig_name, &ascii_data, 0)?;
        }

        Ok(())
    }

    /// Compute the length of a contig from its segment descriptors without decompression.
    ///
    /// The total length is: sum of (raw_length - kmer_length) for each segment after the first,
    /// plus the full raw_length of the first segment.
    fn compute_contig_length_from_segments(&self, segments: &[SegmentDesc]) -> usize {
        if segments.is_empty() {
            return 0;
        }

        let mut total_len = segments[0].raw_length as usize;
        for seg in &segments[1..] {
            // Each subsequent segment contributes (raw_length - kmer_length) bases
            // because the first kmer_length bases overlap with the previous segment
            total_len += seg.raw_length as usize - self.kmer_length as usize;
        }
        total_len
    }

    /// Compute segment boundaries (start positions) for a contig.
    /// Returns a vector of (segment_start, segment_end) positions in contig coordinates.
    fn compute_segment_boundaries(&self, segments: &[SegmentDesc]) -> Vec<(usize, usize)> {
        let mut boundaries = Vec::with_capacity(segments.len());
        let mut pos = 0usize;

        for (i, seg) in segments.iter().enumerate() {
            let seg_len = if i == 0 {
                seg.raw_length as usize
            } else {
                seg.raw_length as usize - self.kmer_length as usize
            };
            boundaries.push((pos, pos + seg_len));
            pos += seg_len;
        }

        boundaries
    }

    /// Extract a subsequence from a specific contig.
    ///
    /// This method only decompresses the segments that overlap with the requested range,
    /// which can be significantly faster than extracting the full contig for small queries.
    ///
    /// # Arguments
    /// * `sample_name` - The sample name
    /// * `contig_name` - The contig name
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Returns
    /// The subsequence as a Vec<u8> with numeric encoding (A=0, C=1, G=2, T=3)
    ///
    /// # Example
    /// ```no_run
    /// use ragc_core::{Decompressor, DecompressorConfig};
    ///
    /// let mut dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
    /// // Extract bases 100-200 from chr1
    /// let subseq = dec.get_contig_range("sample1", "chr1", 100, 200)?;
    /// assert_eq!(subseq.len(), 100);
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn get_contig_range(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Contig> {
        if start > end {
            anyhow::bail!(
                "Start position {} exceeds end position {}",
                start,
                end
            );
        }

        self.ensure_contig_batches_loaded()?;

        // Get contig segment descriptors
        let segments = self
            .collection
            .get_contig_desc(sample_name, contig_name)
            .ok_or_else(|| anyhow!("Contig not found: {sample_name}/{contig_name}"))?;

        // Compute total contig length and validate range
        let contig_len = self.compute_contig_length_from_segments(&segments);
        if end > contig_len {
            anyhow::bail!(
                "End position {} exceeds contig length {}",
                end,
                contig_len
            );
        }

        // For small contigs or when requesting most of the contig, just decompress everything
        // The overhead of partial decompression isn't worth it for small data
        if contig_len < 100_000 || (end - start) > contig_len / 2 {
            let contig = self.get_contig(sample_name, contig_name)?;
            return Ok(contig[start..end].to_vec());
        }

        // Compute segment boundaries
        let boundaries = self.compute_segment_boundaries(&segments);

        // Find which segments overlap with [start, end)
        let mut first_seg_idx = None;
        let mut last_seg_idx = None;

        for (i, &(seg_start, seg_end)) in boundaries.iter().enumerate() {
            if seg_end > start && first_seg_idx.is_none() {
                first_seg_idx = Some(i);
            }
            if seg_start < end {
                last_seg_idx = Some(i);
            }
        }

        let first_seg_idx = first_seg_idx.ok_or_else(|| anyhow!("No segments found for range"))?;
        let last_seg_idx = last_seg_idx.ok_or_else(|| anyhow!("No segments found for range"))?;

        // Decompress only the overlapping segments
        let mut result = Contig::new();

        for seg_idx in first_seg_idx..=last_seg_idx {
            let segment_desc = &segments[seg_idx];
            let mut segment_data = self.get_segment(segment_desc)?;

            // Apply reverse complement if needed
            if segment_desc.is_rev_comp {
                segment_data = Self::reverse_complement_segment(&segment_data);
            }

            let (seg_start, _seg_end) = boundaries[seg_idx];

            if seg_idx == first_seg_idx && seg_idx == last_seg_idx {
                // Single segment case
                let local_start = start - seg_start;
                let local_end = end - seg_start;
                if seg_idx == 0 {
                    result.extend_from_slice(&segment_data[local_start..local_end]);
                } else {
                    // Skip kmer_length overlap
                    let adjusted_start = local_start + self.kmer_length as usize;
                    let adjusted_end = local_end + self.kmer_length as usize;
                    result.extend_from_slice(&segment_data[adjusted_start..adjusted_end]);
                }
            } else if seg_idx == first_seg_idx {
                // First segment of range
                let local_start = start - seg_start;
                if seg_idx == 0 {
                    result.extend_from_slice(&segment_data[local_start..]);
                } else {
                    let adjusted_start = local_start + self.kmer_length as usize;
                    result.extend_from_slice(&segment_data[adjusted_start..]);
                }
            } else if seg_idx == last_seg_idx {
                // Last segment of range
                let local_end = end - seg_start + self.kmer_length as usize;
                result.extend_from_slice(&segment_data[self.kmer_length as usize..local_end]);
            } else {
                // Middle segment - take everything except the kmer overlap
                result.extend_from_slice(&segment_data[self.kmer_length as usize..]);
            }
        }

        Ok(result)
    }

    /// Get the length of a contig without full decompression.
    ///
    /// This method computes the length from segment metadata, which is much faster
    /// than decompressing the entire contig.
    ///
    /// # Example
    /// ```no_run
    /// use ragc_core::{Decompressor, DecompressorConfig};
    ///
    /// let mut dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
    /// let len = dec.get_contig_length("sample1", "chr1")?;
    /// println!("Contig length: {} bp", len);
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    pub fn get_contig_length(&mut self, sample_name: &str, contig_name: &str) -> Result<usize> {
        self.ensure_contig_batches_loaded()?;

        // Get contig segment descriptors
        let segments = self
            .collection
            .get_contig_desc(sample_name, contig_name)
            .ok_or_else(|| anyhow!("Contig not found: {sample_name}/{contig_name}"))?;

        Ok(self.compute_contig_length_from_segments(&segments))
    }

    /// Close the archive
    pub fn close(mut self) -> Result<()> {
        self.archive.close()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Compressor, CompressorConfig};
    use std::fs;

    #[test]
    fn test_decompressor_basic() {
        let archive_path = "/tmp/test_decompress.agc";
        let _ = fs::remove_file(archive_path);

        // Create archive
        {
            let config = CompressorConfig::default();
            let mut compressor = Compressor::new(archive_path, config).unwrap();

            let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
            compressor
                .add_contig("sample1", "chr1", seq1.clone())
                .unwrap();

            compressor.finalize().unwrap();
        }

        // Decompress
        {
            let config = DecompressorConfig { verbosity: 2 }; // Enable debug output
            let mut decompressor = Decompressor::open(archive_path, config).unwrap();

            eprintln!("=== Test: list_samples ===");
            // List samples
            let samples = decompressor.list_samples();
            eprintln!("Samples: {samples:?}");
            assert_eq!(samples.len(), 1);
            assert_eq!(samples[0], "sample1");

            eprintln!("=== Test: list_contigs ===");
            // List contigs
            let contigs = decompressor.list_contigs("sample1").unwrap();
            eprintln!("Contigs: {contigs:?}");
            assert_eq!(contigs.len(), 1);
            assert_eq!(contigs[0], "chr1");

            // Extract contig
            let extracted = decompressor.get_contig("sample1", "chr1").unwrap();
            assert_eq!(extracted, vec![0, 1, 2, 3, 0, 1, 2, 3]);

            decompressor.close().unwrap();
        }

        fs::remove_file(archive_path).unwrap();
    }

    #[test]
    fn test_decompressor_multiple_contigs() {
        let archive_path = "/tmp/test_decompress_multi.agc";
        let _ = fs::remove_file(archive_path);

        // Create archive
        {
            let config = CompressorConfig::default();
            let mut compressor = Compressor::new(archive_path, config).unwrap();

            let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3];
            let seq2 = vec![3, 2, 1, 0, 3, 2, 1, 0];

            compressor
                .add_contig("sample1", "chr1", seq1.clone())
                .unwrap();
            compressor
                .add_contig("sample1", "chr2", seq2.clone())
                .unwrap();

            compressor.finalize().unwrap();
        }

        // Decompress
        {
            let config = DecompressorConfig { verbosity: 0 };
            let mut decompressor = Decompressor::open(archive_path, config).unwrap();

            // Extract full sample
            let sample_contigs = decompressor.get_sample("sample1").unwrap();
            assert_eq!(sample_contigs.len(), 2);

            // Verify data
            assert_eq!(sample_contigs[0].0, "chr1");
            assert_eq!(sample_contigs[0].1, vec![0, 1, 2, 3, 0, 1, 2, 3]);

            assert_eq!(sample_contigs[1].0, "chr2");
            assert_eq!(sample_contigs[1].1, vec![3, 2, 1, 0, 3, 2, 1, 0]);

            decompressor.close().unwrap();
        }

        fs::remove_file(archive_path).unwrap();
    }

    #[test]
    fn test_get_contig_range() {
        let archive_path = "/tmp/test_range.agc";
        let _ = fs::remove_file(archive_path);

        // Create archive with a known sequence
        {
            let config = CompressorConfig::default();
            let mut compressor = Compressor::new(archive_path, config).unwrap();

            // Create a sequence: ACGTACGTACGT (12 bases)
            let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
            compressor
                .add_contig("sample1", "chr1", seq1.clone())
                .unwrap();

            compressor.finalize().unwrap();
        }

        // Test range extraction
        {
            let config = DecompressorConfig { verbosity: 0 };
            let mut decompressor = Decompressor::open(archive_path, config).unwrap();

            // Test various ranges
            let range1 = decompressor.get_contig_range("sample1", "chr1", 0, 4).unwrap();
            assert_eq!(range1, vec![0, 1, 2, 3]); // ACGT

            let range2 = decompressor.get_contig_range("sample1", "chr1", 4, 8).unwrap();
            assert_eq!(range2, vec![0, 1, 2, 3]); // ACGT

            let range3 = decompressor.get_contig_range("sample1", "chr1", 2, 10).unwrap();
            assert_eq!(range3, vec![2, 3, 0, 1, 2, 3, 0, 1]); // GTACGTAC

            // Test full range
            let full = decompressor.get_contig_range("sample1", "chr1", 0, 12).unwrap();
            assert_eq!(full, vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]);

            // Test single base
            let single = decompressor.get_contig_range("sample1", "chr1", 5, 6).unwrap();
            assert_eq!(single, vec![1]); // C

            decompressor.close().unwrap();
        }

        fs::remove_file(archive_path).unwrap();
    }

    #[test]
    fn test_get_contig_length() {
        let archive_path = "/tmp/test_length.agc";
        let _ = fs::remove_file(archive_path);

        // Create archive
        {
            let config = CompressorConfig::default();
            let mut compressor = Compressor::new(archive_path, config).unwrap();

            let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3]; // 8 bases
            let seq2 = vec![0; 100]; // 100 bases

            compressor
                .add_contig("sample1", "chr1", seq1.clone())
                .unwrap();
            compressor
                .add_contig("sample1", "chr2", seq2.clone())
                .unwrap();

            compressor.finalize().unwrap();
        }

        // Test length computation
        {
            let config = DecompressorConfig { verbosity: 0 };
            let mut decompressor = Decompressor::open(archive_path, config).unwrap();

            let len1 = decompressor.get_contig_length("sample1", "chr1").unwrap();
            assert_eq!(len1, 8);

            let len2 = decompressor.get_contig_length("sample1", "chr2").unwrap();
            assert_eq!(len2, 100);

            decompressor.close().unwrap();
        }

        fs::remove_file(archive_path).unwrap();
    }

    #[test]
    fn test_get_contig_range_errors() {
        let archive_path = "/tmp/test_range_errors.agc";
        let _ = fs::remove_file(archive_path);

        // Create archive
        {
            let config = CompressorConfig::default();
            let mut compressor = Compressor::new(archive_path, config).unwrap();

            let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3]; // 8 bases
            compressor
                .add_contig("sample1", "chr1", seq1.clone())
                .unwrap();

            compressor.finalize().unwrap();
        }

        // Test error cases
        {
            let config = DecompressorConfig { verbosity: 0 };
            let mut decompressor = Decompressor::open(archive_path, config).unwrap();

            // Start > end
            let result = decompressor.get_contig_range("sample1", "chr1", 5, 3);
            assert!(result.is_err());

            // End > length
            let result = decompressor.get_contig_range("sample1", "chr1", 0, 100);
            assert!(result.is_err());

            // Non-existent contig
            let result = decompressor.get_contig_range("sample1", "chrX", 0, 5);
            assert!(result.is_err());

            decompressor.close().unwrap();
        }

        fs::remove_file(archive_path).unwrap();
    }
}
