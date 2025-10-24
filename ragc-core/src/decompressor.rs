// AGC Decompressor
// Extracts genomes from AGC archives

use crate::{
    genome_io::GenomeWriter, kmer::reverse_complement, lz_diff::LZDiff,
    segment_compression::decompress_segment,
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
pub struct Decompressor {
    config: DecompressorConfig,
    archive: Archive,
    collection: CollectionV3,

    // Cached segment data (group_id -> reference segment)
    segment_cache: HashMap<u32, Contig>,

    // Archive parameters
    _segment_size: u32,
    kmer_length: u32,
    min_match_len: u32,
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

    /// Get list of contigs for a specific sample
    pub fn list_contigs(&mut self, sample_name: &str) -> Result<Vec<String>> {
        if self.config.verbosity > 1 {
            eprintln!("list_contigs called for: {sample_name}");
            eprintln!(
                "get_no_contigs before load: {:?}",
                self.collection.get_no_contigs(sample_name)
            );
        }

        // Load contig batch if not already loaded (either None or Some(0))
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            if self.config.verbosity > 1 {
                eprintln!("Loading contig batch for sample: {sample_name}");
            }
            self.collection.load_contig_batch(&mut self.archive, 0)?;

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

    /// Extract a specific contig from a sample
    pub fn get_contig(&mut self, sample_name: &str, contig_name: &str) -> Result<Contig> {
        // Load contig batch if needed (either None or Some(0))
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            let _num_samples = self.collection.get_no_samples();
            self.collection.load_contig_batch(&mut self.archive, 0)?;
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

    /// Extract all contigs from a sample
    pub fn get_sample(&mut self, sample_name: &str) -> Result<Vec<(String, Contig)>> {
        // Load contig batch if needed (either None or Some(0))
        if self
            .collection
            .get_no_contigs(sample_name)
            .is_none_or(|count| count == 0)
        {
            let _num_samples = self.collection.get_no_samples();
            self.collection.load_contig_batch(&mut self.archive, 0)?;
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
                    let _marker = ref_data.pop().unwrap();
                    decompress_segment(&ref_data)?
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
            let delta_stream_id = self
                .archive
                .get_stream_id(&delta_stream_name)
                .ok_or_else(|| anyhow!("Delta stream not found: {delta_stream_name}"))?;

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
                let _marker = delta_data.pop().unwrap();
                decompress_segment(&delta_data)?
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
                let _marker = data.pop().unwrap();
                decompress_segment(&data)?
            };

            // Unpack the raw segment
            let contig_data = Self::unpack_contig(&decompressed_pack, position_in_pack)?;

            if self.config.verbosity > 1 {
                eprintln!("Unpacked raw segment: length={}", contig_data.len());
            }

            Ok(contig_data)
        }
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
}
