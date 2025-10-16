// AGC Decompressor
// Extracts genomes from AGC archives

use crate::{genome_io::GenomeWriter, lz_diff::LZDiff, segment_compression::decompress_segment};
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
}

impl Decompressor {
    /// Open an existing archive for decompression
    pub fn open(archive_path: &str, config: DecompressorConfig) -> Result<Self> {
        let mut archive = Archive::new_reader();
        archive
            .open(archive_path)
            .context("Failed to open archive for reading")?;

        let mut collection = CollectionV3::new();

        // Load segment_size and kmer_length from params stream
        let (segment_size, kmer_length) = Self::load_params(&mut archive)?;

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
        })
    }

    /// Load archive parameters from the params stream
    fn load_params(archive: &mut Archive) -> Result<(u32, u32)> {
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

        // Parse 4 little-endian u32 values:
        // kmer_length, min_match_len, pack_cardinality, segment_size
        if data.len() < 16 {
            anyhow::bail!(
                "params stream too short: {} bytes (expected 16)",
                data.len()
            );
        }

        let kmer_length = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
        let _min_match_len = u32::from_le_bytes([data[4], data[5], data[6], data[7]]);
        let _pack_cardinality = u32::from_le_bytes([data[8], data[9], data[10], data[11]]);
        let segment_size = u32::from_le_bytes([data[12], data[13], data[14], data[15]]);

        Ok((segment_size, kmer_length))
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

    /// Reconstruct a contig from its segment descriptors
    fn reconstruct_contig(&mut self, segments: &[SegmentDesc]) -> Result<Contig> {
        let mut contig = Contig::new();

        for (i, segment_desc) in segments.iter().enumerate() {
            let segment_data = self.get_segment(segment_desc)?;

            if i == 0 {
                // First segment: add completely
                contig.extend_from_slice(&segment_data);
            } else {
                // Subsequent segments: skip first kmer_length bases (overlap with previous segment)
                if segment_data.len() < self.kmer_length as usize {
                    anyhow::bail!("Corrupted archive: segment too short");
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
    fn get_segment(&mut self, desc: &SegmentDesc) -> Result<Contig> {
        const PACK_CARDINALITY: usize = 50; // C++ default

        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;

        if self.config.verbosity > 1 {
            eprintln!(
                "get_segment: group_id={}, in_group_id={}, raw_length={}",
                desc.group_id, desc.in_group_id, desc.raw_length
            );
        }

        // Calculate which pack and position within pack
        let pack_id = desc.in_group_id as usize / PACK_CARDINALITY;
        let position_in_pack = desc.in_group_id as usize % PACK_CARDINALITY;

        if self.config.verbosity > 1 {
            eprintln!("  pack_id={pack_id}, position_in_pack={position_in_pack}");
        }

        // Get the pack containing our segment
        let stream_name = stream_delta_name(archive_version, desc.group_id);
        let stream_id = self.archive.get_stream_id(&stream_name).ok_or_else(|| {
            // Try reference stream as fallback
            let ref_stream_name = stream_ref_name(archive_version, desc.group_id);
            if let Some(_ref_id) = self.archive.get_stream_id(&ref_stream_name) {
                return anyhow!(
                    "Found ref stream but not delta stream for group {}",
                    desc.group_id
                );
            }
            anyhow!("Delta stream not found: {stream_name}")
        })?;

        // Fetch pack at pack_id
        let (mut compressed, _) = self.archive.get_part_by_id(stream_id, pack_id)?;

        // Remove marker byte
        if compressed.is_empty() {
            anyhow::bail!("Empty compressed data for segment");
        }
        let _marker = compressed.pop().unwrap();

        // Decompress the pack
        let decompressed_pack = decompress_segment(&compressed)?;

        if self.config.verbosity > 1 {
            eprintln!(
                "Decompressed pack (group {}, pack {}): length={}",
                desc.group_id,
                pack_id,
                decompressed_pack.len()
            );
            eprintln!(
                "Pack bytes: {:?}",
                &decompressed_pack[..decompressed_pack.len().min(50)]
            );
        }

        // Unpack to extract the specific contig
        let contig_data = Self::unpack_contig(&decompressed_pack, position_in_pack)?;

        if self.config.verbosity > 1 {
            eprintln!(
                "Unpacked contig at position {}: length={}, first 20 bytes: {:?}",
                position_in_pack,
                contig_data.len(),
                &contig_data[..contig_data.len().min(20)]
            );
        }

        // Determine if this needs LZ decoding
        // Key insights:
        // 1. C++ has 16 "raw-only" groups (0-15) where ALL segments are raw (no LZ encoding)
        // 2. Groups 16+ can use LZ encoding with first segment as reference
        // 3. C++ vs Rust archives differ:
        //    - C++: Each contig often in different groups, with in_group_id=1 (raw)
        //    - Rust: Multiple contigs in same group
        //
        // Detection strategy:
        // - If group_id < 16: Always raw (C++ raw-only groups)
        // - If in_group_id=0: Always raw (reference), cache it
        // - If in_group_id>=1 in groups 16+: Check if group already has cached reference
        //   - Has reference: This is LZ-encoded (Rust format with LZ)
        //   - No reference: This is raw data (C++ format), cache it

        const NO_RAW_GROUPS: u32 = 16; // C++ constant

        if desc.group_id < NO_RAW_GROUPS {
            // Groups 0-15 are raw-only, never apply LZ decoding
            // Cache first segment as reference for consistency
            if desc.in_group_id == 0 {
                self.segment_cache
                    .insert(desc.group_id, contig_data.clone());
            }
            Ok(contig_data)
        } else if desc.in_group_id == 0 {
            // Position 0 in groups 16+ is raw reference data
            self.segment_cache
                .insert(desc.group_id, contig_data.clone());
            Ok(contig_data)
        } else if let std::collections::hash_map::Entry::Vacant(e) =
            self.segment_cache.entry(desc.group_id)
        {
            // No cached reference for this group, so this must be raw data (C++ format)
            // Cache it as the reference
            e.insert(contig_data.clone());
            Ok(contig_data)
        } else {
            // We have a reference cached for this group (16+), so this must be LZ-encoded
            let reference = self.segment_cache.get(&desc.group_id).ok_or_else(|| {
                anyhow!("Reference segment not loaded for group {}", desc.group_id)
            })?;

            if self.config.verbosity > 1 {
                eprintln!(
                    "Applying LZ decoding with reference: length={}, first 20 bytes: {:?}",
                    reference.len(),
                    &reference[..reference.len().min(20)]
                );
                eprintln!(
                    "Encoded data: length={}, first 20 bytes: {:?}",
                    contig_data.len(),
                    &contig_data[..contig_data.len().min(20)]
                );
            }

            // Apply LZ diff decoding
            let mut lz_diff = LZDiff::new(15);
            lz_diff.prepare(reference);
            let reconstructed = lz_diff.decode(&contig_data);

            if self.config.verbosity > 1 {
                eprintln!(
                    "Reconstructed: length={}, first 20 bytes: {:?}",
                    reconstructed.len(),
                    &reconstructed[..reconstructed.len().min(20)]
                );
            }

            Ok(reconstructed)
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
