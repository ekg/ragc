// AGC Compressor
// Orchestrates the compression pipeline: FASTA → segments → compression → archive

use crate::{
    genome_io::GenomeIO,
    kmer::{Kmer, KmerMode},
    lz_diff::LZDiff,
    segment_compression::compress_segment,
};
use anyhow::{Context, Result};
use ragc_common::{
    stream_delta_name, Archive, CollectionV3, Contig, AGC_FILE_MAJOR, AGC_FILE_MINOR,
    CONTIG_SEPARATOR,
};
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

/// Configuration for the compressor
#[derive(Debug, Clone)]
pub struct CompressorConfig {
    pub kmer_length: u32,
    pub segment_size: u32,
    pub min_match_len: u32,
    pub verbosity: u32,
}

impl Default for CompressorConfig {
    fn default() -> Self {
        CompressorConfig {
            kmer_length: 21,
            segment_size: 1000,
            min_match_len: 15,
            verbosity: 1,
        }
    }
}

/// A segment with its metadata
#[derive(Debug, Clone)]
struct SegmentInfo {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    data: Contig,
    #[allow(dead_code)]
    kmer_front: u64,
    #[allow(dead_code)]
    kmer_back: u64,
    is_rev_comp: bool,
}

/// Segment group identified by flanking k-mers
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct SegmentGroupKey {
    kmer_front: u64,
    kmer_back: u64,
}

/// AGC Compressor
pub struct Compressor {
    config: CompressorConfig,
    archive: Archive,
    collection: CollectionV3,

    // Segment storage
    segment_groups: HashMap<SegmentGroupKey, Vec<SegmentInfo>>,
    group_references: HashMap<SegmentGroupKey, Vec<u8>>, // First segment in each group

    // Statistics
    total_bases_processed: usize,
    total_segments: usize,
}

impl Compressor {
    /// Create a new compressor
    pub fn new(archive_path: &str, config: CompressorConfig) -> Result<Self> {
        let mut archive = Archive::new_writer();
        archive
            .open(archive_path)
            .context("Failed to create archive")?;

        let mut collection = CollectionV3::new();
        collection.set_config(config.segment_size, config.kmer_length, None);
        collection.prepare_for_compression(&mut archive)?;

        Ok(Compressor {
            config,
            archive,
            collection,
            segment_groups: HashMap::new(),
            group_references: HashMap::new(),
            total_bases_processed: 0,
            total_segments: 0,
        })
    }

    /// Add a FASTA file to the archive
    pub fn add_fasta_file(&mut self, sample_name: &str, fasta_path: &Path) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Processing sample: {sample_name} from {fasta_path:?}");
        }

        let mut reader =
            GenomeIO::<Box<dyn Read>>::open(fasta_path).context("Failed to open FASTA file")?;

        // Read contigs with conversion (ASCII -> numeric)
        while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
            self.add_contig(sample_name, &contig_name, sequence)?;
        }

        Ok(())
    }

    /// Add a contig to the archive
    pub fn add_contig(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        sequence: Contig,
    ) -> Result<()> {
        if sequence.is_empty() {
            return Ok(());
        }

        self.total_bases_processed += sequence.len();

        // Register in collection
        self.collection
            .register_sample_contig(sample_name, contig_name)?;

        // For now, treat entire contig as a single segment
        // TODO: Add splitter-based segmentation later
        self.add_segment(sample_name, contig_name, 0, sequence)?;

        Ok(())
    }

    /// Add a segment to the compressor
    #[allow(clippy::needless_range_loop)]
    fn add_segment(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_part_no: usize,
        segment: Contig,
    ) -> Result<()> {
        // Extract flanking k-mers
        let k = self.config.kmer_length;

        let (kmer_front, kmer_back, is_rev_comp) = if segment.len() >= (k * 2) as usize {
            // Extract front k-mer
            let mut front = Kmer::new(k, KmerMode::Canonical);
            for i in 0..(k as usize) {
                if segment[i] > 3 {
                    // Invalid base, use 0
                    front.reset();
                    break;
                }
                front.insert(segment[i] as u64);
            }
            let front_kmer_val = if front.is_full() { front.data() } else { 0 };

            // Extract back k-mer
            let mut back = Kmer::new(k, KmerMode::Canonical);
            let start = segment.len() - (k as usize);
            for i in 0..(k as usize) {
                if segment[start + i] > 3 {
                    // Invalid base, use 0
                    back.reset();
                    break;
                }
                back.insert(segment[start + i] as u64);
            }
            let back_kmer_val = if back.is_full() { back.data() } else { 0 };

            (front_kmer_val, back_kmer_val, false)
        } else {
            // Too short for proper k-mers
            (0u64, 0u64, false)
        };

        self.add_segment_with_kmers(
            sample_name,
            contig_name,
            seg_part_no,
            segment,
            kmer_front,
            kmer_back,
            is_rev_comp,
        )
    }

    /// Add a segment with known k-mers
    #[allow(clippy::too_many_arguments)]
    fn add_segment_with_kmers(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        seg_part_no: usize,
        segment: Contig,
        kmer_front: u64,
        kmer_back: u64,
        is_rev_comp: bool,
    ) -> Result<()> {
        let key = SegmentGroupKey {
            kmer_front,
            kmer_back,
        };

        let seg_info = SegmentInfo {
            sample_name: sample_name.to_string(),
            contig_name: contig_name.to_string(),
            seg_part_no,
            data: segment,
            kmer_front,
            kmer_back,
            is_rev_comp,
        };

        self.segment_groups.entry(key).or_default().push(seg_info);
        self.total_segments += 1;

        Ok(())
    }

    /// Store params stream (C++ compatibility)
    /// Format: kmer_length (u32) + min_match_len (u32) + pack_cardinality (u32) + segment_size (u32)
    fn store_params_stream(&mut self) -> Result<()> {
        let mut params_data = Vec::new();

        // Append u32 values in little-endian format (matching C++ append function)
        let append_u32 = |data: &mut Vec<u8>, value: u32| {
            data.extend_from_slice(&value.to_le_bytes());
        };

        append_u32(&mut params_data, self.config.kmer_length);
        append_u32(&mut params_data, self.config.min_match_len);
        append_u32(&mut params_data, 50); // pack_cardinality - default 50 (TODO: make configurable)
        append_u32(&mut params_data, self.config.segment_size);

        let stream_id = self.archive.register_stream("params");
        self.archive
            .add_part(stream_id, &params_data, params_data.len() as u64)?;

        Ok(())
    }

    /// Store empty splitters stream (C++ compatibility)
    fn store_splitters_stream(&mut self) -> Result<()> {
        // For now, store empty splitters (no segmentation yet)
        let splitters_data = Vec::new();
        let stream_id = self.archive.register_stream("splitters");
        self.archive.add_part(stream_id, &splitters_data, 0)?;
        Ok(())
    }

    /// Store empty segment-splitters stream (C++ compatibility)
    fn store_segment_splitters_stream(&mut self) -> Result<()> {
        // For now, store empty segment-splitters (no segmentation yet)
        let seg_splitters_data = Vec::new();
        let stream_id = self.archive.register_stream("segment-splitters");
        self.archive.add_part(stream_id, &seg_splitters_data, 0)?;
        Ok(())
    }

    /// Store file_type_info stream (C++ compatibility)
    /// Format: null-terminated string pairs (key, value, key, value, ...)
    /// C++ reads this first to determine archive version
    fn store_file_type_info(&mut self) -> Result<()> {
        let mut data = Vec::new();

        // Helper to append null-terminated string
        let append_str = |data: &mut Vec<u8>, s: &str| {
            data.extend_from_slice(s.as_bytes());
            data.push(0); // null terminator
        };

        // Add key-value pairs (C++ expects these specific keys)
        append_str(&mut data, "producer");
        append_str(&mut data, "agc-rust");

        append_str(&mut data, "producer_version_major");
        append_str(&mut data, &AGC_FILE_MAJOR.to_string());

        append_str(&mut data, "producer_version_minor");
        append_str(&mut data, &AGC_FILE_MINOR.to_string());

        append_str(&mut data, "producer_version_build");
        append_str(&mut data, "0"); // TODO: Add actual build version

        append_str(&mut data, "file_version_major");
        append_str(&mut data, &AGC_FILE_MAJOR.to_string());

        append_str(&mut data, "file_version_minor");
        append_str(&mut data, &AGC_FILE_MINOR.to_string());

        append_str(&mut data, "comment");
        append_str(
            &mut data,
            &format!("AGC (Rust implementation) v.{AGC_FILE_MAJOR}.{AGC_FILE_MINOR}"),
        );

        let stream_id = self.archive.register_stream("file_type_info");
        // raw_size = number of key-value pairs (7)
        self.archive.add_part(stream_id, &data, 7)?;

        Ok(())
    }

    /// Finalize compression and write all segments
    pub fn finalize(&mut self) -> Result<()> {
        if self.config.verbosity > 0 {
            println!("Finalizing compression...");
            println!("Total segments: {}", self.total_segments);
            println!("Segment groups: {}", self.segment_groups.len());
        }

        // Process each segment group with packed-contig mode
        let mut group_id = 0u32;
        const PACK_CARDINALITY: usize = 50; // C++ default
        const NO_RAW_GROUPS: u32 = 16; // C++ expects first 16 groups to be raw-only

        for (key, segments) in &self.segment_groups {
            if segments.is_empty() {
                continue;
            }

            // Use first segment as reference (only if using LZ encoding)
            let reference = &segments[0].data;
            self.group_references.insert(key.clone(), reference.clone());

            let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
            let stream_name = stream_delta_name(archive_version, group_id);
            let stream_id = self.archive.register_stream(&stream_name);

            // Determine if this group should use LZ encoding
            // Groups 0-15 must be raw-only for C++ compatibility
            let use_lz_encoding = group_id >= NO_RAW_GROUPS;

            // Pack contigs into batches of PACK_CARDINALITY (50)
            for (pack_idx, pack_segments) in segments.chunks(PACK_CARDINALITY).enumerate() {
                let mut packed_data = Vec::new();

                // Process segments (all raw for groups 0-15, LZ for groups 16+)
                for (idx_in_pack, seg_info) in pack_segments.iter().enumerate() {
                    let global_in_group_id = pack_idx * PACK_CARDINALITY + idx_in_pack;

                    let contig_data = if !use_lz_encoding || global_in_group_id == 0 {
                        // Raw segment (either in raw-only group, or first segment in LZ group)
                        seg_info.data.clone()
                    } else {
                        // LZ-encoded segment (only for groups 16+ and not first segment)
                        let mut lz_diff = LZDiff::new(self.config.min_match_len);
                        lz_diff.prepare(reference);
                        lz_diff.encode(&seg_info.data)
                    };

                    // Add contig data
                    packed_data.extend_from_slice(&contig_data);
                    // Add separator (0xFF)
                    packed_data.push(CONTIG_SEPARATOR);

                    // Register in collection
                    self.collection.add_segment_placed(
                        &seg_info.sample_name,
                        &seg_info.contig_name,
                        seg_info.seg_part_no,
                        group_id,
                        global_in_group_id as u32,
                        seg_info.is_rev_comp,
                        seg_info.data.len() as u32,
                    )?;
                }

                // Compress and store the packed data
                // The raw size must include separators (one per segment)
                let total_raw_size = packed_data.len() as u64;

                let mut compressed = compress_segment(&packed_data)?;
                compressed.push(0); // Marker byte

                self.archive
                    .add_part(stream_id, &compressed, total_raw_size)?;
            }

            group_id += 1;
        }

        // Store params stream (C++ compatibility)
        self.store_params_stream()?;

        // Store empty splitters stream (C++ compatibility)
        self.store_splitters_stream()?;

        // Store empty segment-splitters stream (C++ compatibility)
        self.store_segment_splitters_stream()?;

        // Store collection metadata
        self.collection
            .store_batch_sample_names(&mut self.archive)?;

        // Store all samples in one batch
        let num_samples = self.collection.get_no_samples();
        if num_samples > 0 {
            self.collection
                .store_contig_batch(&mut self.archive, 0, num_samples)?;
        }

        // Store file_type_info stream (must be last, C++ reads this first)
        self.store_file_type_info()?;

        // Close archive
        self.archive.close()?;

        if self.config.verbosity > 0 {
            println!("Compression complete!");
            println!("Total bases: {}", self.total_bases_processed);
            println!("Groups created: {group_id}");
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_compressor_basic() {
        let archive_path = "/tmp/test_compress.agc";
        let _ = fs::remove_file(archive_path); // Clean up if exists

        let config = CompressorConfig::default();
        let mut compressor = Compressor::new(archive_path, config).unwrap();

        // Add a simple contig (numeric encoding: A=0, C=1, G=2, T=3)
        let sequence = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        compressor.add_contig("sample1", "chr1", sequence).unwrap();

        compressor.finalize().unwrap();

        // Verify archive was created
        assert!(Path::new(archive_path).exists());

        fs::remove_file(archive_path).unwrap();
    }

    #[test]
    fn test_compressor_multiple_samples() {
        let archive_path = "/tmp/test_compress_multi.agc";
        let _ = fs::remove_file(archive_path);

        let config = CompressorConfig::default();
        let mut compressor = Compressor::new(archive_path, config).unwrap();

        // Add multiple samples and contigs
        let seq1 = vec![0, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        let seq2 = vec![3, 2, 1, 0, 3, 2, 1, 0]; // TGCATGCA

        compressor
            .add_contig("sample1", "chr1", seq1.clone())
            .unwrap();
        compressor
            .add_contig("sample1", "chr2", seq2.clone())
            .unwrap();
        compressor.add_contig("sample2", "chr1", seq1).unwrap();

        compressor.finalize().unwrap();

        assert!(Path::new(archive_path).exists());

        // Verify we can read the collection back
        let mut archive = Archive::new_reader();
        archive.open(archive_path).unwrap();

        let mut collection = CollectionV3::new();
        collection.set_config(1000, 21, None);
        collection.prepare_for_decompression(&archive).unwrap();
        collection.load_batch_sample_names(&mut archive).unwrap();

        assert_eq!(collection.get_no_samples(), 2);
        let samples = collection.get_samples_list(false);
        assert_eq!(samples.len(), 2);
        assert!(samples.contains(&"sample1".to_string()));
        assert!(samples.contains(&"sample2".to_string()));

        fs::remove_file(archive_path).unwrap();
    }

    #[test]
    fn test_compressor_from_fasta() {
        let archive_path = "/tmp/test_compress_fasta.agc";
        let _ = fs::remove_file(archive_path);

        let config = CompressorConfig::default();
        let mut compressor = Compressor::new(archive_path, config).unwrap();

        // Use the test FASTA file
        let fasta_path = Path::new("../test-data/test_simple.fasta");

        if fasta_path.exists() {
            compressor
                .add_fasta_file("test_sample", fasta_path)
                .unwrap();
            compressor.finalize().unwrap();

            assert!(Path::new(archive_path).exists());

            // Verify collection metadata
            let mut archive = Archive::new_reader();
            archive.open(archive_path).unwrap();

            let mut collection = CollectionV3::new();
            collection.set_config(1000, 21, None);
            collection.prepare_for_decompression(&archive).unwrap();
            collection.load_batch_sample_names(&mut archive).unwrap();

            assert_eq!(collection.get_no_samples(), 1);
            let samples = collection.get_samples_list(false);
            assert_eq!(samples, vec!["test_sample"]);

            collection.load_contig_batch(&mut archive, 0).unwrap();
            assert_eq!(collection.get_no_contigs("test_sample"), Some(2));

            let contigs = collection.get_contig_list("test_sample").unwrap();
            assert_eq!(contigs.len(), 2);
            assert!(contigs.contains(&"chr1".to_string()));
            assert!(contigs.contains(&"chr2".to_string()));

            fs::remove_file(archive_path).unwrap();
        }
    }
}
