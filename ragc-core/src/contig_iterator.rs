// Contig Iterator abstraction for unified input handling
//
// Provides a common interface for reading contigs from:
// - Single pansn-format FASTA file (sample#hap#chr headers)
//   - With sample ordering detection and indexed random access support
// - Multiple per-sample FASTA files
//
// This allows the streaming compressor to handle both input formats identically.

use crate::genome_io::GenomeIO;
use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;
use ragc_common::Contig;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

#[cfg(feature = "indexed-fasta")]
use faigz_rs::{FastaIndex, FastaReader, FastaFormat};

/// Trait for iterating over contigs from various input sources
pub trait ContigIterator {
    /// Get the next contig. Returns None when no more contigs.
    /// Returns (sample_name, contig_name, sequence)
    fn next_contig(&mut self) -> Result<Option<(String, String, Contig)>>;

    /// Reset the iterator to the beginning (for second pass)
    fn reset(&mut self) -> Result<()>;
}

/// Information about sample ordering in a PanSN file
#[derive(Debug)]
struct SampleOrderInfo {
    /// Map of sample -> list of (contig_name, file_byte_offset)
    sample_contigs: HashMap<String, Vec<String>>,
    /// Samples in the order they should be processed
    sample_order: Vec<String>,
    /// Whether samples appear in contiguous blocks (all of sample A, then all of sample B)
    is_contiguous: bool,
}

impl SampleOrderInfo {
    /// Scan a PanSN file to determine sample ordering
    /// This does a fast header-only pass through the file
    fn analyze_file(file_path: &Path) -> Result<Self> {
        let file = File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;

        let reader: Box<dyn Read> = if file_path.extension().and_then(|s| s.to_str()) == Some("gz") {
            Box::new(MultiGzDecoder::new(BufReader::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut genome_io = GenomeIO::new(reader);
        let mut sample_contigs: HashMap<String, Vec<String>> = HashMap::new();
        let mut sample_order = Vec::new();
        let mut last_sample = None;

        // Scan headers
        while let Some((full_header, sample_name, _contig_name, _sequence)) =
            genome_io.read_contig_with_sample()? {

            // Track which sample this contig belongs to
            sample_contigs
                .entry(sample_name.clone())
                .or_insert_with(Vec::new)
                .push(full_header);

            // Track when samples change (for contiguity check)
            if last_sample.as_ref() != Some(&sample_name) {
                sample_order.push(sample_name.clone());
                last_sample = Some(sample_name);
            }
        }

        // Check if samples are contiguous (each sample appears only once in order)
        let unique_samples: std::collections::HashSet<_> = sample_order.iter().collect();
        let is_contiguous = unique_samples.len() == sample_order.len();

        // If not contiguous, use sorted order for processing
        let final_sample_order = if is_contiguous {
            sample_order
        } else {
            let mut sorted: Vec<_> = sample_contigs.keys().cloned().collect();
            sorted.sort();
            sorted
        };

        Ok(SampleOrderInfo {
            sample_contigs,
            sample_order: final_sample_order,
            is_contiguous,
        })
    }
}

/// Iterator for a single pansn-format FASTA file
/// Parses sample names from headers like: >sample#hap#chromosome
pub struct PansnFileIterator {
    file_path: PathBuf,
    reader: Option<GenomeIO<Box<dyn Read>>>,
}

impl PansnFileIterator {
    /// Create a new iterator for a single pansn-format FASTA file
    ///
    /// This iterator requires samples to appear in contiguous blocks for C++ AGC compatibility.
    /// If your file has samples in non-contiguous order (e.g., all chr1, then all chr2),
    /// use one of these solutions:
    /// 1. Use IndexedPansnFileIterator with a bgzip+indexed file
    /// 2. Reorder the file by sample using `ragc sort-fasta`
    /// 3. Split into per-sample files
    pub fn new(file_path: &Path) -> Result<Self> {
        // Check if samples are contiguously ordered
        let order_info = SampleOrderInfo::analyze_file(file_path)?;

        if !order_info.is_contiguous {
            // Check if an index exists
            let index_path = format!("{}.fai", file_path.display());
            let has_index = std::path::Path::new(&index_path).exists();

            #[cfg(feature = "indexed-fasta")]
            {
                if has_index {
                    return Err(anyhow!(
                        "Samples are not contiguously ordered in file: {}\n\
                        \n\
                        For C++ AGC compatibility, samples must appear in contiguous blocks.\n\
                        \n\
                        This file has an index, so you can use random access.\n\
                        Use IndexedPansnFileIterator::new() instead of PansnFileIterator::new()\n\
                        Or enable automatic detection in your code.",
                        file_path.display()
                    ));
                }
            }

            return Err(anyhow!(
                "Samples are not contiguously ordered in file: {}\n\
                \n\
                For C++ AGC compatibility, samples must appear in contiguous blocks.\n\
                \n\
                Your file has samples scattered throughout (e.g., all chr1, then all chr2).\n\
                \n\
                Solutions:\n\
                1. Reorder file by sample:\n   \
                   ragc sort-fasta {} -o sorted.fa.gz\n\
                \n\
                2. Use bgzip compression with index for random access:\n   \
                   gunzip {}\n   \
                   bgzip {}\n   \
                   samtools faidx {}.gz\n   \
                   # Then use IndexedPansnFileIterator\n\
                \n\
                3. Split into per-sample files and compress together\n\
                \n\
                Found {} samples in non-contiguous order.",
                file_path.display(),
                file_path.display(),
                file_path.display(),
                file_path.file_stem().and_then(|s| s.to_str()).unwrap_or("input"),
                file_path.file_stem().and_then(|s| s.to_str()).unwrap_or("input"),
                order_info.sample_contigs.len()
            ));
        }

        let reader = Self::open_reader(file_path)?;
        Ok(PansnFileIterator {
            file_path: file_path.to_path_buf(),
            reader: Some(reader),
        })
    }

    fn open_reader(file_path: &Path) -> Result<GenomeIO<Box<dyn Read>>> {
        let file = File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;

        let reader: Box<dyn Read> = if file_path.extension().and_then(|s| s.to_str()) == Some("gz")
        {
            Box::new(MultiGzDecoder::new(BufReader::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(GenomeIO::new(reader))
    }
}

impl ContigIterator for PansnFileIterator {
    fn next_contig(&mut self) -> Result<Option<(String, String, Contig)>> {
        let reader = match &mut self.reader {
            Some(r) => r,
            None => return Ok(None),
        };

        match reader.read_contig_with_sample()? {
            Some((full_header, sample_name, _contig_name, sequence)) => {
                // IMPORTANT: Use full_header as contig name to preserve complete PanSN format
                // Sample name comes from parsing (e.g., "AAA#0" from "AAA#0#chrI")
                // Contig name is the full header (e.g., "AAA#0#chrI")
                // This groups contigs by sample while preserving full headers
                Ok(Some((sample_name, full_header, sequence)))
            }
            None => Ok(None),
        }
    }

    fn reset(&mut self) -> Result<()> {
        self.reader = Some(Self::open_reader(&self.file_path)?);
        Ok(())
    }
}

/// Iterator for indexed PanSN FASTA files with random access
/// Uses faigz-rs to read contigs in sample-grouped order even if file is out-of-order
#[cfg(feature = "indexed-fasta")]
pub struct IndexedPansnFileIterator {
    file_path: PathBuf,
    index: std::sync::Arc<FastaIndex>,
    order_info: SampleOrderInfo,
    current_sample_idx: usize,
    current_contig_idx: usize,
}

#[cfg(feature = "indexed-fasta")]
impl IndexedPansnFileIterator {
    /// Create a new indexed iterator for a PanSN FASTA file
    /// Requires the file to be bgzip-compressed with a .fai index
    pub fn new(file_path: &Path) -> Result<Self> {
        // Check for index file
        let index_path = format!("{}.fai", file_path.display());
        if !std::path::Path::new(&index_path).exists() {
            return Err(anyhow!(
                "Index file not found: {}\n\
                To use indexed random access, create an index with:\n  \
                samtools faidx {}",
                index_path,
                file_path.display()
            ));
        }

        // Load the index
        let index = FastaIndex::new(
            file_path.to_str().ok_or_else(|| anyhow!("Invalid path"))?,
            FastaFormat::Fasta,
        )
        .with_context(|| format!("Failed to load FASTA index for {}", file_path.display()))?;

        // Analyze sample ordering
        let order_info = SampleOrderInfo::analyze_file(file_path)?;

        Ok(IndexedPansnFileIterator {
            file_path: file_path.to_path_buf(),
            index: std::sync::Arc::new(index),
            order_info,
            current_sample_idx: 0,
            current_contig_idx: 0,
        })
    }
}

#[cfg(feature = "indexed-fasta")]
impl ContigIterator for IndexedPansnFileIterator {
    fn next_contig(&mut self) -> Result<Option<(String, String, Contig)>> {
        // Check if we've exhausted all samples
        if self.current_sample_idx >= self.order_info.sample_order.len() {
            return Ok(None);
        }

        let sample_name = &self.order_info.sample_order[self.current_sample_idx];
        let contigs = self
            .order_info
            .sample_contigs
            .get(sample_name)
            .ok_or_else(|| anyhow!("Sample not found: {}", sample_name))?;

        // Check if we've exhausted contigs for this sample
        if self.current_contig_idx >= contigs.len() {
            // Move to next sample
            self.current_sample_idx += 1;
            self.current_contig_idx = 0;
            return self.next_contig(); // Recursively get first contig of next sample
        }

        // Fetch the contig using random access
        let full_header = &contigs[self.current_contig_idx];
        self.current_contig_idx += 1;

        // Use faigz-rs to fetch the sequence
        let reader = FastaReader::new(&self.index)
            .with_context(|| "Failed to create FASTA reader")?;

        let sequence = reader
            .fetch_seq_all(full_header)
            .with_context(|| format!("Failed to fetch sequence for {}", full_header))?;

        // Convert to Contig
        let contig = Contig::from(sequence.as_bytes());

        Ok(Some((sample_name.clone(), full_header.clone(), contig)))
    }

    fn reset(&mut self) -> Result<()> {
        self.current_sample_idx = 0;
        self.current_contig_idx = 0;
        Ok(())
    }
}

/// Iterator for multiple FASTA files (one per sample)
/// Each file contains all contigs for that sample
pub struct MultiFileIterator {
    file_paths: Vec<PathBuf>,
    current_file_idx: usize,
    current_reader: Option<GenomeIO<Box<dyn Read>>>,
    current_sample_name: String,
}

impl MultiFileIterator {
    /// Create a new iterator for multiple FASTA files
    /// Each file is treated as a separate sample
    pub fn new(file_paths: Vec<PathBuf>) -> Result<Self> {
        let mut iterator = MultiFileIterator {
            file_paths,
            current_file_idx: 0,
            current_reader: None,
            current_sample_name: String::new(),
        };

        // Open first file
        if !iterator.file_paths.is_empty() {
            iterator.advance_to_next_file()?;
        }

        Ok(iterator)
    }

    fn open_file(&self, file_path: &Path) -> Result<(GenomeIO<Box<dyn Read>>, String)> {
        let file = File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;

        let reader: Box<dyn Read> = if file_path.extension().and_then(|s| s.to_str()) == Some("gz")
        {
            Box::new(MultiGzDecoder::new(BufReader::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        // Extract sample name from filename
        let sample_name = file_path
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| {
                // Remove .fa or .fasta extensions if present
                s.trim_end_matches(".fa")
                    .trim_end_matches(".fasta")
                    .to_string()
            })
            .unwrap_or_else(|| "unknown".to_string());

        Ok((GenomeIO::new(reader), sample_name))
    }

    fn advance_to_next_file(&mut self) -> Result<bool> {
        if self.current_file_idx >= self.file_paths.len() {
            self.current_reader = None;
            return Ok(false);
        }

        let file_path = &self.file_paths[self.current_file_idx];
        let (reader, sample_name) = self.open_file(file_path)?;

        self.current_reader = Some(reader);
        self.current_sample_name = sample_name;
        self.current_file_idx += 1;

        Ok(true)
    }
}

impl ContigIterator for MultiFileIterator {
    fn next_contig(&mut self) -> Result<Option<(String, String, Contig)>> {
        loop {
            let reader = match &mut self.current_reader {
                Some(r) => r,
                None => return Ok(None),
            };

            match reader.read_contig_with_sample()? {
                Some((full_header, sample_from_header, _contig_name, sequence)) => {
                    // Use sample name from header (PanSN format) if available, else fallback to filename
                    let sample_name = if sample_from_header != "unknown" {
                        sample_from_header
                    } else {
                        self.current_sample_name.clone()
                    };
                    // IMPORTANT: Use full_header as contig name to match C++ AGC expectations
                    return Ok(Some((sample_name, full_header, sequence)));
                }
                None => {
                    // Current file exhausted, move to next file
                    if !self.advance_to_next_file()? {
                        return Ok(None);
                    }
                    // Continue loop to read from new file
                }
            }
        }
    }

    fn reset(&mut self) -> Result<()> {
        self.current_file_idx = 0;
        self.current_reader = None;
        self.current_sample_name.clear();

        if !self.file_paths.is_empty() {
            self.advance_to_next_file()?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_pansn_file_iterator() {
        // Create a temporary pansn-format FASTA file
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, ">sample1#1#chr1").unwrap();
        writeln!(temp_file, "ACGTACGT").unwrap();
        writeln!(temp_file, ">sample1#1#chr2").unwrap();
        writeln!(temp_file, "TGCATGCA").unwrap();
        writeln!(temp_file, ">sample2#0#chr1").unwrap();
        writeln!(temp_file, "AAAACCCC").unwrap();
        temp_file.flush().unwrap();

        let mut iterator = PansnFileIterator::new(temp_file.path()).unwrap();

        // First contig (contig name is now the full header for C++ AGC compatibility)
        let (sample, contig, _seq) = iterator.next_contig().unwrap().unwrap();
        assert_eq!(sample, "sample1#1");
        assert_eq!(contig, "sample1#1#chr1");

        // Second contig
        let (sample, contig, _seq) = iterator.next_contig().unwrap().unwrap();
        assert_eq!(sample, "sample1#1");
        assert_eq!(contig, "sample1#1#chr2");

        // Third contig
        let (sample, contig, _seq) = iterator.next_contig().unwrap().unwrap();
        assert_eq!(sample, "sample2#0");
        assert_eq!(contig, "sample2#0#chr1");

        // No more contigs
        assert!(iterator.next_contig().unwrap().is_none());

        // Test reset
        iterator.reset().unwrap();
        let (sample, contig, _seq) = iterator.next_contig().unwrap().unwrap();
        assert_eq!(sample, "sample1#1");
        assert_eq!(contig, "sample1#1#chr1");
    }
}
