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
use ragc_core::Contig;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

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
    /// Analyze sample ordering from FAI index (fast - no decompression needed)
    fn analyze_from_index(index_path: &str) -> Result<Self> {
        use std::io::BufRead;

        let file = File::open(index_path)
            .with_context(|| format!("Failed to open index: {index_path}"))?;
        let reader = BufReader::new(file);

        let mut sample_contigs: HashMap<String, Vec<String>> = HashMap::new();
        let mut sample_order = Vec::new();
        let mut last_sample = None;

        // Read headers from first column of FAI file
        for line in reader.lines() {
            let line = line?;
            let full_header = line
                .split('\t')
                .next()
                .ok_or_else(|| anyhow!("Invalid FAI format"))?
                .to_string();

            // Parse sample name (sample#hap format)
            let sample_name = if let Some(parts) = full_header.split('#').nth(0) {
                if let Some(hap) = full_header.split('#').nth(1) {
                    format!("{parts}#{hap}")
                } else {
                    full_header.clone()
                }
            } else {
                full_header.clone()
            };

            // Track which sample this contig belongs to
            sample_contigs
                .entry(sample_name.clone())
                .or_default()
                .push(full_header);

            // Track when samples change (for contiguity check)
            if last_sample.as_ref() != Some(&sample_name) {
                sample_order.push(sample_name.clone());
                last_sample = Some(sample_name);
            }
        }

        // Check if samples are contiguous
        let unique_samples: std::collections::HashSet<_> = sample_order.iter().collect();
        let is_contiguous = unique_samples.len() == sample_order.len();

        // If not contiguous, use sorted order
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

    /// Scan a PanSN file to determine sample ordering
    /// This does a fast header-only pass through the file
    fn analyze_file(file_path: &Path) -> Result<Self> {
        // Check if FAI index exists - if so, use it for instant header reading
        let index_path = format!("{}.fai", file_path.display());
        if std::path::Path::new(&index_path).exists() {
            return Self::analyze_from_index(&index_path);
        }

        // Otherwise fall back to decompressing (slow)
        let file = File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;

        let reader: Box<dyn Read> = if file_path.extension().and_then(|s| s.to_str()) == Some("gz")
        {
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
            genome_io.read_contig_with_sample()?
        {
            // Track which sample this contig belongs to
            sample_contigs
                .entry(sample_name.clone())
                .or_default()
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

            let _ = has_index; // suppress unused warning
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
                file_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("input"),
                file_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("input"),
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

/// Iterator that buffers entire file in memory and outputs in sample-grouped order
/// Much faster than random access for files with non-contiguous sample ordering
pub struct BufferedPansnFileIterator {
    // Map from sample name to vector of (contig_name, sequence)
    sample_contigs: HashMap<String, Vec<(String, Contig)>>,
    sample_order: Vec<String>,
    current_sample_idx: usize,
    current_contig_idx: usize,
}

impl BufferedPansnFileIterator {
    /// Create a new buffered iterator that reads entire file into memory
    pub fn new(file_path: &Path) -> Result<Self> {
        eprintln!("Reading entire file into memory for reordering...");

        // Open reader
        let file = File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
        let reader: Box<dyn Read> = if file_path.to_string_lossy().ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        let mut genome_io = GenomeIO::new(reader);

        // Read all contigs and group by sample
        let mut sample_contigs: HashMap<String, Vec<(String, Contig)>> = HashMap::new();
        let mut sample_order = Vec::new();
        let mut seen_samples = std::collections::HashSet::new();

        while let Some((header, contig)) = genome_io.read_contig_converted()? {
            // Parse sample name from header (sample#hap#chr format)
            let sample_name = if let Some(parts) = header.split('#').nth(0) {
                if let Some(hap) = header.split('#').nth(1) {
                    format!("{parts}#{hap}")
                } else {
                    header.clone()
                }
            } else {
                header.clone()
            };

            // Track sample order (first occurrence)
            if !seen_samples.contains(&sample_name) {
                sample_order.push(sample_name.clone());
                seen_samples.insert(sample_name.clone());
            }

            // Store contig
            sample_contigs
                .entry(sample_name)
                .or_default()
                .push((header, contig));
        }

        eprintln!("Loaded {} samples into memory", sample_order.len());

        Ok(BufferedPansnFileIterator {
            sample_contigs,
            sample_order,
            current_sample_idx: 0,
            current_contig_idx: 0,
        })
    }
}

impl ContigIterator for BufferedPansnFileIterator {
    fn next_contig(&mut self) -> Result<Option<(String, String, Contig)>> {
        // Check if we've exhausted all samples
        if self.current_sample_idx >= self.sample_order.len() {
            return Ok(None);
        }

        let sample_name = &self.sample_order[self.current_sample_idx];
        let contigs = self
            .sample_contigs
            .get(sample_name)
            .ok_or_else(|| anyhow!("Sample not found: {sample_name}"))?;

        // Check if we've exhausted contigs for this sample
        if self.current_contig_idx >= contigs.len() {
            // Move to next sample
            self.current_sample_idx += 1;
            self.current_contig_idx = 0;
            return self.next_contig(); // Recursively get first contig of next sample
        }

        // Get the contig
        let (contig_name, contig) = &contigs[self.current_contig_idx];
        self.current_contig_idx += 1;

        Ok(Some((
            sample_name.clone(),
            contig_name.clone(),
            contig.clone(),
        )))
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
