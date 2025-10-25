// Contig Iterator abstraction for unified input handling
//
// Provides a common interface for reading contigs from:
// - Single pansn-format FASTA file (sample#hap#chr headers)
// - Multiple per-sample FASTA files
//
// This allows the streaming compressor to handle both input formats identically.

use crate::genome_io::GenomeIO;
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use ragc_common::Contig;
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

/// Iterator for a single pansn-format FASTA file
/// Parses sample names from headers like: >sample#hap#chromosome
pub struct PansnFileIterator {
    file_path: PathBuf,
    reader: Option<GenomeIO<Box<dyn Read>>>,
}

impl PansnFileIterator {
    /// Create a new iterator for a single pansn-format FASTA file
    pub fn new(file_path: &Path) -> Result<Self> {
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
                // IMPORTANT: Use full_header as contig name to match C++ AGC expectations
                // For PanSN format, this preserves the full "sample#hap#chr" in contig name
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
