// FASTA/FASTQ genome I/O
// Rust equivalent of genome_io.cpp

#![allow(dead_code)]

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use ragc_common::Contig;

/// Nucleotide conversion table matching C++ cnv_num
/// Maps ASCII characters to numeric representation
/// Index 0-15 maps numeric encoding to ASCII: CNV_NUM[0]='A', CNV_NUM[1]='C', etc.
/// Index 65+ maps ASCII to numeric: CNV_NUM['A']='0', CNV_NUM['C']='1', etc.
pub const CNV_NUM: [u8; 128] = [
    // 0-15: Control characters -> ACGT... (matching C++)
    b'A', b'C', b'G', b'T', b'N', b'R', b'Y', b'S', b'W', b'K', b'M', b'B', b'D', b'H', b'V', b'U',
    // 16-63: Unused
    b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
    b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
    b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
    // 64-79: Upper case letters starting at 'A' (65)
    b' ', 0, 11, 1, 12, 30, 30, 2, 13, 30, 30, 9, 30, 10, 4, 30,
    // 80-95: Upper case continued
    30, 30, 5, 7, 3, 15, 14, 8, 30, 6, 30, 30, 30, 30, 30, 30,
    // 96-111: Lower case letters starting at 'a' (97)
    b' ', 0, 11, 1, 12, 30, 30, 2, 13, 30, 30, 9, 30, 10, 4, 30,
    // 112-127: Lower case continued
    30, 30, 5, 7, 3, 15, 14, 8, 30, 6, 30, 30, 30, 30, 30, 30,
];

/// Parse sample name and contig name from FASTA header
///
/// Format: >sample#haplotype#chromosome_part or >sample#haplotype#chromosome
/// Examples:
///   "S288C#1#chrI" -> ("S288C#1", "chrI")
///   "AAA#0#chrI_1" -> ("AAA#0", "chrI_1")
///   "simple_name" -> ("unknown", "simple_name")
///
pub fn parse_sample_from_header(header: &str) -> (String, String) {
    // Split by '#'
    let parts: Vec<&str> = header.split('#').collect();

    if parts.len() >= 3 {
        // Format: sample#haplotype#chromosome
        let sample_name = format!("{}#{}", parts[0], parts[1]);
        let contig_name = parts[2..].join("#"); // Join remaining parts (handles chr#extra)
        (sample_name, contig_name)
    } else {
        // No multi-sample format detected, use entire header as contig name
        ("unknown".to_string(), header.to_string())
    }
}

/// FASTA/FASTQ parser
pub struct GenomeIO<R> {
    reader: Option<BufReader<R>>,
    buffer: Vec<u8>,
    next_header: Option<Vec<u8>>, // Buffered header line when we read ahead
}

impl<R: Read> GenomeIO<R> {
    /// Create a new genome reader
    pub fn new(reader: R) -> Self {
        GenomeIO {
            reader: Some(BufReader::with_capacity(4 << 20, reader)),
            buffer: Vec::with_capacity(4 << 20),
            next_header: None,
        }
    }

    /// Read next contig without conversion (preserves ASCII)
    pub fn read_contig(&mut self) -> io::Result<Option<(String, Contig)>> {
        self.read_contig_impl(false)
    }

    /// Read next contig with nucleotide conversion (ASCII -> numeric)
    pub fn read_contig_converted(&mut self) -> io::Result<Option<(String, Contig)>> {
        self.read_contig_impl(true)
    }

    /// Read next contig with parsed sample/contig names (for multi-sample FASTAs)
    /// Returns: (full_header, sample_name, contig_name, sequence)
    ///
    /// Parses headers in format: >sample#haplotype#chromosome
    /// Example: >S288C#1#chrI -> sample="S288C#1", contig="chrI"
    pub fn read_contig_with_sample(
        &mut self,
    ) -> io::Result<Option<(String, String, String, Contig)>> {
        let (full_header, sequence) = match self.read_contig_impl(true)? {
            Some(result) => result,
            None => return Ok(None),
        };

        // Parse sample name from header
        let (sample_name, contig_name) = parse_sample_from_header(&full_header);

        Ok(Some((full_header, sample_name, contig_name, sequence)))
    }

    /// Read next contig preserving raw format (including newlines)
    pub fn read_contig_raw(&mut self) -> io::Result<Option<(String, Contig)>> {
        let reader = match &mut self.reader {
            Some(r) => r,
            None => return Ok(None),
        };

        let mut contig = Contig::new();

        // Read ID line (starts with '>')
        // Check if we have a buffered header from previous read
        let header_line = if let Some(buffered) = self.next_header.take() {
            buffered
        } else {
            self.buffer.clear();
            let bytes_read = reader.read_until(b'\n', &mut self.buffer)?;
            if bytes_read == 0 {
                return Ok(None);
            }
            self.buffer.clone()
        };

        // Extract ID (skip '>' and trim whitespace)
        let id_line = String::from_utf8_lossy(&header_line);
        let id = id_line.trim_start_matches('>').trim().to_string();

        // Read sequence data until next '>' or EOF
        loop {
            self.buffer.clear();
            let bytes_read = reader.read_until(b'\n', &mut self.buffer)?;

            if bytes_read == 0 {
                // EOF reached
                break;
            }

            // Check if this is the start of a new contig
            if !self.buffer.is_empty() && self.buffer[0] == b'>' {
                // Save this header for the next read
                self.next_header = Some(self.buffer.clone());
                break;
            }

            // Append sequence data
            contig.extend_from_slice(&self.buffer);
        }

        if id.is_empty() || contig.is_empty() {
            return Ok(None);
        }

        Ok(Some((id, contig)))
    }

    /// Internal implementation of contig reading with optional conversion
    fn read_contig_impl(&mut self, converted: bool) -> io::Result<Option<(String, Contig)>> {
        // First read raw
        let (id, raw_contig) = match self.read_contig_raw()? {
            Some(result) => result,
            None => return Ok(None),
        };

        // Filter and optionally convert
        let mut contig = Contig::with_capacity(raw_contig.len());

        for &c in &raw_contig {
            // Only keep characters > 64 (ASCII 'A' is 65)
            if c > 64 && (c as usize) < CNV_NUM.len() {
                if converted {
                    contig.push(CNV_NUM[c as usize]);
                } else {
                    contig.push(c);
                }
            }
        }

        Ok(Some((id, contig)))
    }
}

impl GenomeIO<File> {
    /// Open a FASTA file for reading (uncompressed files only)
    fn open_raw<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(GenomeIO::new(file))
    }
}

impl GenomeIO<Box<dyn Read>> {
    /// Open a FASTA file for reading (supports .gz files)
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;

        // Check if file is gzipped by extension
        let reader: Box<dyn Read> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
            // Use MultiGzDecoder to handle BGZIP (multi-member gzip) files
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };

        Ok(GenomeIO::new(reader))
    }
}

/// Writer for FASTA files
pub struct GenomeWriter<W> {
    writer: W,
}

impl<W: Write> GenomeWriter<W> {
    /// Create a new genome writer
    pub fn new(writer: W) -> Self {
        GenomeWriter { writer }
    }

    /// Save a contig directly (without line wrapping)
    /// If gzip_level > 0, the header is gzip-compressed but sequence is raw
    pub fn save_contig_directly(
        &mut self,
        id: &str,
        contig: &Contig,
        _gzip_level: u32,
    ) -> io::Result<()> {
        // For now, ignore gzip_level and write directly
        // Full gzip support would require the refresh library or similar
        writeln!(self.writer, ">{id}")?;

        // Write sequence with line wrapping at 80 characters (matching C++ AGC behavior)
        const LINE_WIDTH: usize = 80;
        for chunk in contig.chunks(LINE_WIDTH) {
            self.writer.write_all(chunk)?;
            writeln!(self.writer)?;
        }

        Ok(())
    }
}

impl GenomeWriter<File> {
    /// Open a FASTA file for writing
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(GenomeWriter::new(file))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_simple_fasta() {
        let fasta_data = b">seq1\nACGT\n>seq2\nTGCA\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta_data));

        // Read first contig
        let (id, seq) = reader.read_contig().unwrap().unwrap();
        assert_eq!(id, "seq1");
        assert_eq!(seq, b"ACGT");

        // Read second contig
        let (id, seq) = reader.read_contig().unwrap().unwrap();
        assert_eq!(id, "seq2");
        assert_eq!(seq, b"TGCA");

        // No more contigs
        assert!(reader.read_contig().unwrap().is_none());
    }

    #[test]
    fn test_read_converted_fasta() {
        let fasta_data = b">test\nACGT\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta_data));

        let (id, seq) = reader.read_contig_converted().unwrap().unwrap();
        assert_eq!(id, "test");
        // A=0, C=1, G=2, T=3
        assert_eq!(seq, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_cnv_num_table() {
        // Test uppercase
        assert_eq!(CNV_NUM[b'A' as usize], 0);
        assert_eq!(CNV_NUM[b'C' as usize], 1);
        assert_eq!(CNV_NUM[b'G' as usize], 2);
        assert_eq!(CNV_NUM[b'T' as usize], 3);

        // Test lowercase
        assert_eq!(CNV_NUM[b'a' as usize], 0);
        assert_eq!(CNV_NUM[b'c' as usize], 1);
        assert_eq!(CNV_NUM[b'g' as usize], 2);
        assert_eq!(CNV_NUM[b't' as usize], 3);
    }
}
