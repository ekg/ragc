#![allow(clippy::all)]
// Integration test for FASTA parsing
// Should produce output identical to C++ test_fasta_parser

use ragc_core::GenomeIO;
use std::io::Cursor;

fn print_bytes(data: &[u8], max_len: usize) {
    for (i, &byte) in data.iter().enumerate().take(max_len) {
        if i > 0 {
            print!(" ");
        }
        print!("{byte:02x}");
    }
}

fn main() {
    // Test 1: Simple FASTA with single contig
    println!("# Test 1: Simple FASTA single contig");
    {
        let fasta = b">seq1\nACGT\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta));

        if let Some((id, contig)) = reader.read_contig().unwrap() {
            println!("ID: {id}");
            println!("Length: {}", contig.len());
            print!("Sequence: ");
            for &c in &contig {
                print!("{}", c as char);
            }
            println!();
            print!("Bytes: ");
            print_bytes(&contig, 100);
            println!();
        }
    }

    // Test 2: Multiple contigs
    println!("\n# Test 2: Multiple contigs");
    {
        let fasta = b">seq1\nACGT\n>seq2\nTGCA\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta));

        let mut count = 0;
        while let Some((id, contig)) = reader.read_contig().unwrap() {
            count += 1;
            println!("Contig {count}:");
            println!("  ID: {id}");
            println!("  Length: {}", contig.len());
            print!("  Sequence: ");
            for &c in &contig {
                print!("{}", c as char);
            }
            println!();
        }
    }

    // Test 3: Converted mode (nucleotide to numeric)
    println!("\n# Test 3: Converted mode");
    {
        let fasta = b">test\nACGT\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta));

        if let Some((id, contig)) = reader.read_contig_converted().unwrap() {
            println!("ID: {id}");
            println!("Length: {}", contig.len());
            print!("Numeric: ");
            print_bytes(&contig, 100);
            println!();
        }
    }

    // Test 4: Multi-line sequence
    println!("\n# Test 4: Multi-line sequence");
    {
        let fasta = b">multiline\nACGT\nTGCA\nAAAA\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta));

        if let Some((id, contig)) = reader.read_contig().unwrap() {
            println!("ID: {id}");
            println!("Length: {}", contig.len());
            print!("Sequence: ");
            for &c in &contig {
                print!("{}", c as char);
            }
            println!();
        }
    }

    // Test 5: IUPAC codes (extended nucleotides)
    println!("\n# Test 5: IUPAC codes");
    {
        let fasta = b">iupac\nACGTNRYSWKMBDHVU\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta));

        if let Some((id, contig)) = reader.read_contig_converted().unwrap() {
            println!("ID: {id}");
            println!("Length: {}", contig.len());
            print!("Numeric: ");
            print_bytes(&contig, 100);
            println!();
        }
    }

    // Test 6: Lowercase nucleotides
    println!("\n# Test 6: Lowercase nucleotides");
    {
        let fasta = b">lowercase\nacgt\n";
        let mut reader = GenomeIO::new(Cursor::new(fasta));

        if let Some((id, contig)) = reader.read_contig_converted().unwrap() {
            println!("ID: {id}");
            println!("Length: {}", contig.len());
            print!("Numeric: ");
            print_bytes(&contig, 100);
            println!();
        }
    }
}
