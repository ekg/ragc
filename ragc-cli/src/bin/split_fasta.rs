use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn parse_sample_name(header: &str) -> Option<String> {
    let parts: Vec<&str> = header.split('#').collect();
    if parts.len() >= 2 {
        Some(format!("{}#{}", parts[0], parts[1]))
    } else {
        None
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.fa[.gz]> <output_dir>", args[0]);
        eprintln!("Splits a PanSN-formatted FASTA file into one file per sample#haplotype");
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_dir = &args[2];

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Open input file (with optional gzip decompression)
    let input_file = File::open(input_path)?;
    let reader: Box<dyn BufRead> = if input_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(input_file)))
    } else {
        Box::new(BufReader::new(input_file))
    };

    let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();
    let mut current_sample: Option<String> = None;
    let mut header_count = 0;

    for line_result in reader.lines() {
        let line = line_result?;

        if line.starts_with('>') {
            header_count += 1;
            if header_count <= 10 {
                eprintln!("DEBUG: Header {header_count}: {line}");
            }

            // Parse sample name from header
            let header = line.trim_start_matches('>');
            if let Some(sample_name) = parse_sample_name(header) {
                if header_count <= 10 {
                    eprintln!("DEBUG: Parsed sample: {sample_name}");
                }
                current_sample = Some(sample_name.clone());

                // Create writer for this sample if it doesn't exist
                if !writers.contains_key(&sample_name) {
                    let safe_name = sample_name.replace(['/', '#'], "_");
                    let output_path = Path::new(output_dir).join(format!("{safe_name}.fa"));
                    let file = File::create(&output_path)?;
                    eprintln!("Creating file: {}", output_path.display());
                    writers.insert(sample_name.clone(), BufWriter::new(file));
                }
            } else {
                eprintln!("Warning: Cannot parse sample from header: {line}");
                current_sample = None;
            }
        }

        // Write line to current sample's file
        if let Some(ref sample) = current_sample {
            if let Some(writer) = writers.get_mut(sample) {
                writeln!(writer, "{line}")?;
            }
        }
    }

    // Flush all writers
    for (sample, mut writer) in writers {
        writer.flush()?;
        eprintln!("Wrote sample: {sample}");
    }

    eprintln!("Splitting complete!");
    Ok(())
}
