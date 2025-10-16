// AGC CLI - Rust implementation
// Compatible with C++ AGC format

use ragc_core::{Compressor, CompressorConfig, Decompressor, DecompressorConfig};
use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};
use std::io::{self, Write};

#[derive(Parser, Debug)]
#[command(name = "agc")]
#[command(version, about = "Assembled Genomes Compressor", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Create a new AGC archive from FASTA files
    Create {
        /// Output archive file path
        #[arg(short = 'o', long)]
        output: PathBuf,

        /// Input FASTA files (can specify multiple)
        #[arg(required = true)]
        inputs: Vec<PathBuf>,

        /// K-mer length for splitter identification
        #[arg(short = 'k', long, default_value_t = 21)]
        kmer_length: u32,

        /// Segment size for splitting contigs
        #[arg(short = 's', long, default_value_t = 1000)]
        segment_size: u32,

        /// Minimum match length for LZ encoding
        #[arg(short = 'm', long, default_value_t = 15)]
        min_match_len: u32,

        /// Verbosity level (0=quiet, 1=normal, 2=verbose)
        #[arg(short = 'v', long, default_value_t = 1)]
        verbosity: u32,
    },

    /// Display information about an AGC archive
    Info {
        /// Input archive file path
        archive: PathBuf,
    },

    /// Extract sample(s) from archive
    Getset {
        /// Input archive file path
        archive: PathBuf,

        /// Sample name(s) to extract
        #[arg(required = true)]
        samples: Vec<String>,

        /// Output file (default: stdout)
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,

        /// Verbosity level (0=quiet, 1=normal, 2=verbose)
        #[arg(short = 'v', long, default_value_t = 0)]
        verbosity: u32,
    },

    /// List sample names in archive
    Listset {
        /// Input archive file path
        archive: PathBuf,

        /// Output file (default: stdout)
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },

    /// List contig names for sample(s) in archive
    Listctg {
        /// Input archive file path
        archive: PathBuf,

        /// Sample name(s) to list contigs for
        #[arg(required = true)]
        samples: Vec<String>,

        /// Output file (default: stdout)
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Create {
            output,
            inputs,
            kmer_length,
            segment_size,
            min_match_len,
            verbosity,
        } => create_archive(output, inputs, kmer_length, segment_size, min_match_len, verbosity)?,

        Commands::Info { archive } => {
            eprintln!("Info command not yet implemented for archive: {:?}", archive);
            eprintln!("This will be implemented in a future version.");
        }

        Commands::Getset {
            archive,
            samples,
            output,
            verbosity,
        } => getset_command(archive, samples, output, verbosity)?,

        Commands::Listset { archive, output } => listset_command(archive, output)?,

        Commands::Listctg {
            archive,
            samples,
            output,
        } => listctg_command(archive, samples, output)?,
    }

    Ok(())
}

fn create_archive(
    output: PathBuf,
    inputs: Vec<PathBuf>,
    kmer_length: u32,
    segment_size: u32,
    min_match_len: u32,
    verbosity: u32,
) -> Result<()> {
    if verbosity > 0 {
        eprintln!("Creating AGC archive: {:?}", output);
        eprintln!("Input files: {} FASTA file(s)", inputs.len());
        eprintln!("Parameters:");
        eprintln!("  k-mer length: {}", kmer_length);
        eprintln!("  segment size: {}", segment_size);
        eprintln!("  min match length: {}", min_match_len);
    }

    let config = CompressorConfig {
        kmer_length,
        segment_size,
        min_match_len,
        verbosity,
    };

    let output_str = output
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid output path"))?;

    let mut compressor = Compressor::new(output_str, config)?;

    // Process each input file
    for input_path in &inputs {
        if !input_path.exists() {
            eprintln!("Warning: Input file not found: {:?}", input_path);
            continue;
        }

        // Extract sample name from file path (without extension)
        let sample_name = input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        if verbosity > 0 {
            eprintln!("Processing: {} <- {:?}", sample_name, input_path);
        }

        compressor.add_fasta_file(sample_name, input_path)?;
    }

    // Finalize the archive
    if verbosity > 0 {
        eprintln!("Finalizing archive...");
    }

    compressor.finalize()?;

    if verbosity > 0 {
        eprintln!("Archive created successfully: {:?}", output);
    }

    Ok(())
}

fn getset_command(
    archive: PathBuf,
    samples: Vec<String>,
    output: Option<PathBuf>,
    verbosity: u32,
) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity };
    let mut decompressor = Decompressor::open(archive_str, config)?;

    // If output file specified, extract to file
    // Otherwise, extract to stdout (via temp file for simplicity)
    if let Some(output_path) = output {
        // Extract each sample to the output file (append mode)
        for sample_name in &samples {
            if verbosity > 0 {
                eprintln!("Extracting sample: {}", sample_name);
            }
            decompressor.write_sample_fasta(sample_name, &output_path)?;
        }
    } else {
        // Extract to temp file then write to stdout
        let temp_path = std::env::temp_dir().join(format!("agc_extract_{}.fasta", std::process::id()));
        for sample_name in &samples {
            if verbosity > 0 {
                eprintln!("Extracting sample: {}", sample_name);
            }
            decompressor.write_sample_fasta(sample_name, &temp_path)?;
        }
        // Write temp file to stdout
        let contents = std::fs::read(&temp_path)?;
        io::stdout().write_all(&contents)?;
        std::fs::remove_file(&temp_path)?;
    }

    decompressor.close()?;
    Ok(())
}

fn listset_command(archive: PathBuf, output: Option<PathBuf>) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity: 0 };
    let decompressor = Decompressor::open(archive_str, config)?;

    let samples = decompressor.list_samples();

    if let Some(output_path) = output {
        let mut file = std::fs::File::create(output_path)?;
        for sample in samples {
            writeln!(file, "{}", sample)?;
        }
    } else {
        for sample in samples {
            println!("{}", sample);
        }
    }

    Ok(())
}

fn listctg_command(
    archive: PathBuf,
    samples: Vec<String>,
    output: Option<PathBuf>,
) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity: 0 };
    let mut decompressor = Decompressor::open(archive_str, config)?;

    let mut output_lines = Vec::new();

    for sample_name in &samples {
        let contigs = decompressor.list_contigs(sample_name)?;
        for contig_name in contigs {
            output_lines.push(format!("{}\t{}", sample_name, contig_name));
        }
    }

    if let Some(output_path) = output {
        let mut file = std::fs::File::create(output_path)?;
        for line in output_lines {
            writeln!(file, "{}", line)?;
        }
    } else {
        for line in output_lines {
            println!("{}", line);
        }
    }

    decompressor.close()?;
    Ok(())
}
