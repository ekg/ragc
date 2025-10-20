// AGC CLI - Rust implementation
// Compatible with C++ AGC format

use anyhow::Result;
use clap::{Parser, Subcommand};
use ragc_core::{StreamingCompressor, StreamingCompressorConfig, Decompressor, DecompressorConfig};
use std::io::{self, Write};
use std::path::PathBuf;

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
        #[arg(short = 'k', long, default_value_t = 31)]
        kmer_length: u32,

        /// Segment size for splitting contigs
        #[arg(short = 's', long, default_value_t = 60000)]
        segment_size: u32,

        /// Minimum match length for LZ encoding
        #[arg(short = 'm', long, default_value_t = 20)]
        min_match_len: u32,

        /// ZSTD compression level (1-22, higher = better compression but slower)
        #[arg(short = 'c', long, default_value_t = 17)]
        compression_level: i32,

        /// Verbosity level (0=quiet, 1=normal, 2=verbose)
        #[arg(short = 'v', long, default_value_t = 1)]
        verbosity: u32,

        /// Adaptive mode: find new splitters for samples with poor segmentation
        /// (matches C++ AGC -a flag for pangenome-aware compression)
        #[arg(short = 'a', long)]
        adaptive: bool,
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

/// Extract sample name from file path by stripping all genomic file extensions
/// Examples:
///   scerevisiae8.fa.gz -> scerevisiae8
///   genome.fasta       -> genome
///   data.fa            -> data
fn extract_sample_name(path: &PathBuf) -> String {
    let mut name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    // Strip known genomic file extensions (in order)
    let extensions = [".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna", ".gz"];
    for ext in &extensions {
        if name.ends_with(ext) {
            name = name[..name.len() - ext.len()].to_string();
            break;
        }
    }

    name
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
            compression_level,
            verbosity,
            adaptive,
        } => create_archive(
            output,
            inputs,
            kmer_length,
            segment_size,
            min_match_len,
            compression_level,
            verbosity,
            adaptive,
        )?,

        Commands::Info { archive } => {
            eprintln!("Info command not yet implemented for archive: {archive:?}");
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
    compression_level: i32,
    verbosity: u32,
    adaptive: bool,
) -> Result<()> {
    if verbosity > 0 {
        eprintln!("Creating AGC archive: {output:?}");
        eprintln!("Input files: {} FASTA file(s)", inputs.len());
        eprintln!("Parameters:");
        eprintln!("  k-mer length: {kmer_length}");
        eprintln!("  segment size: {segment_size}");
        eprintln!("  min match length: {min_match_len}");
        eprintln!("  compression level: {compression_level}");
        if adaptive {
            eprintln!("  adaptive mode: enabled (pangenome-aware splitters)");
        }
        eprintln!();
    }

    // Use default config (which auto-detects num_threads) and override CLI params
    let config = StreamingCompressorConfig {
        kmer_length,
        segment_size,
        min_match_len,
        compression_level,
        verbosity,
        adaptive_mode: adaptive,
        ..StreamingCompressorConfig::default()
    };

    let output_str = output
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid output path"))?;

    let mut compressor = StreamingCompressor::new(output_str, config)?;

    // Detect if we have multi-sample FASTAs (sample names in headers)
    // or separate files (sample names from filenames)
    if inputs.len() == 1 {
        // Single input file - check if it's multi-sample
        let input_path = &inputs[0];
        if !input_path.exists() {
            anyhow::bail!("Input file not found: {input_path:?}");
        }

        let is_multi_sample = StreamingCompressor::detect_multi_sample_fasta(input_path)?;

        if is_multi_sample {
            // Multi-sample FASTA: group by sample names in headers
            if verbosity > 0 {
                eprintln!("Detected multi-sample FASTA format (sample#haplotype#chromosome)");
                eprintln!("Will group contigs by sample names extracted from headers");
                eprintln!();
            }
            compressor.add_multi_sample_fasta_with_splitters(input_path)?;
        } else {
            // Single-sample file: use filename as sample name, with splitter-based segmentation
            if verbosity > 0 {
                eprintln!("Processing as single-sample file (sample name from filename)");
                eprintln!("Using splitter-based segmentation (matching C++ AGC behavior)");
                eprintln!();
            }
            let sample_name = extract_sample_name(input_path);
            compressor.add_fasta_files_with_splitters(&[(sample_name, input_path.as_path())])?;
        }
    } else {
        // Multiple input files: use filenames as sample names
        if verbosity > 0 {
            eprintln!("Processing multiple files (sample names from filenames)");
            eprintln!();
        }

        let mut fasta_files = Vec::new();
        for input_path in &inputs {
            if !input_path.exists() {
                eprintln!("Warning: Input file not found: {input_path:?}");
                continue;
            }

            let sample_name = extract_sample_name(input_path);
            fasta_files.push((sample_name, input_path.as_path()));
        }

        compressor.add_fasta_files_with_splitters(&fasta_files)?;
    }

    // Finalize the archive
    if verbosity > 0 {
        eprintln!();
        eprintln!("Finalizing archive...");
    }

    compressor.finalize()?;

    if verbosity > 0 {
        eprintln!("Archive created successfully: {output:?}");
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
                eprintln!("Extracting sample: {sample_name}");
            }
            decompressor.write_sample_fasta(sample_name, &output_path)?;
        }
    } else {
        // Extract to temp file then write to stdout
        let temp_path =
            std::env::temp_dir().join(format!("agc_extract_{}.fasta", std::process::id()));
        for sample_name in &samples {
            if verbosity > 0 {
                eprintln!("Extracting sample: {sample_name}");
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
            writeln!(file, "{sample}")?;
        }
    } else {
        for sample in samples {
            println!("{sample}");
        }
    }

    Ok(())
}

fn listctg_command(archive: PathBuf, samples: Vec<String>, output: Option<PathBuf>) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity: 0 };
    let mut decompressor = Decompressor::open(archive_str, config)?;

    let mut output_lines = Vec::new();

    for sample_name in &samples {
        let contigs = decompressor.list_contigs(sample_name)?;
        for contig_name in contigs {
            output_lines.push(format!("{sample_name}\t{contig_name}"));
        }
    }

    if let Some(output_path) = output {
        let mut file = std::fs::File::create(output_path)?;
        for line in output_lines {
            writeln!(file, "{line}")?;
        }
    } else {
        for line in output_lines {
            println!("{line}");
        }
    }

    decompressor.close()?;
    Ok(())
}
