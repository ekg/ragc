// AGC CLI - Rust implementation
// Compatible with C++ AGC format

// Use jemalloc if enabled (reduces memory overhead and fragmentation)
#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod inspect;

use anyhow::Result;
use clap::{Parser, Subcommand};
use ragc_core::{
    contig_iterator::ContigIterator, Decompressor, DecompressorConfig, MultiFileIterator,
    StreamingQueueCompressor, StreamingQueueConfig,
};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

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

        /// Concatenated genomes mode: treat all input as one continuous sample
        /// Disables split detection. Useful for pangenomes viewed as concatenated sequences.
        #[arg(long)]
        concatenated: bool,

        /// Number of threads for parallel compression (default: auto-detect)
        /// Set to 1 for single-threaded compression
        #[arg(short = 't', long)]
        threads: Option<usize>,

        /// Use legacy batch mode instead of streaming queue (for debugging/comparison)
        /// By default, uses streaming queue mode for constant memory usage.
        #[arg(long)]
        batch: bool,

        /// Queue capacity in bytes for streaming mode (default: 2GB)
        /// Only applies to streaming queue mode (not --batch).
        /// Accepts suffixes: K, M, G (e.g., "512M", "2G")
        #[arg(long, default_value = "2G")]
        queue_capacity: String,
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

        /// Sample name(s) to extract (not required if --prefix is used)
        samples: Vec<String>,

        /// Extract all samples matching this prefix
        #[arg(short = 'p', long)]
        prefix: Option<String>,

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

    /// Inspect archive structure (groups, segments, compression details)
    Inspect {
        /// Input archive file path
        archive: PathBuf,

        /// Show detailed segment information
        #[arg(short = 's', long)]
        segments: bool,

        /// Filter by specific group ID
        #[arg(short = 'g', long)]
        group_id: Option<u32>,

        /// Show only single-segment groups (for debugging fragmentation)
        #[arg(long)]
        single_groups: bool,

        /// Show segment layout in CSV format for comparison
        #[arg(long)]
        segment_layout: bool,

        /// Look up segment by sample name
        #[arg(long)]
        sample: Option<String>,

        /// Look up segment by contig name (requires --sample)
        #[arg(long)]
        contig: Option<String>,

        /// Look up segment by index (requires --sample and --contig)
        #[arg(long)]
        index: Option<usize>,

        /// Verbosity level (0=quiet, 1=normal, 2=verbose)
        #[arg(short = 'v', long, default_value_t = 0)]
        verbosity: u32,
    },
}

/// Parse capacity string with suffixes (K, M, G) into bytes
/// Examples: "512M" -> 536870912, "2G" -> 2147483648
fn parse_capacity(s: &str) -> Result<usize> {
    let s = s.trim().to_uppercase();

    if let Some(num_str) = s.strip_suffix('K') {
        let num: usize = num_str.parse()?;
        Ok(num * 1024)
    } else if let Some(num_str) = s.strip_suffix('M') {
        let num: usize = num_str.parse()?;
        Ok(num * 1024 * 1024)
    } else if let Some(num_str) = s.strip_suffix('G') {
        let num: usize = num_str.parse()?;
        Ok(num * 1024 * 1024 * 1024)
    } else {
        // No suffix, parse as bytes
        Ok(s.parse()?)
    }
}

/// Extract sample name from file path by stripping all genomic file extensions
/// Examples:
///   scerevisiae8.fa.gz -> scerevisiae8
///   genome.fasta       -> genome
///   data.fa            -> data
fn extract_sample_name(path: &Path) -> String {
    let mut name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    // Strip known genomic file extensions (in order)
    let extensions = [
        ".fa.gz",
        ".fasta.gz",
        ".fna.gz",
        ".fa",
        ".fasta",
        ".fna",
        ".gz",
    ];
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
            concatenated,
            threads,
            batch,
            queue_capacity,
        } => create_archive(
            output,
            inputs,
            kmer_length,
            segment_size,
            min_match_len,
            compression_level,
            verbosity,
            adaptive,
            concatenated,
            threads,
            batch,
            &queue_capacity,
        )?,

        Commands::Info { archive } => {
            eprintln!("Info command not yet implemented for archive: {archive:?}");
            eprintln!("This will be implemented in a future version.");
        }

        Commands::Getset {
            archive,
            samples,
            prefix,
            output,
            verbosity,
        } => getset_command(archive, samples, prefix, output, verbosity)?,

        Commands::Listset { archive, output } => listset_command(archive, output)?,

        Commands::Listctg {
            archive,
            samples,
            output,
        } => listctg_command(archive, samples, output)?,

        Commands::Inspect {
            archive,
            segments,
            group_id,
            single_groups,
            segment_layout,
            sample,
            contig,
            index,
            verbosity,
        } => {
            let config = inspect::InspectConfig {
                verbosity,
                show_groups: true,
                show_segments: segments,
                group_id_filter: group_id,
                sample_filter: sample,
                contig_filter: contig,
                segment_index: index,
                show_single_segment_groups: single_groups,
                show_segment_layout: segment_layout,
            };
            inspect::inspect_archive(archive, config)?
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn create_archive(
    output: PathBuf,
    inputs: Vec<PathBuf>,
    kmer_length: u32,
    segment_size: u32,
    min_match_len: u32,
    compression_level: i32,
    verbosity: u32,
    adaptive: bool,
    concatenated: bool,
    threads: Option<usize>,
    batch: bool,
    queue_capacity_str: &str,
) -> Result<()> {
    // Determine thread count (use provided or auto-detect)
    let num_threads = threads.unwrap_or_else(|| {
        let num_cpus = num_cpus::get();
        if num_cpus < 8 {
            num_cpus
        } else {
            num_cpus - 1
        }
    });

    if verbosity > 0 {
        eprintln!("Creating AGC archive: {output:?}");
        eprintln!("Input files: {} FASTA file(s)", inputs.len());
        eprintln!("Parameters:");
        eprintln!("  k-mer length: {kmer_length}");
        eprintln!("  segment size: {segment_size}");
        eprintln!("  min match length: {min_match_len}");
        eprintln!("  compression level: {compression_level}");
        eprintln!("  threads: {num_threads}");
        if !batch {
            let capacity = parse_capacity(queue_capacity_str)?;
            eprintln!("  mode: streaming queue (constant memory, default)");
            eprintln!(
                "  queue capacity: {} bytes ({:.2} GB)",
                capacity,
                capacity as f64 / 1024.0 / 1024.0 / 1024.0
            );
        } else {
            eprintln!("  mode: batch (legacy)");
        }
        if adaptive {
            eprintln!("  adaptive mode: enabled (pangenome-aware splitters)");
        }
        if concatenated {
            eprintln!("  concatenated mode: enabled (all contigs as one sample)");
        }
        eprintln!();
    }

    // Branch: Streaming queue mode (default) vs. batch mode (legacy)
    if !batch {
        // Streaming queue mode: constant memory with bounded queue (DEFAULT)
        if adaptive || concatenated {
            anyhow::bail!(
                "Streaming queue mode does not support --adaptive or --concatenated flags yet"
            );
        }

        let output_str = output
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid output path"))?;

        let queue_capacity = parse_capacity(queue_capacity_str)?;
        let config = StreamingQueueConfig {
            k: kmer_length as usize,
            segment_size: segment_size as usize,
            min_match_len: min_match_len as usize,
            queue_capacity,
            num_threads,
            verbosity: verbosity as usize,
            adaptive_mode: adaptive,
            ..StreamingQueueConfig::default()
        };

        // Detect splitters from first input file
        if inputs.is_empty() {
            anyhow::bail!("No input files provided");
        }

        if verbosity > 0 {
            eprintln!("Detecting splitters from first input file: {:?}", inputs[0]);
        }

        let splitters = ragc_core::determine_splitters_streaming(
            &inputs[0],
            kmer_length as usize,
            segment_size as usize,
        )?
        .0;

        if verbosity > 0 {
            eprintln!("Found {} splitters", splitters.len());
            eprintln!();
        }

        // Create compressor with splitters from first file
        let mut compressor =
            StreamingQueueCompressor::with_splitters(output_str, config, splitters)?;

        // CRITICAL: Process reference sample (first file) completely BEFORE other samples
        // This matches C++ AGC architecture and ensures splits work correctly.
        // C++ AGC loads the entire reference sample first, creating all its segment groups
        // and terminators, then processes remaining samples. This allows splits to find
        // target groups that already exist in map_segments.

        if verbosity > 0 {
            eprintln!("Processing reference sample (first input): {:?}", inputs[0]);
        }

        // Process ONLY the first file (reference sample)
        let mut ref_iterator = MultiFileIterator::new(vec![inputs[0].clone()])?;
        while let Some((sample_name, contig_name, sequence)) = ref_iterator.next_contig()? {
            if !sequence.is_empty() {
                compressor.push(sample_name, contig_name, sequence)?;
            }
        }

        // Wait for reference sample to be completely processed
        // This ensures all reference groups and terminators are created before
        // we start processing other samples
        if verbosity > 0 {
            eprintln!("Reference sample queued - waiting for processing to complete...");
        }
        compressor.drain()?;

        if verbosity > 0 {
            eprintln!("Reference sample complete! Processing remaining {} samples...", inputs.len() - 1);
            eprintln!();
        }

        // Now process remaining samples
        // At this point, all reference groups exist, so splits can happen
        if inputs.len() > 1 {
            let mut remaining_iterator = MultiFileIterator::new(inputs[1..].to_vec())?;
            while let Some((sample_name, contig_name, sequence)) = remaining_iterator.next_contig()? {
                if !sequence.is_empty() {
                    compressor.push(sample_name, contig_name, sequence)?;
                }
            }
        }

        compressor.finalize()?;

        if verbosity > 0 {
            eprintln!("\nArchive created successfully: {output:?}");
        }

        return Ok(());
    }

    Ok(())
}

fn getset_command(
    archive: PathBuf,
    samples: Vec<String>,
    prefix: Option<String>,
    output: Option<PathBuf>,
    verbosity: u32,
) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity };
    let mut decompressor = Decompressor::open(archive_str, config)?;

    // Determine which samples to extract
    let samples_to_extract = if let Some(prefix_str) = prefix {
        // Extract by prefix
        if verbosity > 0 {
            eprintln!("Finding samples with prefix: {prefix_str}");
        }
        let matching_samples = decompressor.list_samples_with_prefix(&prefix_str);
        if matching_samples.is_empty() {
            anyhow::bail!("No samples found matching prefix '{prefix_str}'");
        }
        if verbosity > 0 {
            eprintln!("Found {} matching samples", matching_samples.len());
        }
        matching_samples
    } else if !samples.is_empty() {
        // Extract specific samples
        samples
    } else {
        anyhow::bail!("Must specify either sample names or --prefix");
    };

    // If output file specified, extract to file
    // Otherwise, extract to stdout (via temp file for simplicity)
    if let Some(output_path) = output {
        // Extract each sample to the output file (append mode)
        for sample_name in &samples_to_extract {
            if verbosity > 0 {
                eprintln!("Extracting sample: {sample_name}");
            }
            decompressor.write_sample_fasta(sample_name, &output_path)?;
        }
    } else {
        // Extract to temp file then write to stdout
        let temp_path =
            std::env::temp_dir().join(format!("agc_extract_{}.fasta", std::process::id()));
        for sample_name in &samples_to_extract {
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
