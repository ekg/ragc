#![allow(clippy::too_many_arguments)]
#![allow(clippy::redundant_field_names)]
#![allow(clippy::single_match)]
#![allow(clippy::redundant_locals)]
#![allow(unexpected_cfgs)]

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

        /// Pack cardinality: number of contigs between batch flushes
        /// (matches C++ AGC -l flag for max_no_contigs_before_synchronization)
        #[arg(short = 'l', long, default_value_t = 50)]
        pack_cardinality: u32,

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

        /// Fallback minimizers fraction for grouping terminator segments (0.0-1.0)
        /// When > 0, enables fallback k-mer indexing for segments with missing terminators.
        /// Higher values (e.g., 0.1) index more k-mers for better grouping at memory cost.
        /// Default 0.0 disables fallback minimizers (matches C++ AGC default).
        #[arg(long, default_value_t = 0.0)]
        fallback_frac: f64,

        /// Use C++ AGC compression via FFI for byte-identical archives
        /// This delegates the entire compression to the original C++ AGC implementation.
        /// Requires: cargo build --features cpp_agc
        #[cfg_attr(not(feature = "cpp_agc"), arg(hide = true))]
        #[arg(long, default_value_t = false)]
        cpp_agc: bool,
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

    /// Extract a subsequence from a contig (efficient range query)
    Getrange {
        /// Input archive file path
        archive: PathBuf,

        /// Sample name containing the contig
        #[arg(short = 's', long)]
        sample: String,

        /// Contig name
        #[arg(short = 'c', long)]
        contig: String,

        /// Start position (0-based, inclusive)
        #[arg(long)]
        start: usize,

        /// End position (0-based, exclusive). If omitted, extracts to end of contig.
        #[arg(long)]
        end: Option<usize>,

        /// Output file (default: stdout)
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,

        /// Output format: fasta or raw (default: fasta)
        #[arg(short = 'f', long, default_value = "fasta")]
        format: String,

        /// Verbosity level (0=quiet, 1=normal, 2=verbose)
        #[arg(short = 'v', long, default_value_t = 0)]
        verbosity: u32,
    },

    /// Get the length of a contig without extracting it
    Ctglen {
        /// Input archive file path
        archive: PathBuf,

        /// Sample name containing the contig
        #[arg(short = 's', long)]
        sample: String,

        /// Contig name
        #[arg(short = 'c', long)]
        contig: String,
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

        /// Compare with another archive (reports differences in grouping, segments, references, packs)
        #[arg(long)]
        compare: Option<PathBuf>,

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

        /// Show pack layout (how segments are organized into ZSTD packs)
        #[arg(long)]
        pack_layout: bool,

        /// Show compression statistics for each stream (group packs)
        #[arg(short = 'c', long)]
        compression: bool,

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

    /// Debug: extract a segment and two group references for cost verification
    DebugCost {
        /// Input archive file path
        archive: PathBuf,
        /// Sample name
        #[arg(long)]
        sample: String,
        /// Contig name
        #[arg(long)]
        contig: String,
        /// Segment index within contig
        #[arg(long)]
        index: usize,
        /// Left group_id (k1, k_mid). If omitted, uses prev segment's group_id.
        #[arg(long)]
        left_group: Option<u32>,
        /// Right group_id (k_mid, k2). If omitted, uses current segment's group_id.
        #[arg(long)]
        right_group: Option<u32>,
        /// Left costs placed at prefix (true) or suffix (false)
        #[arg(long, default_value_t = true)]
        left_prefix: bool,
        /// Right costs placed at prefix (true) or suffix (false)
        #[arg(long, default_value_t = false)]
        right_prefix: bool,
        /// Verbosity
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
#[cfg(feature = "cpp_agc")]
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
            pack_cardinality,
            compression_level,
            verbosity,
            adaptive,
            concatenated,
            threads,
            batch,
            queue_capacity,
            fallback_frac,
            cpp_agc,
        } => create_archive(
            output,
            inputs,
            kmer_length,
            segment_size,
            min_match_len,
            pack_cardinality,
            compression_level,
            verbosity,
            adaptive,
            concatenated,
            threads,
            batch,
            &queue_capacity,
            fallback_frac,
            cpp_agc,
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

        Commands::Getrange {
            archive,
            sample,
            contig,
            start,
            end,
            output,
            format,
            verbosity,
        } => getrange_command(archive, sample, contig, start, end, output, format, verbosity)?,

        Commands::Ctglen {
            archive,
            sample,
            contig,
        } => ctglen_command(archive, sample, contig)?,

        Commands::Listset { archive, output } => listset_command(archive, output)?,

        Commands::Listctg {
            archive,
            samples,
            output,
        } => listctg_command(archive, samples, output)?,

        Commands::Inspect {
            archive,
            compare,
            segments,
            group_id,
            single_groups,
            segment_layout,
            pack_layout,
            compression,
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
                show_pack_layout: pack_layout,
                show_compression: compression,
                compare_with: compare,
            };
            inspect::inspect_archive(archive, config)?
        }

        Commands::DebugCost {
            archive,
            sample,
            contig,
            index,
            left_group,
            right_group,
            left_prefix,
            right_prefix,
            verbosity,
        } => {
            debug_cost_command(
                archive,
                sample,
                contig,
                index,
                left_group,
                right_group,
                left_prefix,
                right_prefix,
                verbosity,
            )?;
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
    pack_cardinality: u32,
    compression_level: i32,
    verbosity: u32,
    adaptive: bool,
    concatenated: bool,
    threads: Option<usize>,
    batch: bool,
    queue_capacity_str: &str,
    fallback_frac: f64,
    cpp_agc: bool,
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

    // Initialize Rayon's global thread pool with the user's thread count
    // This MUST happen BEFORE any par_iter() is called (e.g., in splitter discovery)
    // Without this, Rayon initializes with all CPUs on first par_iter(), ignoring -t flag
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok(); // Ignore error if already initialized

    // DO NOT SORT: C++ AGC processes files in command-line order
    // The first file in the input list becomes the reference for compression
    // Sorting alphabetically breaks compatibility with C++ AGC
    let inputs = inputs;

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
        #[cfg(feature = "cpp_agc")]
        if cpp_agc {
            eprintln!("  C++ AGC FFI mode: enabled (byte-identical archives)");
        }
        eprintln!();
    }

    // C++ AGC FFI mode: delegate entire compression to C++ AGC for byte-identical archives
    #[cfg(feature = "cpp_agc")]
    if cpp_agc {
        if verbosity > 0 {
            eprintln!("Using C++ AGC compression via FFI...");
        }

        // Build sample_files list: (sample_name, file_path)
        let sample_files: Vec<(String, String)> = inputs
            .iter()
            .map(|p| {
                let sample_name = extract_sample_name(p);
                let file_path = p.to_string_lossy().to_string();
                (sample_name, file_path)
            })
            .collect();

        ragc_core::agc_compress_ffi::compress_with_cpp_agc(
            &output,
            &sample_files,
            kmer_length,
            segment_size,
            min_match_len,
            pack_cardinality,
            concatenated,
            adaptive,
            verbosity,
            num_threads as u32,
            0.0, // fallback_frac
        )?;

        if verbosity > 0 {
            eprintln!("\nArchive created successfully: {output:?}");
        }
        return Ok(());
    }

    // Error if cpp_agc requested but feature not enabled
    #[cfg(not(feature = "cpp_agc"))]
    if cpp_agc {
        anyhow::bail!("--cpp-agc requires building with: cargo build --features cpp_agc");
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

        // FIX 5: Determine if we're in concatenated genomes mode
        // Matches C++ AGC behavior: concatenated mode uses pack_cardinality for sync tokens,
        // non-concatenated mode uses sample boundaries for sync tokens
        let concatenated_genomes = concatenated || inputs.len() == 1;

        let config = StreamingQueueConfig {
            k: kmer_length as usize,
            segment_size: segment_size as usize,
            min_match_len: min_match_len as usize,
            pack_size: pack_cardinality as usize,
            queue_capacity,
            num_threads,
            verbosity: verbosity as usize,
            adaptive_mode: adaptive,
            fallback_frac,
            concatenated_genomes,
            ..StreamingQueueConfig::default()
        };

        // Splitter discovery strategy depends on adaptive mode (matching C++ AGC)
        if inputs.is_empty() {
            anyhow::bail!("No input files provided");
        }

        let phase_start = std::time::Instant::now();
        let splitters = if adaptive {
            // Adaptive mode: Discover splitters from hard contigs across ALL files
            // (matches C++ AGC with -a flag)
            ragc_core::two_pass_splitter_discovery(
                &inputs,
                kmer_length as usize,
                segment_size as usize,
                verbosity as usize,
            )?
        } else if inputs.len() == 1 {
            // Single-file PanSN mode: Use splitters from FIRST SAMPLE only
            // This matches multi-file behavior where first file = reference sample
            let (splitters, _singletons, _duplicates) =
                ragc_core::determine_splitters_streaming_first_sample(
                    &inputs[0],
                    kmer_length as usize,
                    segment_size as usize,
                )?;
            splitters
        } else {
            // Multi-file mode: Use splitters from first file (reference sample)
            // (matches C++ AGC default behavior - agc_compressor.cpp:428-563)
            let (splitters, _singletons, _duplicates) = ragc_core::determine_splitters_streaming(
                &inputs[0],
                kmer_length as usize,
                segment_size as usize,
            )?;
            splitters
        };

        if verbosity > 0 {
            eprintln!("Final splitter count: {}", splitters.len());
            eprintln!(
                "CLI_TIMING: Splitter discovery took {:?}",
                phase_start.elapsed()
            );
            eprintln!();
        }

        let phase_start = std::time::Instant::now();
        // Create compressor with combined splitter set
        let mut compressor =
            StreamingQueueCompressor::with_splitters(output_str, config, splitters)?;

        // CRITICAL: Process reference sample completely BEFORE other samples
        // This matches C++ AGC architecture and ensures splits work correctly.
        // C++ AGC loads the entire reference sample first, creating all its segment groups
        // and terminators, then processes remaining samples. This allows splits to find
        // target groups that already exist in map_segments.

        if inputs.len() == 1 {
            // SINGLE-FILE MODE: Stream through file, treating sample name changes as boundaries
            // This handles PanSN files with multiple samples (e.g., >sample1#0#chr1, >sample2#0#chr1)
            // Requires samples to be sorted (all of sample1, then all of sample2, etc.)
            if verbosity > 0 {
                eprintln!("Single-file PanSN mode: detecting sample boundaries from headers");
                eprintln!("Processing input: {:?}", inputs[0]);
            }

            let mut file_iterator = MultiFileIterator::new(vec![inputs[0].clone()])?;
            let mut current_sample: Option<String> = None;
            let mut seen_samples: std::collections::HashSet<String> =
                std::collections::HashSet::new();
            let mut sample_count = 0;
            let mut is_reference_done = false;

            while let Some((sample_name, contig_name, sequence)) = file_iterator.next_contig()? {
                if sequence.is_empty() {
                    continue;
                }

                // Check if sample changed
                let sample_changed = current_sample.as_ref() != Some(&sample_name);

                if sample_changed {
                    // Verify samples are sorted (no revisiting)
                    if seen_samples.contains(&sample_name) {
                        anyhow::bail!(
                            "Single-file PanSN mode requires samples to be sorted by name.\n\
                             Saw sample '{}' again after seeing other samples.\n\
                             Please sort your FASTA file by sample name (e.g., using: ragc sort-fasta).",
                            sample_name
                        );
                    }

                    // If we've completed the reference sample, drain before continuing
                    if !is_reference_done && current_sample.is_some() {
                        if verbosity > 0 {
                            eprintln!(
                                "Reference sample queued - waiting for processing to complete..."
                            );
                        }
                        compressor.drain()?;
                        is_reference_done = true;
                        if verbosity > 0 {
                            eprintln!("Reference sample complete!");
                        }
                    }

                    // Mark previous sample as seen
                    if let Some(prev) = current_sample.take() {
                        seen_samples.insert(prev);
                    }

                    // Start new sample
                    sample_count += 1;
                    if verbosity > 0 {
                        eprintln!("Processing sample {}: {}", sample_count, sample_name);
                    }
                    current_sample = Some(sample_name.clone());
                }

                compressor.push(sample_name, contig_name, sequence)?;
            }

            if verbosity > 0 {
                eprintln!("Processed {} samples from single file", sample_count);
            }
        } else {
            // MULTI-FILE MODE: Each file is a separate sample
            if verbosity > 0 {
                eprintln!(
                    "CLI_TIMING: Compressor init took {:?}",
                    phase_start.elapsed()
                );
                eprintln!("Processing reference sample (first input): {:?}", inputs[0]);
            }

            let phase_start = std::time::Instant::now();
            // Process ONLY the first file (reference sample)
            let mut ref_iterator = MultiFileIterator::new(vec![inputs[0].clone()])?;
            while let Some((sample_name, contig_name, sequence)) = ref_iterator.next_contig()? {
                if !sequence.is_empty() {
                    compressor.push(sample_name, contig_name, sequence)?;
                }
            }

            // Wait for reference sample to be completely processed
            if verbosity > 0 {
                eprintln!("Reference sample queued - waiting for processing to complete...");
            }
            compressor.drain()?;

            // CRITICAL: Flush reference sample's terminators so they're available for segment splitting
            // in subsequent samples. Without this, segment splitting can't find middle splitters.
            compressor.sync_and_flush("AAA#0_REF")?;

            if verbosity > 0 {
                eprintln!(
                    "CLI_TIMING: Reference sample took {:?}",
                    phase_start.elapsed()
                );
                eprintln!(
                    "Reference sample complete! Processing remaining {} samples...",
                    inputs.len() - 1
                );
                eprintln!();
            }

            let phase_start = std::time::Instant::now();
            // Process remaining samples ONE FILE AT A TIME
            for input_file in &inputs[1..] {
                if verbosity > 0 {
                    eprintln!("Processing sample file: {:?}", input_file);
                }

                let mut file_iterator = MultiFileIterator::new(vec![input_file.clone()])?;
                while let Some((sample_name, contig_name, sequence)) =
                    file_iterator.next_contig()?
                {
                    if !sequence.is_empty() {
                        compressor.push(sample_name, contig_name, sequence)?;
                    }
                }

                if verbosity > 0 {
                    eprintln!("Sample complete!\n");
                }
            }
            if verbosity > 0 {
                eprintln!(
                    "CLI_TIMING: Remaining samples took {:?}",
                    phase_start.elapsed()
                );
            }
        }

        let phase_start = std::time::Instant::now();
        compressor.finalize()?;
        if verbosity > 0 {
            eprintln!("CLI_TIMING: Finalize took {:?}", phase_start.elapsed());
        }

        if verbosity > 0 {
            eprintln!("\nArchive created successfully: {output:?}");
        }

        return Ok(());
    }

    Ok(())
}

fn write_bin<P: AsRef<Path>>(path: P, data: &[u8]) -> Result<()> {
    std::fs::create_dir_all(path.as_ref().parent().unwrap_or(Path::new(".")))?;
    std::fs::write(path, data)?;
    Ok(())
}

fn debug_cost_command(
    archive: PathBuf,
    sample: String,
    contig: String,
    index: usize,
    left_group: Option<u32>,
    right_group: Option<u32>,
    left_prefix: bool,
    right_prefix: bool,
    verbosity: u32,
) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;
    let mut dec = Decompressor::open(archive_str, DecompressorConfig { verbosity })?;

    // Load contig descriptors
    let segments = dec.get_contig_segments_desc(&sample, &contig)?;
    if index >= segments.len() {
        anyhow::bail!(
            "Index {} out of range (contig has {} segments)",
            index,
            segments.len()
        );
    }

    let seg_desc = &segments[index];
    let seg_data = dec.get_segment_data_by_desc(seg_desc)?;

    // Auto-pick neighbor groups if not provided
    let left_gid = if let Some(g) = left_group {
        g
    } else {
        if index == 0 {
            anyhow::bail!("No left group for index 0; provide --left-group explicitly");
        }
        segments[index - 1].group_id
    };
    let right_gid = if let Some(g) = right_group {
        g
    } else {
        segments[index].group_id
    };

    let left_ref = dec.get_reference_segment(left_gid)?;
    let right_ref = dec.get_reference_segment(right_gid)?;

    // Write to tmp folder
    let out_dir = PathBuf::from(format!(
        "./tmp/cost_debug_{}_{}_{}",
        sample.replace('#', "_"),
        contig.replace('#', "_"),
        index
    ));
    std::fs::create_dir_all(&out_dir)?;
    let seg_path = out_dir.join("segment.bin");
    let left_path = out_dir.join("left_ref.bin");
    let right_path = out_dir.join("right_ref.bin");
    write_bin(&seg_path, &seg_data)?;
    write_bin(&left_path, &left_ref)?;
    write_bin(&right_path, &right_ref)?;

    println!("Wrote files:");
    println!(
        "  segment:   {} (len={})",
        seg_path.display(),
        seg_data.len()
    );
    println!(
        "  left_ref:  {} (group_id={} len={})",
        left_path.display(),
        left_gid,
        left_ref.len()
    );
    println!(
        "  right_ref: {} (group_id={} len={})",
        right_path.display(),
        right_gid,
        right_ref.len()
    );

    // Try run cost verifier if built
    let verifier = PathBuf::from("./target/cost_verifier");
    if !verifier.exists() {
        println!(
            "\nCost verifier not found at {}. To build it:",
            verifier.display()
        );
        println!("  scripts/build_cost_verifier.sh");
    } else {
        println!("\nRunning C++ verifier (left, prefix={})…", left_prefix);
        let out = std::process::Command::new(&verifier)
            .args([
                if left_prefix { "1" } else { "0" },
                left_path.to_str().unwrap(),
                seg_path.to_str().unwrap(),
            ])
            .output()?;
        std::io::stdout().write_all(&out.stdout)?;

        println!("\nRunning C++ verifier (right, prefix={})…", right_prefix);
        let out2 = std::process::Command::new(&verifier)
            .args([
                if right_prefix { "1" } else { "0" },
                right_path.to_str().unwrap(),
                seg_path.to_str().unwrap(),
            ])
            .output()?;
        std::io::stdout().write_all(&out2.stdout)?;
    }

    // Compute Rust cost vectors for comparison
    use ragc_core::LZDiff;
    let mut lz_left = LZDiff::new(dec.min_match_len);
    lz_left.prepare(&left_ref);
    let mut lz_right = LZDiff::new(dec.min_match_len);
    lz_right.prepare(&right_ref);

    let v_left_raw = lz_left.get_coding_cost_vector(&seg_data, left_prefix);
    let v_right_raw = lz_right.get_coding_cost_vector(&seg_data, right_prefix);
    let mut v_left = v_left_raw.clone();
    let mut v_right = v_right_raw.clone();

    // Apply cumulative sums as streaming splitter does
    // Left: always forward cumulative sum
    let mut sum = 0u32;
    for c in v_left.iter_mut() {
        sum = sum.saturating_add(*c);
        *c = sum;
    }

    // Right: if suffix costs (prefix=false), we need reverse cumulative sum (right to left)
    // Otherwise (prefix=true), forward cumulative then reverse the vector
    if !right_prefix {
        let mut acc = 0u32;
        for c in v_right.iter_mut().rev() {
            acc = acc.saturating_add(*c);
            *c = acc;
        }
    } else {
        let mut acc = 0u32;
        for c in v_right.iter_mut() {
            acc = acc.saturating_add(*c);
            *c = acc;
        }
        v_right.reverse();
    }

    // Find best split pos by minimizing left[i] + right[i]
    let mut best_sum = u64::MAX;
    let mut best_pos = 0usize;
    for i in 0..v_left.len().min(v_right.len()) {
        let s = (v_left[i] as u64) + (v_right[i] as u64);
        if s < best_sum {
            best_sum = s;
            best_pos = i;
        }
    }

    println!(
        "\n[RUST] RAW left len={} head20={:?}",
        v_left_raw.len(),
        &v_left_raw.iter().take(20).cloned().collect::<Vec<_>>()
    );
    println!(
        "[RUST] RAW left tail20={:?}",
        &v_left_raw
            .iter()
            .rev()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
    );
    println!(
        "[RUST] RAW right len={} head20={:?}",
        v_right_raw.len(),
        &v_right_raw.iter().take(20).cloned().collect::<Vec<_>>()
    );
    println!(
        "[RUST] RAW right tail20={:?}",
        &v_right_raw
            .iter()
            .rev()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
    );

    println!(
        "[RUST] CUM left head20={:?}",
        &v_left.iter().take(20).cloned().collect::<Vec<_>>()
    );
    println!(
        "[RUST] CUM left tail20={:?}",
        &v_left
            .iter()
            .rev()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
    );
    println!(
        "[RUST] CUM right head20={:?}",
        &v_right.iter().take(20).cloned().collect::<Vec<_>>()
    );
    println!(
        "[RUST] CUM right tail20={:?}",
        &v_right
            .iter()
            .rev()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
    );
    println!("[RUST] best_pos={} best_sum={}", best_pos, best_sum);

    // If FFI is available, compute C++ best split and show seg2_start
    #[cfg(feature = "ffi_cost")]
    {
        let flm = {
            // left group is usually (k1, k_mid), so flm indicates if k1 < k_mid
            // We don't have k-mers here; assume forward orientation for display only
            true
        };
        let mlb = true;
        if let Some((ffi_best_pos, ffi_seg2_start)) = ragc_core::ragc_ffi::best_split(
            &left_ref,
            &right_ref,
            &seg_data,
            dec.min_match_len as u32,
            dec.kmer_length as u32,
            flm,
            mlb,
        ) {
            println!(
                "[FFI] best_pos={} seg2_start={}",
                ffi_best_pos, ffi_seg2_start
            );
        } else {
            println!("[FFI] best split unavailable (missing refs?)");
        }
    }

    println!("\nTo tweak inputs: change left/right groups to neighbor segment group_ids or adjust prefix flags to match orientation.");
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

#[allow(clippy::too_many_arguments)]
fn getrange_command(
    archive: PathBuf,
    sample: String,
    contig: String,
    start: usize,
    end: Option<usize>,
    output: Option<PathBuf>,
    format: String,
    verbosity: u32,
) -> Result<()> {
    use ragc_core::CNV_NUM;

    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity };
    let mut decompressor = Decompressor::open(archive_str, config)?;

    // Get end position (use contig length if not specified)
    let end = match end {
        Some(e) => e,
        None => decompressor.get_contig_length(&sample, &contig)?,
    };

    if verbosity > 0 {
        eprintln!(
            "Extracting {}:{} range {}..{} ({} bp)",
            sample,
            contig,
            start,
            end,
            end.saturating_sub(start)
        );
    }

    // Extract the range
    let sequence = decompressor.get_contig_range(&sample, &contig, start, end)?;

    // Convert to ASCII
    let ascii_seq: Vec<u8> = sequence
        .iter()
        .map(|&b| if b < 16 { CNV_NUM[b as usize] } else { b'N' })
        .collect();

    // Output
    let output_data = if format == "raw" {
        ascii_seq
    } else {
        // FASTA format
        let header = format!(">{}:{} {}:{}-{}\n", sample, contig, contig, start, end);
        let mut fasta = header.into_bytes();
        // Wrap at 80 characters
        for chunk in ascii_seq.chunks(80) {
            fasta.extend_from_slice(chunk);
            fasta.push(b'\n');
        }
        fasta
    };

    if let Some(output_path) = output {
        std::fs::write(output_path, &output_data)?;
    } else {
        io::stdout().write_all(&output_data)?;
    }

    decompressor.close()?;
    Ok(())
}

fn ctglen_command(archive: PathBuf, sample: String, contig: String) -> Result<()> {
    let archive_str = archive
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid archive path"))?;

    let config = DecompressorConfig { verbosity: 0 };
    let mut decompressor = Decompressor::open(archive_str, config)?;

    let length = decompressor.get_contig_length(&sample, &contig)?;
    println!("{}", length);

    decompressor.close()?;
    Ok(())
}
