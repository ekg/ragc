// AGC Archive Inspection Tool
// Extracts detailed information about archive structure for debugging and comparison

use anyhow::Result;
use ragc_core::{Decompressor, DecompressorConfig};
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct InspectConfig {
    pub verbosity: u32,
    pub show_groups: bool,
    pub show_segments: bool,
    pub group_id_filter: Option<u32>,
}

impl Default for InspectConfig {
    fn default() -> Self {
        InspectConfig {
            verbosity: 0,
            show_groups: true,
            show_segments: false,
            group_id_filter: None,
        }
    }
}

pub fn inspect_archive(archive_path: PathBuf, config: InspectConfig) -> Result<()> {
    let dec_config = DecompressorConfig {
        verbosity: config.verbosity,
    };

    let mut decompressor = Decompressor::open(
        archive_path.to_str().ok_or_else(|| anyhow::anyhow!("Invalid path"))?,
        dec_config,
    )?;

    println!("=== AGC Archive Inspection ===");
    println!("Archive: {}", archive_path.display());
    println!();

    // Show basic statistics
    let samples = decompressor.list_samples();
    println!("Samples: {}", samples.len());

    if config.show_groups {
        inspect_groups(&mut decompressor, &config)?;
    }

    if config.show_segments {
        inspect_segments(&mut decompressor, &config)?;
    }

    Ok(())
}

fn inspect_groups(decompressor: &mut Decompressor, config: &InspectConfig) -> Result<()> {
    println!("\n=== Group Statistics ===");

    let group_stats = decompressor.get_group_statistics()?;

    if let Some(filter_gid) = config.group_id_filter {
        // Show only filtered group
        if let Some(stats) = group_stats.iter().find(|s| s.0 == filter_gid) {
            println!("Group {}: total={}, refs={}, deltas={}",
                     stats.0, stats.1, stats.2, stats.3);
        } else {
            println!("Group {} not found", filter_gid);
        }
    } else {
        // Show summary
        println!("Total unique groups: {}", group_stats.len());
        println!();
        println!("GroupID  Total  Refs  Deltas");
        println!("-------  -----  ----  ------");
        for (gid, total, refs, deltas) in &group_stats {
            println!("{:7}  {:5}  {:4}  {:6}", gid, total, refs, deltas);
        }

        // Summary statistics
        let total_segments: usize = group_stats.iter().map(|s| s.1).sum();
        let total_refs: usize = group_stats.iter().map(|s| s.2).sum();
        let total_deltas: usize = group_stats.iter().map(|s| s.3).sum();

        println!();
        println!("Summary:");
        println!("  Total segments:    {}", total_segments);
        println!("  Reference segments: {}", total_refs);
        println!("  Delta segments:     {}", total_deltas);
    }

    Ok(())
}

fn inspect_segments(decompressor: &mut Decompressor, config: &InspectConfig) -> Result<()> {
    println!("\n=== Segment Details ===");

    let all_segments = decompressor.get_all_segments()?;

    for (sample_name, contig_name, segments) in &all_segments {
        // Filter by group if requested
        let filtered_segments: Vec<_> = if let Some(filter_gid) = config.group_id_filter {
            segments.iter().filter(|s| s.group_id == filter_gid).collect()
        } else {
            segments.iter().collect()
        };

        if filtered_segments.is_empty() {
            continue;
        }

        println!("\n{}/{}", sample_name, contig_name);
        println!("  Segments: {}", filtered_segments.len());
        for (i, seg) in filtered_segments.iter().enumerate() {
            println!("    [{}] group={}, in_group={}, rc={}, len={}",
                     i, seg.group_id, seg.in_group_id, seg.is_rev_comp, seg.raw_length);
        }
    }

    Ok(())
}
