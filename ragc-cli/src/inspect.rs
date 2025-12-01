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
    pub sample_filter: Option<String>,
    pub contig_filter: Option<String>,
    pub segment_index: Option<usize>,
    pub show_single_segment_groups: bool,
    pub show_segment_layout: bool,
    pub show_compression: bool,
}

impl Default for InspectConfig {
    fn default() -> Self {
        InspectConfig {
            verbosity: 0,
            show_groups: true,
            show_segments: false,
            group_id_filter: None,
            sample_filter: None,
            contig_filter: None,
            segment_index: None,
            show_single_segment_groups: false,
            show_segment_layout: false,
            show_compression: false,
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

    // Check if this is a segment lookup query
    if let (Some(sample), Some(contig), Some(idx)) =
        (&config.sample_filter, &config.contig_filter, config.segment_index) {
        lookup_segment(&mut decompressor, sample, contig, idx)?;
        return Ok(());
    }

    if config.show_single_segment_groups {
        show_single_segment_groups(&mut decompressor)?;
        return Ok(());
    }

    if config.show_segment_layout {
        show_segment_layout(&mut decompressor)?;
        return Ok(());
    }

    if config.show_compression {
        show_compression_stats(&decompressor)?;
        return Ok(());
    }

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

/// Look up which group a specific segment ended up in
fn lookup_segment(decompressor: &mut Decompressor, sample: &str, contig: &str, index: usize) -> Result<()> {
    println!("\n=== Segment Lookup ===");
    println!("Sample: {}", sample);
    println!("Contig: {}", contig);
    println!("Index:  {}", index);
    println!();

    let all_segments = decompressor.get_all_segments()?;

    for (sample_name, contig_name, segments) in &all_segments {
        if sample_name == sample && contig_name == contig {
            if index < segments.len() {
                let seg = &segments[index];
                println!("Found segment:");
                println!("  Group ID:    {}", seg.group_id);
                println!("  In-group ID: {}", seg.in_group_id);
                println!("  Rev comp:    {}", seg.is_rev_comp);
                println!("  Length:      {}", seg.raw_length);
                return Ok(());
            } else {
                println!("ERROR: Index {} out of range (contig has {} segments)", index, segments.len());
                return Ok(());
            }
        }
    }

    println!("ERROR: Segment not found (check sample/contig names)");
    Ok(())
}

/// Show all single-segment groups for comparison
fn show_single_segment_groups(decompressor: &mut Decompressor) -> Result<()> {
    println!("\n=== Single-Segment Groups ===");

    let group_stats = decompressor.get_group_statistics()?;
    let all_segments = decompressor.get_all_segments()?;

    // Find all groups with exactly 1 segment
    let single_groups: Vec<_> = group_stats.iter()
        .filter(|(_, total, _, _)| *total == 1)
        .map(|(gid, _, _, _)| *gid)
        .collect();

    println!("Total single-segment groups: {}", single_groups.len());
    println!();

    for group_id in &single_groups {
        // Find the segment in this group
        for (sample_name, contig_name, segments) in &all_segments {
            for (idx, seg) in segments.iter().enumerate() {
                if seg.group_id == *group_id {
                    println!("Group {}: {}/{} segment[{}] len={} rc={}",
                             group_id, sample_name, contig_name, idx, seg.raw_length, seg.is_rev_comp);
                }
            }
        }
    }

    Ok(())
}

/// Show segment layout in CSV format for comparison between implementations
fn show_segment_layout(decompressor: &mut Decompressor) -> Result<()> {
    let all_segments = decompressor.get_all_segments()?;

    // CSV header
    println!("sample,contig,segment_index,group_id,in_group_id,raw_length");

    for (sample_name, contig_name, segments) in &all_segments {
        for (seg_idx, seg) in segments.iter().enumerate() {
            println!("{},{},{},{},{},{}",
                sample_name, contig_name, seg_idx,
                seg.group_id, seg.in_group_id, seg.raw_length);
        }
    }

    Ok(())
}

/// Show compression statistics for all streams
fn show_compression_stats(decompressor: &Decompressor) -> Result<()> {
    println!("\n=== Compression Statistics ===");

    let stats = decompressor.get_compression_stats();

    // Calculate totals
    let total_raw: u64 = stats.iter().map(|(_, raw, _, _)| raw).sum();
    let total_packed: u64 = stats.iter().map(|(_, _, packed, _)| packed).sum();
    let total_parts: usize = stats.iter().map(|(_, _, _, parts)| parts).sum();

    println!("\nOverall:");
    println!("  Total streams:     {}", stats.len());
    println!("  Total parts:       {}", total_parts);
    println!("  Total raw size:    {} bytes ({:.2} MB)", total_raw, total_raw as f64 / 1_000_000.0);
    println!("  Total packed size: {} bytes ({:.2} MB)", total_packed, total_packed as f64 / 1_000_000.0);
    if total_raw > 0 {
        println!("  Compression ratio: {:.2}:1", total_raw as f64 / total_packed as f64);
        println!("  Space saved:       {:.2}%", (1.0 - (total_packed as f64 / total_raw as f64)) * 100.0);
    }

    // Group streams by type
    let mut ref_streams = Vec::new();
    let mut delta_streams = Vec::new();
    let mut other_streams = Vec::new();

    for (name, raw, packed, parts) in &stats {
        if name.ends_with(".ref") {
            ref_streams.push((name, *raw, *packed, *parts));
        } else if name.ends_with(".delta") {
            delta_streams.push((name, *raw, *packed, *parts));
        } else {
            other_streams.push((name, *raw, *packed, *parts));
        }
    }

    // Show reference streams
    if !ref_streams.is_empty() {
        println!("\n=== Reference Streams ({}) ===", ref_streams.len());
        let ref_raw: u64 = ref_streams.iter().map(|(_, raw, _, _)| raw).sum();
        let ref_packed: u64 = ref_streams.iter().map(|(_, _, packed, _)| packed).sum();
        println!("Total: {} bytes -> {} bytes (ratio: {:.2}:1)",
                 ref_raw, ref_packed,
                 if ref_packed > 0 { ref_raw as f64 / ref_packed as f64 } else { 0.0 });
    }

    // Show delta streams
    if !delta_streams.is_empty() {
        println!("\n=== Delta Streams ({}) ===", delta_streams.len());
        let delta_raw: u64 = delta_streams.iter().map(|(_, raw, _, _)| raw).sum();
        let delta_packed: u64 = delta_streams.iter().map(|(_, _, packed, _)| packed).sum();
        println!("Total: {} bytes -> {} bytes (ratio: {:.2}:1)",
                 delta_raw, delta_packed,
                 if delta_packed > 0 { delta_raw as f64 / delta_packed as f64 } else { 0.0 });
    }

    // Show metadata/other streams
    if !other_streams.is_empty() {
        println!("\n=== Metadata Streams ({}) ===", other_streams.len());
        for (name, raw, packed, _parts) in &other_streams {
            println!("  {}: {} -> {} bytes", name, raw, packed);
        }
    }

    Ok(())
}
