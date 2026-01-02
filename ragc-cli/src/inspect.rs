// AGC Archive Inspection Tool
// Extracts detailed information about archive structure for debugging and comparison

use anyhow::Result;
use ragc_core::{Decompressor, DecompressorConfig};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

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
    pub show_pack_layout: bool,
    pub show_compression: bool,
    pub compare_with: Option<PathBuf>,
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
            show_pack_layout: false,
            show_compression: false,
            compare_with: None,
        }
    }
}

pub fn inspect_archive(archive_path: PathBuf, config: InspectConfig) -> Result<()> {
    let dec_config = DecompressorConfig {
        verbosity: config.verbosity,
    };

    let mut decompressor = Decompressor::open(
        archive_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid path"))?,
        dec_config,
    )?;

    // If comparing with another archive, run comparison and exit
    if let Some(compare_path) = &config.compare_with {
        return compare_archives(&archive_path, compare_path, &config);
    }

    println!("=== AGC Archive Inspection ===");
    println!("Archive: {}", archive_path.display());
    println!();

    // Show basic statistics
    let samples = decompressor.list_samples();
    println!("Samples: {}", samples.len());

    // Check if this is a segment lookup query
    if let (Some(sample), Some(contig), Some(idx)) = (
        &config.sample_filter,
        &config.contig_filter,
        config.segment_index,
    ) {
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

    if config.show_pack_layout {
        show_pack_layout(&mut decompressor)?;
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
            println!(
                "Group {}: total={}, refs={}, deltas={}",
                stats.0, stats.1, stats.2, stats.3
            );
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
            segments
                .iter()
                .filter(|s| s.group_id == filter_gid)
                .collect()
        } else {
            segments.iter().collect()
        };

        if filtered_segments.is_empty() {
            continue;
        }

        println!("\n{}/{}", sample_name, contig_name);
        println!("  Segments: {}", filtered_segments.len());
        for (i, seg) in filtered_segments.iter().enumerate() {
            println!(
                "    [{}] group={}, in_group={}, rc={}, len={}",
                i, seg.group_id, seg.in_group_id, seg.is_rev_comp, seg.raw_length
            );
        }
    }

    Ok(())
}

/// Look up which group a specific segment ended up in
fn lookup_segment(
    decompressor: &mut Decompressor,
    sample: &str,
    contig: &str,
    index: usize,
) -> Result<()> {
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
                println!(
                    "ERROR: Index {} out of range (contig has {} segments)",
                    index,
                    segments.len()
                );
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
    let single_groups: Vec<_> = group_stats
        .iter()
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
                    println!(
                        "Group {}: {}/{} segment[{}] len={} rc={}",
                        group_id, sample_name, contig_name, idx, seg.raw_length, seg.is_rev_comp
                    );
                }
            }
        }
    }

    Ok(())
}

/// Show segment layout in CSV format for comparison between implementations
fn show_segment_layout(decompressor: &mut Decompressor) -> Result<()> {
    let all_segments = decompressor.get_all_segments()?;

    // CSV header - includes is_rev_comp for orientation comparison
    println!("sample,contig,segment_index,group_id,in_group_id,is_rev_comp,raw_length");

    for (sample_name, contig_name, segments) in &all_segments {
        for (seg_idx, seg) in segments.iter().enumerate() {
            println!(
                "{},{},{},{},{},{},{}",
                sample_name,
                contig_name,
                seg_idx,
                seg.group_id,
                seg.in_group_id,
                seg.is_rev_comp,
                seg.raw_length
            );
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
    println!(
        "  Total raw size:    {} bytes ({:.2} MB)",
        total_raw,
        total_raw as f64 / 1_000_000.0
    );
    println!(
        "  Total packed size: {} bytes ({:.2} MB)",
        total_packed,
        total_packed as f64 / 1_000_000.0
    );
    if total_raw > 0 {
        println!(
            "  Compression ratio: {:.2}:1",
            total_raw as f64 / total_packed as f64
        );
        println!(
            "  Space saved:       {:.2}%",
            (1.0 - (total_packed as f64 / total_raw as f64)) * 100.0
        );
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
        println!(
            "Total: {} bytes -> {} bytes (ratio: {:.2}:1)",
            ref_raw,
            ref_packed,
            if ref_packed > 0 {
                ref_raw as f64 / ref_packed as f64
            } else {
                0.0
            }
        );
    }

    // Show delta streams
    if !delta_streams.is_empty() {
        println!("\n=== Delta Streams ({}) ===", delta_streams.len());
        let delta_raw: u64 = delta_streams.iter().map(|(_, raw, _, _)| raw).sum();
        let delta_packed: u64 = delta_streams.iter().map(|(_, _, packed, _)| packed).sum();
        println!(
            "Total: {} bytes -> {} bytes (ratio: {:.2}:1)",
            delta_raw,
            delta_packed,
            if delta_packed > 0 {
                delta_raw as f64 / delta_packed as f64
            } else {
                0.0
            }
        );
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

/// Show pack layout - how segments are distributed across ZSTD packs
fn show_pack_layout(decompressor: &mut Decompressor) -> Result<()> {
    let all_segments = decompressor.get_all_segments()?;
    let compression_stats = decompressor.get_compression_stats();

    // CSV header
    println!("group_id,pack_id,segment_count,sample,contig,segment_indices");

    // Build a map of group_id -> segments
    let mut group_segments: std::collections::HashMap<u32, Vec<(String, String, usize, usize)>> =
        std::collections::HashMap::new();

    for (sample_name, contig_name, segments) in &all_segments {
        for (seg_idx, seg) in segments.iter().enumerate() {
            group_segments.entry(seg.group_id).or_default().push((
                sample_name.clone(),
                contig_name.clone(),
                seg_idx,
                seg.in_group_id as usize,
            ));
        }
    }

    // Find delta streams and their pack counts
    for (stream_name, _raw, _packed, num_parts) in &compression_stats {
        // AGC v3 uses "xGd" (short name with 'd' suffix for delta)
        // AGC v2 uses "seg-1234-delta" (long name with "-delta" suffix)
        let is_delta = stream_name.ends_with('d') && stream_name.starts_with('x');
        if !is_delta {
            continue;
        }

        // Extract group_id from stream name
        // AGC v3: "xGd" -> decode base64 "G" -> group_id
        let group_id: u32 = if stream_name.starts_with('x') {
            // Remove 'x' prefix and 'd' suffix, then decode base64
            let b64_str = stream_name.trim_start_matches('x').trim_end_matches('d');
            base64_decode(b64_str)
        } else {
            // AGC v2: "seg-1234-delta" -> 1234
            stream_name
                .trim_start_matches("seg-")
                .trim_end_matches("-delta")
                .parse()
                .unwrap_or(0)
        };

        if let Some(segments) = group_segments.get(&group_id) {
            // Sort segments by in_group_id
            let mut sorted_segments = segments.clone();
            sorted_segments.sort_by_key(|s| s.3);

            // Group into packs of 50 (PACK_CARDINALITY)
            const PACK_CARDINALITY: usize = 50;
            for pack_id in 0..*num_parts {
                let pack_start = pack_id * PACK_CARDINALITY;
                let pack_end = (pack_start + PACK_CARDINALITY).min(sorted_segments.len());

                if pack_start >= sorted_segments.len() {
                    break;
                }

                let pack_segments = &sorted_segments[pack_start..pack_end];

                // Format segment info: "sample#contig[idx]"
                let seg_info: Vec<String> = pack_segments
                    .iter()
                    .map(|(sample, contig, idx, _)| format!("{}#{}[{}]", sample, contig, idx))
                    .collect();

                println!(
                    "{},{},{},\"{}\",\"{}\",\"{}\"",
                    group_id,
                    pack_id,
                    pack_segments.len(),
                    pack_segments[0].0, // First sample name
                    pack_segments[0].1, // First contig name
                    seg_info.join(";")
                );
            }
        }
    }

    Ok(())
}

/// Decode base64 stream name to group_id (matches C++ int_to_base64)
fn base64_decode(s: &str) -> u32 {
    const DIGITS: &str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_#";

    let mut result = 0u32;
    let mut multiplier = 1u32;

    for c in s.chars() {
        if let Some(digit_value) = DIGITS.find(c) {
            result += (digit_value as u32) * multiplier;
            multiplier *= 64;
        }
    }

    result
}
// Archive comparison functions for debugging differences between implementations
// These functions systematically compare two archives to answer:
// 1. Segment grouping: Which segments are assigned to different groups?
// 2. Segment boundaries: Which segments have different lengths/boundaries?
// 3. Reference selection: Which groups have different reference segments?
// 4. Pack structure: How do ZSTD pack boundaries differ?

/// Main comparison entry point - calls all 4 comparison functions
fn compare_archives(
    archive1_path: &Path,
    archive2_path: &Path,
    config: &InspectConfig,
) -> Result<()> {
    let dec_config = DecompressorConfig {
        verbosity: config.verbosity,
    };

    let mut dec1 = Decompressor::open(
        archive1_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid path"))?,
        dec_config.clone(),
    )?;

    let mut dec2 = Decompressor::open(
        archive2_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("Invalid path"))?,
        dec_config,
    )?;

    println!("=== AGC Archive Comparison ===");
    println!("Archive 1: {}", archive1_path.display());
    println!("Archive 2: {}", archive2_path.display());
    println!();

    // Run all 5 comparisons
    compare_segment_boundaries(&mut dec1, &mut dec2, archive1_path, archive2_path)?;
    compare_grouping(&mut dec1, &mut dec2, archive1_path, archive2_path)?;
    compare_references(&mut dec1, &mut dec2, archive1_path, archive2_path)?;
    compare_orientations(&mut dec1, &mut dec2, archive1_path, archive2_path)?;
    compare_packs(&mut dec1, &mut dec2, archive1_path, archive2_path)?;

    Ok(())
}

/// 1. Compare segment boundaries: Are segments split at the same positions?
fn compare_segment_boundaries(
    dec1: &mut Decompressor,
    dec2: &mut Decompressor,
    arch1_name: &Path,
    arch2_name: &Path,
) -> Result<()> {
    println!("\n=== 1. SEGMENT BOUNDARIES COMPARISON ===");

    let segs1 = dec1.get_all_segments()?;
    let segs2 = dec2.get_all_segments()?;

    // Build map: (sample, contig, seg_idx) -> raw_length
    let mut map1: BTreeMap<(String, String, usize), usize> = BTreeMap::new();
    let mut map2: BTreeMap<(String, String, usize), usize> = BTreeMap::new();

    for (sample, contig, segments) in &segs1 {
        for (idx, seg) in segments.iter().enumerate() {
            map1.insert(
                (sample.clone(), contig.clone(), idx),
                seg.raw_length as usize,
            );
        }
    }

    for (sample, contig, segments) in &segs2 {
        for (idx, seg) in segments.iter().enumerate() {
            map2.insert(
                (sample.clone(), contig.clone(), idx),
                seg.raw_length as usize,
            );
        }
    }

    // Find differences
    let mut diff_lengths = Vec::new();
    let mut only_in_1 = Vec::new();
    let mut only_in_2 = Vec::new();

    let all_keys: BTreeSet<_> = map1.keys().chain(map2.keys()).collect();

    for key in all_keys {
        match (map1.get(key), map2.get(key)) {
            (Some(len1), Some(len2)) => {
                if len1 != len2 {
                    diff_lengths.push((key.clone(), *len1, *len2));
                }
            }
            (Some(len1), None) => {
                only_in_1.push((key.clone(), *len1));
            }
            (None, Some(len2)) => {
                only_in_2.push((key.clone(), *len2));
            }
            (None, None) => unreachable!(),
        }
    }

    println!("Total segments:");
    println!(
        "  {}: {} segments",
        arch1_name.file_name().unwrap().to_str().unwrap(),
        map1.len()
    );
    println!(
        "  {}: {} segments",
        arch2_name.file_name().unwrap().to_str().unwrap(),
        map2.len()
    );
    println!();

    if diff_lengths.is_empty() && only_in_1.is_empty() && only_in_2.is_empty() {
        println!("✅ IDENTICAL segment boundaries!");
    } else {
        println!("❌ DIFFERENT segment splitting:");
        println!("  Note: Segments are compared by INDEX, not genomic position.");
        println!(
            "  When one archive splits a segment and the other doesn't, indices become misaligned."
        );
        println!();
        println!("  Segment index mismatches: {}", diff_lengths.len());
        println!("  Extra segments in archive 1: {}", only_in_1.len());
        println!("  Extra segments in archive 2: {}", only_in_2.len());

        // Show first few differences
        if !diff_lengths.is_empty() {
            println!("\nFirst 10 index mismatches (may indicate different split decisions):");
            for ((sample, contig, idx), len1, len2) in diff_lengths.iter().take(10) {
                println!(
                    "  {} {} seg[{}]: {} vs {} bytes (diff: {})",
                    sample,
                    contig,
                    idx,
                    len1,
                    len2,
                    (*len1 as i64 - *len2 as i64).abs()
                );
            }
        }

        if !only_in_1.is_empty() {
            println!("\nFirst 5 segments only in archive 1:");
            for ((sample, contig, idx), len) in only_in_1.iter().take(5) {
                println!("  {} {} seg[{}]: {} bytes", sample, contig, idx, len);
            }
        }

        if !only_in_2.is_empty() {
            println!("\nFirst 5 segments only in archive 2:");
            for ((sample, contig, idx), len) in only_in_2.iter().take(5) {
                println!("  {} {} seg[{}]: {} bytes", sample, contig, idx, len);
            }
        }
    }

    Ok(())
}

/// 2. Compare logical grouping: Are segments assigned to the same groups?
fn compare_grouping(
    dec1: &mut Decompressor,
    dec2: &mut Decompressor,
    arch1_name: &Path,
    arch2_name: &Path,
) -> Result<()> {
    println!("\n=== 2. SEGMENT GROUPING COMPARISON ===");

    let segs1 = dec1.get_all_segments()?;
    let segs2 = dec2.get_all_segments()?;

    // Build map: (sample, contig, seg_idx) -> group_id
    let mut map1: BTreeMap<(String, String, usize), u32> = BTreeMap::new();
    let mut map2: BTreeMap<(String, String, usize), u32> = BTreeMap::new();

    // Also build reverse map: group_id -> Set of segments
    let mut groups1: BTreeMap<u32, BTreeSet<(String, String, usize)>> = BTreeMap::new();
    let mut groups2: BTreeMap<u32, BTreeSet<(String, String, usize)>> = BTreeMap::new();

    for (sample, contig, segments) in &segs1 {
        for (idx, seg) in segments.iter().enumerate() {
            let key = (sample.clone(), contig.clone(), idx);
            map1.insert(key.clone(), seg.group_id);
            groups1.entry(seg.group_id).or_default().insert(key);
        }
    }

    for (sample, contig, segments) in &segs2 {
        for (idx, seg) in segments.iter().enumerate() {
            let key = (sample.clone(), contig.clone(), idx);
            map2.insert(key.clone(), seg.group_id);
            groups2.entry(seg.group_id).or_default().insert(key);
        }
    }

    // Find segments in groups with different membership
    let mut mismatched_segments = Vec::new();

    for (key, gid1) in &map1 {
        if let Some(gid2) = map2.get(key) {
            // Check if these groups have identical membership
            if let (Some(members1), Some(members2)) = (groups1.get(gid1), groups2.get(gid2)) {
                if members1 != members2 {
                    mismatched_segments.push((
                        key.clone(),
                        *gid1,
                        *gid2,
                        members1.len(),
                        members2.len(),
                    ));
                }
            }
        }
    }

    println!("Total groups:");
    println!(
        "  {}: {} groups",
        arch1_name.file_name().unwrap().to_str().unwrap(),
        groups1.len()
    );
    println!(
        "  {}: {} groups",
        arch2_name.file_name().unwrap().to_str().unwrap(),
        groups2.len()
    );
    println!();

    if mismatched_segments.is_empty() {
        println!("✅ IDENTICAL logical grouping!");
        println!("   All segments are assigned to groups with identical membership.");
        println!("   (Group IDs may differ, but that's just labeling)");
    } else {
        println!("❌ DIFFERENT logical grouping:");
        println!(
            "  {} segments are in groups with different membership",
            mismatched_segments.len()
        );

        println!("\nFirst 10 segments with different group membership:");
        for ((sample, contig, idx), gid1, gid2, size1, size2) in mismatched_segments.iter().take(10)
        {
            println!("  {} {} seg[{}]:", sample, contig, idx);
            println!("    Archive 1: group {} ({} segments)", gid1, size1);
            println!("    Archive 2: group {} ({} segments)", gid2, size2);
        }
    }

    Ok(())
}

/// 3. Compare reference selection: Are the same segments chosen as references?
fn compare_references(
    dec1: &mut Decompressor,
    dec2: &mut Decompressor,
    arch1_name: &Path,
    arch2_name: &Path,
) -> Result<()> {
    println!("\n=== 3. REFERENCE SEGMENT SELECTION COMPARISON ===");

    let segs1 = dec1.get_all_segments()?;
    let segs2 = dec2.get_all_segments()?;

    // Build map: group_id -> (sample, contig, seg_idx) of reference (in_group_id=0)
    let mut refs1: BTreeMap<u32, (String, String, usize)> = BTreeMap::new();
    let mut refs2: BTreeMap<u32, (String, String, usize)> = BTreeMap::new();

    for (sample, contig, segments) in &segs1 {
        for (idx, seg) in segments.iter().enumerate() {
            if seg.in_group_id == 0 {
                refs1.insert(seg.group_id, (sample.clone(), contig.clone(), idx));
            }
        }
    }

    for (sample, contig, segments) in &segs2 {
        for (idx, seg) in segments.iter().enumerate() {
            if seg.in_group_id == 0 {
                refs2.insert(seg.group_id, (sample.clone(), contig.clone(), idx));
            }
        }
    }

    // For logical comparison, we need to map groups by their membership
    // Build reverse map: Set of members -> group_id
    let mut group_members1: BTreeMap<u32, BTreeSet<(String, String, usize)>> = BTreeMap::new();
    let mut group_members2: BTreeMap<u32, BTreeSet<(String, String, usize)>> = BTreeMap::new();

    for (sample, contig, segments) in &segs1 {
        for (idx, seg) in segments.iter().enumerate() {
            group_members1.entry(seg.group_id).or_default().insert((
                sample.clone(),
                contig.clone(),
                idx,
            ));
        }
    }

    for (sample, contig, segments) in &segs2 {
        for (idx, seg) in segments.iter().enumerate() {
            group_members2.entry(seg.group_id).or_default().insert((
                sample.clone(),
                contig.clone(),
                idx,
            ));
        }
    }

    // Build membership -> group_id map for archive 2
    let mut members_to_gid2: BTreeMap<BTreeSet<(String, String, usize)>, u32> = BTreeMap::new();
    for (gid, members) in &group_members2 {
        members_to_gid2.insert(members.clone(), *gid);
    }

    // Compare references for groups with identical membership
    let mut diff_refs = Vec::new();

    for (gid1, members1) in &group_members1 {
        if let Some(gid2) = members_to_gid2.get(members1) {
            // Same logical group - compare references
            let ref1 = refs1.get(gid1);
            let ref2 = refs2.get(gid2);

            match (ref1, ref2) {
                (Some(r1), Some(r2)) if r1 != r2 => {
                    diff_refs.push((*gid1, *gid2, r1.clone(), r2.clone(), members1.len()));
                }
                (Some(_), None) => {
                    // Archive 1 has reference, archive 2 doesn't
                }
                (None, Some(_)) => {
                    // Archive 2 has reference, archive 1 doesn't
                }
                _ => {}
            }
        }
    }

    println!("Total groups with references:");
    println!(
        "  {}: {} groups",
        arch1_name.file_name().unwrap().to_str().unwrap(),
        refs1.len()
    );
    println!(
        "  {}: {} groups",
        arch2_name.file_name().unwrap().to_str().unwrap(),
        refs2.len()
    );
    println!();

    if diff_refs.is_empty() {
        println!("✅ IDENTICAL reference selection!");
        println!("   All groups choose the same segment as reference.");
    } else {
        println!("❌ DIFFERENT reference selection:");
        println!(
            "  {} groups have different reference segments",
            diff_refs.len()
        );

        println!("\nFirst 10 groups with different references:");
        for (gid1, gid2, (s1, c1, i1), (s2, c2, i2), size) in diff_refs.iter().take(10) {
            println!("  Group {} / {} ({} segments):", gid1, gid2, size);
            println!("    Archive 1 ref: {} {} seg[{}]", s1, c1, i1);
            println!("    Archive 2 ref: {} {} seg[{}]", s2, c2, i2);
        }
    }

    Ok(())
}

/// 4. Compare segment orientations: Are segments stored in the same orientation (rev comp)?
fn compare_orientations(
    dec1: &mut Decompressor,
    dec2: &mut Decompressor,
    arch1_name: &Path,
    arch2_name: &Path,
) -> Result<()> {
    println!("\n=== 4. SEGMENT ORIENTATION COMPARISON ===");

    let segs1 = dec1.get_all_segments()?;
    let segs2 = dec2.get_all_segments()?;

    // Build map: (sample, contig, seg_idx) -> is_rev_comp
    let mut map1: BTreeMap<(String, String, usize), bool> = BTreeMap::new();
    let mut map2: BTreeMap<(String, String, usize), bool> = BTreeMap::new();

    for (sample, contig, segments) in &segs1 {
        for (idx, seg) in segments.iter().enumerate() {
            map1.insert((sample.clone(), contig.clone(), idx), seg.is_rev_comp);
        }
    }

    for (sample, contig, segments) in &segs2 {
        for (idx, seg) in segments.iter().enumerate() {
            map2.insert((sample.clone(), contig.clone(), idx), seg.is_rev_comp);
        }
    }

    // Find differences
    let mut diff_orientations = Vec::new();

    for (key, rc1) in &map1 {
        if let Some(rc2) = map2.get(key) {
            if rc1 != rc2 {
                diff_orientations.push((key.clone(), *rc1, *rc2));
            }
        }
    }

    // Count orientation stats
    let rc_true_1 = map1.values().filter(|&&rc| rc).count();
    let rc_true_2 = map2.values().filter(|&&rc| rc).count();

    println!("Orientation statistics:");
    println!(
        "  {}: {} rev_comp=true out of {} total",
        arch1_name.file_name().unwrap().to_str().unwrap(),
        rc_true_1,
        map1.len()
    );
    println!(
        "  {}: {} rev_comp=true out of {} total",
        arch2_name.file_name().unwrap().to_str().unwrap(),
        rc_true_2,
        map2.len()
    );
    println!();

    if diff_orientations.is_empty() {
        println!("✅ IDENTICAL segment orientations!");
        println!("   All segments have the same is_rev_comp flag.");
    } else {
        println!("❌ DIFFERENT segment orientations:");
        println!(
            "  {} segments have different orientations",
            diff_orientations.len()
        );

        println!("\nFirst 20 segments with different orientations:");
        for ((sample, contig, idx), rc1, rc2) in diff_orientations.iter().take(20) {
            println!(
                "  {} {} seg[{}]: archive1 rc={} vs archive2 rc={}",
                sample, contig, idx, rc1, rc2
            );
        }

        // Show summary by contig
        let mut contig_diffs: BTreeMap<(String, String), usize> = BTreeMap::new();
        for ((sample, contig, _), _, _) in &diff_orientations {
            *contig_diffs
                .entry((sample.clone(), contig.clone()))
                .or_default() += 1;
        }

        println!("\nOrientation differences by contig:");
        for ((sample, contig), count) in &contig_diffs {
            println!("  {}/{}: {} segments differ", sample, contig, count);
        }
    }

    Ok(())
}

/// 5. Compare pack structure: Are segments organized into the same ZSTD packs?
fn compare_packs(
    dec1: &mut Decompressor,
    dec2: &mut Decompressor,
    arch1_name: &Path,
    arch2_name: &Path,
) -> Result<()> {
    println!("\n=== 4. PACK STRUCTURE COMPARISON ===");

    let stats1 = dec1.get_compression_stats();
    let stats2 = dec2.get_compression_stats();

    // Count packs per group (delta streams)
    let mut group_packs1: BTreeMap<u32, usize> = BTreeMap::new();
    let mut group_packs2: BTreeMap<u32, usize> = BTreeMap::new();

    for (stream_name, _raw, _packed, num_parts) in &stats1 {
        if stream_name.ends_with('d') && stream_name.starts_with('x') {
            let b64_str = stream_name.trim_start_matches('x').trim_end_matches('d');
            let group_id = base64_decode(b64_str);
            group_packs1.insert(group_id, *num_parts);
        }
    }

    for (stream_name, _raw, _packed, num_parts) in &stats2 {
        if stream_name.ends_with('d') && stream_name.starts_with('x') {
            let b64_str = stream_name.trim_start_matches('x').trim_end_matches('d');
            let group_id = base64_decode(b64_str);
            group_packs2.insert(group_id, *num_parts);
        }
    }

    // Compare total compression
    let total_raw1: u64 = stats1.iter().map(|(_, raw, _, _)| raw).sum();
    let total_packed1: u64 = stats1.iter().map(|(_, _, packed, _)| packed).sum();
    let total_raw2: u64 = stats2.iter().map(|(_, raw, _, _)| raw).sum();
    let total_packed2: u64 = stats2.iter().map(|(_, _, packed, _)| packed).sum();

    println!("Archive sizes:");
    println!(
        "  {}: {:.2} MB raw → {:.2} MB packed (ratio: {:.2}:1)",
        arch1_name.file_name().unwrap().to_str().unwrap(),
        total_raw1 as f64 / 1_000_000.0,
        total_packed1 as f64 / 1_000_000.0,
        if total_packed1 > 0 {
            total_raw1 as f64 / total_packed1 as f64
        } else {
            0.0
        }
    );
    println!(
        "  {}: {:.2} MB raw → {:.2} MB packed (ratio: {:.2}:1)",
        arch2_name.file_name().unwrap().to_str().unwrap(),
        total_raw2 as f64 / 1_000_000.0,
        total_packed2 as f64 / 1_000_000.0,
        if total_packed2 > 0 {
            total_raw2 as f64 / total_packed2 as f64
        } else {
            0.0
        }
    );
    println!();

    // Find groups with different pack counts
    let all_groups: BTreeSet<_> = group_packs1.keys().chain(group_packs2.keys()).collect();
    let mut diff_packs = Vec::new();

    for gid in all_groups {
        let packs1 = group_packs1.get(gid).copied().unwrap_or(0);
        let packs2 = group_packs2.get(gid).copied().unwrap_or(0);
        if packs1 != packs2 {
            diff_packs.push((*gid, packs1, packs2));
        }
    }

    if diff_packs.is_empty() {
        println!("✅ IDENTICAL pack structure!");
        println!("   All groups have the same number of ZSTD packs.");
    } else {
        println!("❌ DIFFERENT pack structure:");
        println!("  {} groups have different pack counts", diff_packs.len());

        println!("\nFirst 10 groups with different pack counts:");
        for (gid, packs1, packs2) in diff_packs.iter().take(10) {
            println!("  Group {}: {} packs vs {} packs", gid, packs1, packs2);
        }
    }

    // Show size difference
    let size_diff = (total_packed1 as i64 - total_packed2 as i64).abs() as f64 / 1_000_000.0;
    let size_pct = if total_packed2 > 0 {
        ((total_packed1 as f64 - total_packed2 as f64) / total_packed2 as f64 * 100.0).abs()
    } else {
        0.0
    };

    println!("\nSize difference: {:.2} MB ({:.2}%)", size_diff, size_pct);

    Ok(())
}
