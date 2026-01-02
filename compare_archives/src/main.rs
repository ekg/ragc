use ragc_common::archive::Archive;
use ragc_core::decompress_segment_with_marker;
use sha2::{Digest, Sha256};
use std::collections::HashMap;

fn hash_data(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    format!("{:x}", hasher.finalize())[..16].to_string()
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <cpp_agc_file> <ragc_file>", args[0]);
        std::process::exit(1);
    }

    let cpp_file = &args[1];
    let ragc_file = &args[2];

    println!("Opening archives...");

    // Open C++ AGC archive
    let mut cpp_archive = Archive::new_reader();
    cpp_archive.open(cpp_file)?;

    // Open RAGC archive
    let mut ragc_archive = Archive::new_reader();
    ragc_archive.open(ragc_file)?;

    // Get all stream names
    let cpp_streams = cpp_archive.get_stream_names();
    let ragc_streams = ragc_archive.get_stream_names();

    println!("C++ AGC streams: {} total", cpp_streams.len());
    println!("RAGC streams: {} total", ragc_streams.len());

    // Filter for segment streams and create maps
    let is_segment_stream = |s: &&String| {
        s.starts_with("seg-") || (s.starts_with("x") && (s.ends_with("r") || s.ends_with("d")))
    };

    let cpp_seg_streams: Vec<String> = cpp_streams
        .iter()
        .filter(is_segment_stream)
        .cloned()
        .collect();
    let ragc_seg_streams: Vec<String> = ragc_streams
        .iter()
        .filter(is_segment_stream)
        .cloned()
        .collect();

    println!("\nSegment streams:");
    println!("  C++ AGC: {} streams", cpp_seg_streams.len());
    println!("  RAGC: {} streams", ragc_seg_streams.len());

    // Find common streams
    let cpp_set: std::collections::HashSet<_> = cpp_seg_streams.iter().collect();
    let ragc_set: std::collections::HashSet<_> = ragc_seg_streams.iter().collect();

    let mut common_streams: Vec<_> = cpp_set
        .intersection(&ragc_set)
        .map(|s| s.to_string())
        .collect();
    common_streams.sort();

    println!("  Common streams: {}", common_streams.len());

    if common_streams.is_empty() {
        println!("\nNo common segment streams found!");
        println!("\nFirst 10 C++ AGC streams:");
        for stream in cpp_seg_streams.iter().take(10) {
            println!("  {}", stream);
        }
        println!("\nFirst 10 RAGC streams:");
        for stream in ragc_seg_streams.iter().take(10) {
            println!("  {}", stream);
        }
        return Ok(());
    }

    // Process all blocks from common streams
    println!("\nDecompressing and hashing blocks from common streams...");

    let mut cpp_hashes = Vec::new();
    let mut ragc_hashes = Vec::new();

    // Process C++ AGC blocks
    print!("  C++ AGC: ");
    std::io::Write::flush(&mut std::io::stdout()).ok();
    let mut total_blocks = 0;
    for stream_name in &common_streams {
        if let Some(stream_id) = cpp_archive.get_stream_id(stream_name) {
            let num_parts = cpp_archive.get_num_parts(stream_id);
            for part_id in 0..num_parts {
                match cpp_archive.get_part_by_id(stream_id, part_id) {
                    Ok((compressed, metadata)) => {
                        let marker = (metadata & 0xFF) as u8;
                        match decompress_segment_with_marker(&compressed, marker) {
                            Ok(decompressed) => {
                                let hash = hash_data(&decompressed);
                                cpp_hashes.push((
                                    stream_name.clone(),
                                    part_id,
                                    hash,
                                    compressed.len(),
                                    decompressed.len(),
                                    marker,
                                ));
                            }
                            Err(e) => {
                                cpp_hashes.push((
                                    stream_name.clone(),
                                    part_id,
                                    format!("DECOMPRESS_FAILED: {}", e),
                                    compressed.len(),
                                    0,
                                    marker,
                                ));
                            }
                        }
                    }
                    Err(e) => {
                        cpp_hashes.push((
                            stream_name.clone(),
                            part_id,
                            format!("READ_FAILED: {}", e),
                            0,
                            0,
                            0,
                        ));
                    }
                }
                total_blocks += 1;
                if total_blocks % 100 == 0 {
                    print!("{}...", total_blocks);
                    std::io::Write::flush(&mut std::io::stdout()).ok();
                }
            }
        }
    }
    println!("Done! ({} blocks)", cpp_hashes.len());

    // Process RAGC blocks
    print!("  RAGC: ");
    std::io::Write::flush(&mut std::io::stdout()).ok();
    total_blocks = 0;
    for stream_name in &common_streams {
        if let Some(stream_id) = ragc_archive.get_stream_id(stream_name) {
            let num_parts = ragc_archive.get_num_parts(stream_id);
            for part_id in 0..num_parts {
                match ragc_archive.get_part_by_id(stream_id, part_id) {
                    Ok((compressed, metadata)) => {
                        let marker = (metadata & 0xFF) as u8;
                        match decompress_segment_with_marker(&compressed, marker) {
                            Ok(decompressed) => {
                                let hash = hash_data(&decompressed);
                                ragc_hashes.push((
                                    stream_name.clone(),
                                    part_id,
                                    hash,
                                    compressed.len(),
                                    decompressed.len(),
                                    marker,
                                ));
                            }
                            Err(e) => {
                                ragc_hashes.push((
                                    stream_name.clone(),
                                    part_id,
                                    format!("DECOMPRESS_FAILED: {}", e),
                                    compressed.len(),
                                    0,
                                    marker,
                                ));
                            }
                        }
                    }
                    Err(e) => {
                        ragc_hashes.push((
                            stream_name.clone(),
                            part_id,
                            format!("READ_FAILED: {}", e),
                            0,
                            0,
                            0,
                        ));
                    }
                }
                total_blocks += 1;
                if total_blocks % 100 == 0 {
                    print!("{}...", total_blocks);
                    std::io::Write::flush(&mut std::io::stdout()).ok();
                }
            }
        }
    }
    println!("Done! ({} blocks)", ragc_hashes.len());

    // Compare hashes - match by stream name and part ID
    println!("\n{}", "=".repeat(80));
    println!("COMPARISON RESULTS");
    println!("{}", "=".repeat(80));

    // Build lookup maps
    let mut cpp_map: HashMap<(String, usize), usize> = HashMap::new();
    for (idx, entry) in cpp_hashes.iter().enumerate() {
        cpp_map.insert((entry.0.clone(), entry.1), idx);
    }

    let mut ragc_map: HashMap<(String, usize), usize> = HashMap::new();
    for (idx, entry) in ragc_hashes.iter().enumerate() {
        ragc_map.insert((entry.0.clone(), entry.1), idx);
    }

    let mut identical = 0;
    let mut different_content = 0;
    let mut different_size = 0;
    let mut failures = 0;
    let mut only_in_cpp = 0;
    let mut only_in_ragc = 0;

    let mut differences = Vec::new();

    // Compare blocks that exist in both archives
    for (key, cpp_idx) in &cpp_map {
        if let Some(&ragc_idx) = ragc_map.get(key) {
            let cpp = &cpp_hashes[*cpp_idx];
            let ragc = &ragc_hashes[ragc_idx];

            let cpp_hash = &cpp.2;
            let ragc_hash = &ragc.2;

            if cpp_hash.starts_with("DECOMPRESS_FAILED")
                || cpp_hash.starts_with("READ_FAILED")
                || ragc_hash.starts_with("DECOMPRESS_FAILED")
                || ragc_hash.starts_with("READ_FAILED")
            {
                failures += 1;
                differences.push((*cpp_idx, cpp, ragc, "FAILURE"));
            } else if cpp_hash == ragc_hash {
                identical += 1;
            } else {
                if cpp.4 == ragc.4 {
                    different_content += 1;
                    differences.push((*cpp_idx, cpp, ragc, "CONTENT_DIFFERS"));
                } else {
                    different_size += 1;
                    differences.push((*cpp_idx, cpp, ragc, "SIZE_DIFFERS"));
                }
            }
        } else {
            only_in_cpp += 1;
        }
    }

    // Count blocks only in RAGC
    for key in ragc_map.keys() {
        if !cpp_map.contains_key(key) {
            only_in_ragc += 1;
        }
    }

    let total_compared = cpp_map.len().min(ragc_map.len());

    println!("\nTotal blocks in C++ AGC: {}", cpp_hashes.len());
    println!("Total blocks in RAGC: {}", ragc_hashes.len());
    println!("Blocks in both (compared): {}", total_compared);
    println!("Blocks only in C++ AGC: {}", only_in_cpp);
    println!("Blocks only in RAGC: {}", only_in_ragc);

    if total_compared > 0 {
        println!("\nComparison of common blocks:");
        println!(
            "  Identical decompressed content: {} ({:.1}%)",
            identical,
            100.0 * identical as f64 / total_compared as f64
        );
        println!(
            "  Different content (same size): {} ({:.1}%)",
            different_content,
            100.0 * different_content as f64 / total_compared as f64
        );
        println!(
            "  Different decompressed size: {} ({:.1}%)",
            different_size,
            100.0 * different_size as f64 / total_compared as f64
        );
        println!(
            "  Decompression failures: {} ({:.1}%)",
            failures,
            100.0 * failures as f64 / total_compared as f64
        );
    }

    if !differences.is_empty() {
        println!("\nFirst 20 differing blocks:");
        println!(
            "{:<6} {:<15} {:<5} {:<18} {:<18} {:<12} {:<12} {:<8}",
            "Index", "Stream", "Part", "CPP Hash", "RAGC Hash", "Comp Size", "Decomp Size", "Type"
        );
        println!("{}", "-".repeat(120));

        for (idx, cpp, ragc, diff_type) in differences.iter().take(20) {
            let comp_info = format!("{} vs {}", cpp.3, ragc.3);
            let decomp_info = format!("{} vs {}", cpp.4, ragc.4);

            println!(
                "{:<6} {:<15} {:<5} {:<18} {:<18} {:<12} {:<12} {:<8}",
                idx,
                cpp.0.chars().take(13).collect::<String>(),
                cpp.1,
                cpp.2.chars().take(16).collect::<String>(),
                ragc.2.chars().take(16).collect::<String>(),
                comp_info,
                decomp_info,
                diff_type
            );
        }
    }

    // Save detailed comparison
    let output_file = "/tmp/native_block_comparison.csv";
    let mut csv_output = String::new();
    csv_output.push_str("stream_name,part_id,cpp_hash,cpp_comp_size,cpp_decomp_size,cpp_marker,ragc_hash,ragc_comp_size,ragc_decomp_size,ragc_marker,status\n");

    for (key, cpp_idx) in &cpp_map {
        let cpp = &cpp_hashes[*cpp_idx];

        let (ragc_hash, ragc_comp, ragc_decomp, ragc_marker, status) =
            if let Some(&ragc_idx) = ragc_map.get(key) {
                let ragc = &ragc_hashes[ragc_idx];
                let status = if cpp.2.starts_with("DECOMPRESS_FAILED")
                    || cpp.2.starts_with("READ_FAILED")
                    || ragc.2.starts_with("DECOMPRESS_FAILED")
                    || ragc.2.starts_with("READ_FAILED")
                {
                    "FAILURE"
                } else if cpp.2 == ragc.2 {
                    "IDENTICAL"
                } else if cpp.4 == ragc.4 {
                    "CONTENT_DIFFERS"
                } else {
                    "SIZE_DIFFERS"
                };
                (ragc.2.clone(), ragc.3, ragc.4, ragc.5, status)
            } else {
                ("ONLY_IN_CPP".to_string(), 0, 0, 0, "ONLY_IN_CPP")
            };

        csv_output.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{}\n",
            cpp.0,
            cpp.1,
            cpp.2,
            cpp.3,
            cpp.4,
            cpp.5,
            ragc_hash,
            ragc_comp,
            ragc_decomp,
            ragc_marker,
            status
        ));
    }

    // Add blocks only in RAGC
    for (key, ragc_idx) in &ragc_map {
        if !cpp_map.contains_key(key) {
            let ragc = &ragc_hashes[*ragc_idx];
            csv_output.push_str(&format!(
                "{},{},ONLY_IN_RAGC,0,0,0,{},{},{},{},ONLY_IN_RAGC\n",
                ragc.0, ragc.1, ragc.2, ragc.3, ragc.4, ragc.5
            ));
        }
    }

    std::fs::write(output_file, csv_output)?;
    println!("\nFull comparison saved to {}", output_file);

    Ok(())
}
