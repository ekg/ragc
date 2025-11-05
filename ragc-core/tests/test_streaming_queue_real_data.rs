// Test streaming queue compressor with first 10 yeast samples (sequential multi-file mode)

use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[test]
#[ignore] // Requires external yeast data files in ~/scrapy/yeast235_samples
fn test_streaming_queue_yeast10() -> anyhow::Result<()> {
    let samples_dir = std::env::var("HOME").unwrap() + "/scrapy/yeast235_samples";

    // First 10 yeast samples (sequential processing for multi-file mode)
    let sample_files = vec![
        "AAA#0.fa", "AAB#0.fa", "AAC#0.fa", "AAR#0.fa", "ABA#0.fa", "ABH#0.fa", "ACA#0.fa",
        "ACH#0.fa", "ADE#0.fa", "ADI#0.fa",
    ];

    let config = StreamingQueueConfig {
        k: 21,
        segment_size: 10000,
        min_match_len: 20,
        queue_capacity: 2 * 1024 * 1024 * 1024,
        num_threads: 4,
        verbosity: 2,
        compression_level: 17,
        adaptive_mode: false, // Test without adaptive mode first
    };

    // Detect splitters from first sample (matching batch mode behavior)
    let first_sample_path = format!("{}/{}", samples_dir, sample_files[0]);
    println!("Detecting splitters from first sample: {}", sample_files[0]);
    let (splitters, _, _) = ragc_core::determine_splitters_streaming(
        Path::new(&first_sample_path),
        config.k as usize,
        config.segment_size as usize,
    )?;
    println!("Found {} splitters", splitters.len());

    let mut compressor = StreamingQueueCompressor::with_splitters(
        "/tmp/test_yeast10_streaming_queue.agc",
        config,
        splitters,
    )?;

    // Process each sample file sequentially (multi-file mode)
    for sample_file in sample_files {
        let file_path = format!("{}/{}", samples_dir, sample_file);
        println!("Processing {}...", sample_file);

        let file = File::open(&file_path)?;
        let reader = BufReader::new(file);

        let mut current_sample = String::new();
        let mut current_contig = String::new();
        let mut current_data = Vec::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                // Flush previous contig if any
                if !current_data.is_empty() {
                    compressor.push(
                        current_sample.clone(),
                        current_contig.clone(),
                        current_data.clone(),
                    )?;
                    current_data.clear();
                }

                // Parse header: >AAA#0#chr1
                let parts: Vec<&str> = line[1..].split('#').collect();
                if parts.len() >= 2 {
                    current_sample = format!("{}#{}", parts[0], parts[1]);
                    current_contig = if parts.len() >= 3 {
                        parts[2].to_string()
                    } else {
                        "unknown".to_string()
                    };
                }
            } else {
                // Add sequence data
                for byte in line.bytes() {
                    let base = match byte {
                        b'A' | b'a' => 0u8,
                        b'C' | b'c' => 1u8,
                        b'G' | b'g' => 2u8,
                        b'T' | b't' => 3u8,
                        _ => 0u8, // Treat N as A
                    };
                    current_data.push(base);
                }
            }
        }

        // Flush last contig from this sample
        if !current_data.is_empty() {
            compressor.push(current_sample, current_contig, current_data)?;
        }
    }

    compressor.finalize()?;

    println!("Archive created successfully!");
    Ok(())
}
