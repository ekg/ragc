// Test streaming queue compressor with first 10 yeast samples

use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
use std::collections::HashSet;
use flate2::read::GzDecoder;
use std::io::{BufRead, BufReader};
use std::fs::File;

#[test]
fn test_streaming_queue_yeast10() -> anyhow::Result<()> {
    let input_path = std::env::var("HOME").unwrap() + "/scrapy/yeast235.fa.gz";
    
    let config = StreamingQueueConfig {
        k: 21,
        segment_size: 10000,
        min_match_len: 20,
        queue_capacity: 2 * 1024 * 1024 * 1024,
        num_threads: 4,
        verbosity: 1,
        compression_level: 17,
    };
    
    let mut compressor = StreamingQueueCompressor::with_splitters(
        "/tmp/test_yeast10_streaming_queue.agc",
        config,
        HashSet::new(),
    )?;
    
    // Read and push first 10 samples from yeast235
    let file = File::open(&input_path)?;
    let gz = GzDecoder::new(file);
    let reader = BufReader::new(gz);
    
    let mut current_sample = String::new();
    let mut current_contig = String::new();
    let mut current_data = Vec::new();
    let mut sample_count = 0;
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Flush previous contig if any
            if !current_data.is_empty() && sample_count < 10 {
                compressor.push(current_sample.clone(), current_contig.clone(), current_data.clone())?;
                current_data.clear();
            }
            
            // Parse header: >AAA#0#chr1
            let parts: Vec<&str> = line[1..].split('#').collect();
            if parts.len() >= 2 {
                let new_sample = format!("{}#{}", parts[0], parts[1]);
                if new_sample != current_sample {
                    if !current_sample.is_empty() {
                        sample_count += 1;
                        if sample_count >= 10 {
                            break;
                        }
                    }
                    current_sample = new_sample;
                }
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
    
    // Flush last contig
    if !current_data.is_empty() && sample_count < 10 {
        compressor.push(current_sample, current_contig, current_data)?;
    }
    
    compressor.finalize()?;
    
    println!("Archive created successfully!");
    Ok(())
}
