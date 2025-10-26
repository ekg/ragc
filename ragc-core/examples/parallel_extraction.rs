// Parallel Sample Extraction Example
//
// This example demonstrates:
// - Thread-safe concurrent access to AGC archives
// - Using clone_for_thread() to create independent readers
// - Extracting multiple samples in parallel

use anyhow::Result;
use ragc_core::{Decompressor, DecompressorConfig};
use std::sync::Arc;
use std::thread;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <archive.agc> [num_threads]", args[0]);
        eprintln!();
        eprintln!("Examples:");
        eprintln!("  {} data.agc       # Use default threads (num_cpus)", args[0]);
        eprintln!("  {} data.agc 4     # Use 4 threads", args[0]);
        std::process::exit(1);
    }

    let archive_path = &args[1];
    let num_threads: usize = args
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(num_cpus::get);

    println!("Opening archive: {}", archive_path);
    println!("Using {} threads\n", num_threads);

    // Open the main decompressor
    let decompressor = Decompressor::open(archive_path, DecompressorConfig::default())?;

    // Get all samples
    let samples = decompressor.list_samples();
    println!("Found {} samples", samples.len());

    if samples.is_empty() {
        println!("No samples to extract!");
        return Ok(());
    }

    // Wrap samples in Arc for sharing across threads
    let samples_arc = Arc::new(samples);

    // Spawn worker threads
    let start_time = std::time::Instant::now();
    let mut handles = Vec::new();

    for thread_id in 0..num_threads {
        let samples_clone = Arc::clone(&samples_arc);
        let thread_decompressor = decompressor.clone_for_thread()?;

        let handle = thread::spawn(move || {
            extract_worker(thread_id, thread_decompressor, samples_clone, num_threads)
        });

        handles.push(handle);
    }

    // Wait for all threads to complete and collect results
    let mut total_contigs = 0;
    let mut total_bases: u64 = 0;

    for handle in handles {
        let (contigs, bases) = handle.join().expect("Thread panicked")?;
        total_contigs += contigs;
        total_bases += bases;
    }

    let elapsed = start_time.elapsed();

    // Display results
    println!("\n--- Results ---");
    println!("Samples processed: {}", samples_arc.len());
    println!("Total contigs: {}", total_contigs);
    println!("Total bases: {} bp", total_bases);
    println!("Time elapsed: {:.2}s", elapsed.as_secs_f64());
    println!(
        "Throughput: {:.2} MB/s",
        (total_bases as f64 / 1_000_000.0) / elapsed.as_secs_f64()
    );

    Ok(())
}

/// Worker function that processes samples assigned to this thread
fn extract_worker(
    thread_id: usize,
    mut decompressor: Decompressor,
    samples: Arc<Vec<String>>,
    num_threads: usize,
) -> Result<(usize, u64)> {
    let mut contigs_extracted = 0;
    let mut bases_extracted: u64 = 0;

    // Process samples assigned to this thread (round-robin)
    for (i, sample_name) in samples.iter().enumerate() {
        // Skip samples not assigned to this thread
        if i % num_threads != thread_id {
            continue;
        }

        // Extract the sample
        match decompressor.get_sample(sample_name) {
            Ok(contigs) => {
                let sample_bases: usize = contigs.iter().map(|(_, seq)| seq.len()).sum();

                println!(
                    "[Thread {}] Extracted {}: {} contigs, {} bp",
                    thread_id,
                    sample_name,
                    contigs.len(),
                    sample_bases
                );

                contigs_extracted += contigs.len();
                bases_extracted += sample_bases as u64;
            }
            Err(e) => {
                eprintln!("[Thread {}] Error extracting {}: {}", thread_id, sample_name, e);
            }
        }
    }

    // Close this thread's decompressor
    decompressor.close()?;

    Ok((contigs_extracted, bases_extracted))
}
