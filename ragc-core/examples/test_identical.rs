// Test for byte-identical archive comparison
use ragc_core::agc_compress_ffi::compress_with_cpp_agc;

fn main() -> anyhow::Result<()> {
    // Same sample name as C++ AGC would extract from filename
    let sample_files = vec![("test_identical".to_string(), "/tmp/test_identical.fa".to_string())];

    compress_with_cpp_agc(
        "/tmp/ffi.agc",
        &sample_files,
        21,      // kmer_length
        10000,   // segment_size
        20,      // min_match_length
        50,      // pack_cardinality (C++ AGC default)
        false,   // concatenated_genomes
        false,   // adaptive_compression
        0,       // verbosity (quiet)
        1,       // no_threads
        0.0,     // fallback_frac
    )?;

    println!("FFI archive created: /tmp/ffi.agc");
    Ok(())
}
