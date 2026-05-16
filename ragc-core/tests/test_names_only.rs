/// Integration tests for names-only contig loading.
///
/// Requires AGC_TEST_FILE env var pointing to an AGC archive.
/// Run with: AGC_TEST_FILE=/path/to/file.agc cargo test -p ragc-core --test test_names_only -- --nocapture
use ragc_core::{Decompressor, DecompressorConfig};

fn open_agc() -> Option<(Decompressor, String)> {
    let path = std::env::var("AGC_TEST_FILE").ok()?;
    let config = DecompressorConfig { verbosity: 0 };
    let dec = Decompressor::open(&path, config).ok()?;
    Some((dec, path))
}

/// Test 1: list_contigs_names_only() returns the same names as list_contigs()
#[test]
fn names_only_matches_full_load() {
    let Some((mut dec, _path)) = open_agc() else {
        eprintln!("Skipping: AGC_TEST_FILE not set");
        return;
    };

    let samples = dec.list_samples();
    assert!(!samples.is_empty(), "AGC file should have samples");

    // Pick first 3 samples to keep test fast
    for sample in samples.iter().take(3) {
        let names_only = dec.list_contigs_names_only(sample).unwrap();
        assert!(!names_only.is_empty(), "Sample {sample} should have contigs");

        // Open a fresh decompressor for full load comparison
        let config = DecompressorConfig { verbosity: 0 };
        let mut dec_full = Decompressor::open(&_path, config).unwrap();
        let full = dec_full.list_contigs(sample).unwrap();

        assert_eq!(
            names_only, full,
            "list_contigs_names_only() and list_contigs() differ for sample {sample}"
        );
        eprintln!("  OK: {sample} — {} contigs match", names_only.len());
    }
}

/// Test 2: get_contig_range() works after names-only loading
#[test]
fn get_contig_range_after_names_only() {
    let Some((mut dec, _path)) = open_agc() else {
        eprintln!("Skipping: AGC_TEST_FILE not set");
        return;
    };

    let samples = dec.list_samples();
    let sample = &samples[0];

    // Do names-only load first
    let contigs = dec.list_contigs_names_only(sample).unwrap();
    let contig = &contigs[0];

    // Now fetch a small range — this should trigger full detail reload internally
    let seq = dec.get_contig_range(sample, contig, 0, 100).unwrap();
    assert_eq!(seq.len(), 100, "Should get exactly 100 bases");
    assert!(seq.iter().all(|&b| b <= 3), "All bases should be 0-3 encoded");

    // Verify against a fresh full-load decompressor
    let config = DecompressorConfig { verbosity: 0 };
    let mut dec_full = Decompressor::open(&_path, config).unwrap();
    let seq_full = dec_full.get_contig_range(sample, contig, 0, 100).unwrap();
    assert_eq!(seq, seq_full, "Sequence should match between names-only and full load paths");

    eprintln!("  OK: get_contig_range({sample}, {contig}, 0, 100) — {} bases match", seq.len());
}

/// Test 3: get_contig_length() works after names-only loading
#[test]
fn get_contig_length_after_names_only() {
    let Some((mut dec, _path)) = open_agc() else {
        eprintln!("Skipping: AGC_TEST_FILE not set");
        return;
    };

    let samples = dec.list_samples();
    let sample = &samples[0];

    // Do names-only load first
    let contigs = dec.list_contigs_names_only(sample).unwrap();
    let contig = &contigs[0];

    // get_contig_length should trigger full detail reload internally
    let len = dec.get_contig_length(sample, contig).unwrap();
    assert!(len > 0, "Contig length should be positive");

    // Verify against a fresh full-load decompressor
    let config = DecompressorConfig { verbosity: 0 };
    let mut dec_full = Decompressor::open(&_path, config).unwrap();
    let len_full = dec_full.get_contig_length(sample, contig).unwrap();
    assert_eq!(len, len_full, "Length should match between names-only and full load paths");

    eprintln!("  OK: get_contig_length({sample}, {contig}) = {len}");
}

/// Test 4: two consecutive misses on `list_contigs_names_only` must not panic.
///
/// Regression for the codex-flagged P1: without rewinding `samples_loaded`
/// between calls, a second miss would feed a stale offset into
/// `deserialize_contig_names` and panic with an out-of-bounds index into
/// `sample_desc`. Expected behavior is a clean `Sample not found` error
/// on every miss.
#[test]
fn repeated_missing_samples_do_not_panic() {
    let Some((mut dec, _path)) = open_agc() else {
        eprintln!("Skipping: AGC_TEST_FILE not set");
        return;
    };

    for label in ["__missing_one__", "__missing_two__", "__missing_three__"] {
        let err = dec
            .list_contigs_names_only(label)
            .expect_err("missing sample should return an error, not data");
        assert!(
            err.to_string().contains("Sample not found"),
            "unexpected error for {label}: {err}"
        );
    }
    eprintln!("  OK: three consecutive misses returned errors without panicking");
}
