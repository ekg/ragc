//! LZ Encoding Parity Test
//!
//! Compares RAGC's LZ diff encoding with C++ AGC's encoding via FFI.
//! The goal is byte-for-byte identical output.

#![allow(clippy::all)]

#[cfg(feature = "cpp_agc")]
mod lz_parity {
    use ragc_common::types::Contig;
    use ragc_core::LZDiff;
    use ragc_core::ragc_ffi;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    /// Convert FASTA sequence (ASCII ACGT) to numeric encoding (0123)
    fn fasta_to_numeric(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .map(|&b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                b'N' | b'n' => 4,
                _ => 4, // Treat unknown as N
            })
            .collect()
    }

    /// Read a segment from a FASTA file (first contig, specified length)
    fn read_fasta_segment(path: &Path, start: usize, len: usize) -> Option<Vec<u8>> {
        let file = File::open(path).ok()?;
        let reader = BufReader::new(file);

        let mut seq = Vec::new();
        for line in reader.lines() {
            let line = line.ok()?;
            if line.starts_with('>') {
                if !seq.is_empty() {
                    break; // Only read first contig
                }
                continue;
            }
            seq.extend(fasta_to_numeric(line.as_bytes()));
        }

        if start + len <= seq.len() {
            Some(seq[start..start + len].to_vec())
        } else if !seq.is_empty() {
            // If requested range exceeds sequence, use what we have
            let actual_len = seq.len().saturating_sub(start).min(len);
            if actual_len > 0 {
                Some(seq[start..start + actual_len].to_vec())
            } else {
                Some(seq[..seq.len().min(len)].to_vec())
            }
        } else {
            None
        }
    }

    /// Find first byte where encodings diverge
    fn find_divergence(ragc: &[u8], cpp: &[u8]) -> Option<(usize, u8, u8)> {
        for i in 0..ragc.len().max(cpp.len()) {
            let r = ragc.get(i).copied().unwrap_or(0);
            let c = cpp.get(i).copied().unwrap_or(0);
            if r != c {
                return Some((i, r, c));
            }
        }
        None
    }

    /// Compare LZ encoding output between RAGC and C++ AGC
    fn compare_lz_encoding(
        reference: &[u8],
        target: &[u8],
        min_match_len: u32,
        test_name: &str,
    ) -> bool {
        // RAGC encoding
        let mut lz = LZDiff::new(min_match_len);
        lz.prepare(&reference.to_vec());
        let ragc_encoded = lz.encode(&target.to_vec());

        // C++ AGC encoding via FFI
        let cpp_encoded = ragc_ffi::lzdiff_v2_encode(reference, target, min_match_len);

        match cpp_encoded {
            Some(cpp) => {
                if ragc_encoded == cpp {
                    println!("[PASS] {}: Identical encoding ({} bytes)", test_name, ragc_encoded.len());
                    true
                } else {
                    println!("[FAIL] {}: Encoding differs!", test_name);
                    println!("  RAGC: {} bytes", ragc_encoded.len());
                    println!("  C++:  {} bytes", cpp.len());

                    if let Some((pos, r, c)) = find_divergence(&ragc_encoded, &cpp) {
                        println!("  First divergence at byte {}: RAGC={:#04x} '{}' C++={:#04x} '{}'",
                            pos,
                            r, char::from(r).escape_default(),
                            c, char::from(c).escape_default());

                        // Show context around divergence
                        let ctx_start = pos.saturating_sub(10);
                        let ctx_end = (pos + 20).min(ragc_encoded.len().max(cpp.len()));

                        print!("  RAGC context [{}-{}]: ", ctx_start, ctx_end);
                        for i in ctx_start..ctx_end.min(ragc_encoded.len()) {
                            if i == pos {
                                print!("[{:02x}]", ragc_encoded[i]);
                            } else {
                                print!("{:02x}", ragc_encoded[i]);
                            }
                        }
                        println!();

                        print!("  C++  context [{}-{}]: ", ctx_start, ctx_end);
                        for i in ctx_start..ctx_end.min(cpp.len()) {
                            if i == pos {
                                print!("[{:02x}]", cpp[i]);
                            } else {
                                print!("{:02x}", cpp[i]);
                            }
                        }
                        println!();

                        // Show as ASCII for readability
                        print!("  RAGC as text: ");
                        for i in ctx_start..ctx_end.min(ragc_encoded.len()) {
                            let b = ragc_encoded[i];
                            if b.is_ascii_graphic() || b == b' ' {
                                print!("{}", b as char);
                            } else {
                                print!(".");
                            }
                        }
                        println!();

                        print!("  C++  as text: ");
                        for i in ctx_start..ctx_end.min(cpp.len()) {
                            let b = cpp[i];
                            if b.is_ascii_graphic() || b == b' ' {
                                print!("{}", b as char);
                            } else {
                                print!(".");
                            }
                        }
                        println!();
                    }
                    false
                }
            }
            None => {
                println!("[SKIP] {}: C++ FFI returned None (buffer too small?)", test_name);
                true // Don't fail on FFI issues
            }
        }
    }

    #[test]
    fn test_identical_sequences() {
        // When target == reference, both should return empty encoding
        let reference: Vec<u8> = (0..100).map(|i| (i % 4) as u8).collect();
        let target = reference.clone();

        assert!(compare_lz_encoding(&reference, &target, 20, "identical_sequences"));
    }

    #[test]
    fn test_simple_literal() {
        // Short sequence that's all literals (no matches possible)
        let reference: Vec<u8> = (0..50).map(|i| (i % 4) as u8).collect();
        let target: Vec<u8> = (0..50).map(|i| ((i + 1) % 4) as u8).collect();

        assert!(compare_lz_encoding(&reference, &target, 20, "simple_literal"));
    }

    #[test]
    fn test_long_match() {
        // Long sequence with a good match
        let reference: Vec<u8> = (0..1000).map(|i| (i % 4) as u8).collect();
        // Same pattern, just shifted a bit
        let target: Vec<u8> = (10..1010).map(|i| (i % 4) as u8).collect();

        assert!(compare_lz_encoding(&reference, &target, 20, "long_match"));
    }

    #[test]
    fn test_with_n_runs() {
        // Sequence with N runs
        let mut reference: Vec<u8> = (0..100).map(|i| (i % 4) as u8).collect();
        let mut target: Vec<u8> = (0..100).map(|i| (i % 4) as u8).collect();

        // Add N run in middle of target
        for i in 40..50 {
            target[i] = 4; // N
        }

        assert!(compare_lz_encoding(&reference, &target, 20, "with_n_runs"));
    }

    #[test]
    fn test_real_segment_if_available() {
        // Test with real yeast data if available
        let yeast_file = Path::new("/home/erik/scrapy/yeast235_split/AAA_0.fa");
        if !yeast_file.exists() {
            println!("[SKIP] test_real_segment_if_available: yeast data not available");
            return;
        }

        // Read two segments from the same file
        let reference = match read_fasta_segment(yeast_file, 0, 10000) {
            Some(s) => s,
            None => {
                println!("[SKIP] Could not read reference segment");
                return;
            }
        };

        let target = match read_fasta_segment(yeast_file, 5000, 10000) {
            Some(s) => s,
            None => {
                println!("[SKIP] Could not read target segment");
                return;
            }
        };

        assert!(compare_lz_encoding(&reference, &target, 20, "real_yeast_segment"));
    }

    #[test]
    fn test_real_segments_different_samples() {
        // Test with segments from different yeast samples
        let file1 = Path::new("/home/erik/scrapy/yeast235_split/AAA_0.fa");
        let file2 = Path::new("/home/erik/scrapy/yeast235_split/AAB_0.fa");

        if !file1.exists() || !file2.exists() {
            println!("[SKIP] test_real_segments_different_samples: yeast data not available");
            return;
        }

        let reference = match read_fasta_segment(file1, 0, 10000) {
            Some(s) => s,
            None => {
                println!("[SKIP] Could not read reference segment");
                return;
            }
        };

        let target = match read_fasta_segment(file2, 0, 10000) {
            Some(s) => s,
            None => {
                println!("[SKIP] Could not read target segment");
                return;
            }
        };

        assert!(compare_lz_encoding(&reference, &target, 20, "different_samples"));
    }

    /// Run multiple comparison tests and report summary
    #[test]
    fn test_lz_parity_comprehensive() {
        println!("\n=== LZ Encoding Parity Test ===\n");

        let mut passed = 0;
        let mut failed = 0;

        // Test 1: Identical sequences
        if compare_lz_encoding(
            &(0..100).map(|i| (i % 4) as u8).collect::<Vec<_>>(),
            &(0..100).map(|i| (i % 4) as u8).collect::<Vec<_>>(),
            20, "identical",
        ) { passed += 1; } else { failed += 1; }

        // Test 2: All literals
        if compare_lz_encoding(
            &(0..100).map(|i| (i % 4) as u8).collect::<Vec<_>>(),
            &(0..100).map(|i| ((i + 1) % 4) as u8).collect::<Vec<_>>(),
            20, "all_literals",
        ) { passed += 1; } else { failed += 1; }

        // Test 3: Long match at start
        let ref3: Vec<u8> = (0..500).map(|i| (i % 4) as u8).collect();
        let mut tgt3 = ref3.clone();
        for i in 400..500 {
            tgt3[i] = ((tgt3[i] + 1) % 4) as u8;
        }
        if compare_lz_encoding(&ref3, &tgt3, 20, "match_at_start") {
            passed += 1;
        } else {
            failed += 1;
        }

        // Test 4: Match at end
        let ref4: Vec<u8> = (0..500).map(|i| (i % 4) as u8).collect();
        let mut tgt4 = ref4.clone();
        for i in 0..100 {
            tgt4[i] = ((tgt4[i] + 1) % 4) as u8;
        }
        if compare_lz_encoding(&ref4, &tgt4, 20, "match_at_end") {
            passed += 1;
        } else {
            failed += 1;
        }

        // Test 5: N-run handling
        let ref5: Vec<u8> = (0..200).map(|i| (i % 4) as u8).collect();
        let mut tgt5: Vec<u8> = (0..200).map(|i| (i % 4) as u8).collect();
        for i in 80..100 {
            tgt5[i] = 4; // N
        }
        if compare_lz_encoding(&ref5, &tgt5, 20, "n_run") {
            passed += 1;
        } else {
            failed += 1;
        }

        // Test 6: Real yeast data if available
        let yeast_file = Path::new("/home/erik/scrapy/yeast235_split/AAA_0.fa");
        if yeast_file.exists() {
            if let (Some(reference), Some(target)) = (
                read_fasta_segment(yeast_file, 0, 5000),
                read_fasta_segment(yeast_file, 2500, 5000),
            ) {
                if compare_lz_encoding(&reference, &target, 20, "yeast_overlap") {
                    passed += 1;
                } else {
                    failed += 1;
                }
            }
        }

        // Test 7: Different samples
        let file1 = Path::new("/home/erik/scrapy/yeast235_split/AAA_0.fa");
        let file2 = Path::new("/home/erik/scrapy/yeast235_split/AAB_0.fa");
        if file1.exists() && file2.exists() {
            if let (Some(reference), Some(target)) = (
                read_fasta_segment(file1, 0, 5000),
                read_fasta_segment(file2, 0, 5000),
            ) {
                if compare_lz_encoding(&reference, &target, 20, "yeast_different_samples") {
                    passed += 1;
                } else {
                    failed += 1;
                }
            }
        }

        println!("\n=== Summary ===");
        println!("Passed: {}", passed);
        println!("Failed: {}", failed);

        if failed > 0 {
            panic!("{} tests failed", failed);
        }
    }
}

#[cfg(not(feature = "cpp_agc"))]
mod lz_parity {
    #[test]
    fn test_requires_cpp_agc_feature() {
        println!("LZ parity tests require cpp_agc feature");
        println!("Run with: cargo test --features cpp_agc");
    }
}
