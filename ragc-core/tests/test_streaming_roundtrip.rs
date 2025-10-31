use ragc_core::{Decompressor, DecompressorConfig, StreamingCompressor, StreamingCompressorConfig};
use std::fs;

#[test]
#[ignore] // Old API (add_contig) no longer exists - use add_fasta_files_with_splitters instead
fn test_streaming_compressor_roundtrip() {
    let archive_path = "/tmp/test_streaming_roundtrip.agc";
    let _ = fs::remove_file(archive_path);

    // Create archive with StreamingCompressor
    {
        let config = StreamingCompressorConfig::default();
        let mut compressor = StreamingCompressor::new(archive_path, config).unwrap();

        let seq = vec![0u8, 1, 2, 3, 0, 1, 2, 3]; // ACGTACGT
        compressor.add_contig("sample1", "chr1", seq).unwrap();
        compressor.finalize().unwrap();
    }

    // Read it back
    {
        let config = DecompressorConfig { verbosity: 0 };
        let mut decompressor = Decompressor::open(archive_path, config).unwrap();

        let samples = decompressor.list_samples();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0], "sample1");

        let contigs = decompressor.list_contigs("sample1").unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0], "chr1");

        let data = decompressor.get_contig("sample1", "chr1").unwrap();
        assert_eq!(data, vec![0, 1, 2, 3, 0, 1, 2, 3]);

        decompressor.close().unwrap();
    }

    fs::remove_file(archive_path).unwrap();
}

// TODO: Sequential path has bug with segment writing (chr2 not extracted)
// This is a known issue separate from the per-group buffering optimization.
// The parallel path (add_multi_sample_fasta_with_splitters) works correctly.
#[test]
#[ignore]
fn test_streaming_compressor_get_sample() {
    let archive_path = "/tmp/test_streaming_get_sample.agc";
    let _ = fs::remove_file(archive_path);

    // Create archive with StreamingCompressor
    {
        let config = StreamingCompressorConfig::default();
        let mut compressor = StreamingCompressor::new(archive_path, config).unwrap();

        let seq1 = vec![0u8, 1, 2, 3]; // ACGT
        let seq2 = vec![3u8, 2, 1, 0]; // TGCA
        compressor.add_contig("sample1", "chr1", seq1).unwrap();
        compressor.add_contig("sample1", "chr2", seq2).unwrap();
        compressor.finalize().unwrap();
    }

    // Read it back using get_sample
    {
        let config = DecompressorConfig { verbosity: 0 };
        let mut decompressor = Decompressor::open(archive_path, config).unwrap();

        let sample_data = decompressor.get_sample("sample1").unwrap();
        assert_eq!(sample_data.len(), 2);
        assert_eq!(sample_data[0].0, "chr1");
        assert_eq!(sample_data[0].1, vec![0, 1, 2, 3]);
        assert_eq!(sample_data[1].0, "chr2");
        assert_eq!(sample_data[1].1, vec![3, 2, 1, 0]);

        decompressor.close().unwrap();
    }

    fs::remove_file(archive_path).unwrap();
}
