// Rust FFI bindings to C++ AGC complete compression

use std::ffi::CString;
use std::os::raw::c_char;
use std::path::Path;

#[repr(C)]
struct AgcCompressParams {
    out_archive_name: *const c_char,
    pack_cardinality: u32,
    kmer_length: u32,
    reference_file: *const c_char,
    segment_size: u32,
    min_match_length: u32,
    concatenated_genomes: bool,
    adaptive_compression: bool,
    verbosity: u32,
    no_threads: u32,
    fallback_frac: f64,
}

#[repr(C)]
struct AgcSampleFile {
    sample_name: *const c_char,
    file_path: *const c_char,
}

extern "C" {
    fn agc_compress_create(
        params: *const AgcCompressParams,
        samples: *const AgcSampleFile,
        n_samples: usize,
    ) -> bool;

    fn agc_compress_create_with_rust_splitters(
        params: *const AgcCompressParams,
        samples: *const AgcSampleFile,
        n_samples: usize,
    ) -> bool;
}

/// Compress FASTA files using C++ AGC's complete compression pipeline
///
/// This produces byte-identical archives to C++ AGC by delegating the entire
/// compression process to the C++ implementation.
pub fn compress_with_cpp_agc(
    out_archive: impl AsRef<Path>,
    sample_files: &[(String, String)], // (sample_name, file_path)
    kmer_length: u32,
    segment_size: u32,
    min_match_length: u32,
    pack_cardinality: u32,
    concatenated_genomes: bool,
    adaptive_compression: bool,
    verbosity: u32,
    no_threads: u32,
    fallback_frac: f64,
) -> anyhow::Result<()> {
    if sample_files.is_empty() {
        return Err(anyhow::anyhow!("No sample files provided"));
    }

    // Convert strings to CStrings (must keep alive during FFI call)
    let c_out_archive = CString::new(out_archive.as_ref().to_string_lossy().as_ref())?;
    let c_ref_file = CString::new(sample_files[0].1.as_str())?; // First file is reference

    let c_sample_names: Vec<CString> = sample_files
        .iter()
        .map(|(name, _)| CString::new(name.as_str()).unwrap())
        .collect();

    let c_file_paths: Vec<CString> = sample_files
        .iter()
        .map(|(_, path)| CString::new(path.as_str()).unwrap())
        .collect();

    let c_samples: Vec<AgcSampleFile> = c_sample_names
        .iter()
        .zip(c_file_paths.iter())
        .map(|(name, path)| AgcSampleFile {
            sample_name: name.as_ptr(),
            file_path: path.as_ptr(),
        })
        .collect();

    let params = AgcCompressParams {
        out_archive_name: c_out_archive.as_ptr(),
        pack_cardinality,
        kmer_length,
        reference_file: c_ref_file.as_ptr(),
        segment_size,
        min_match_length,
        concatenated_genomes,
        adaptive_compression,
        verbosity,
        no_threads,
        fallback_frac,
    };

    unsafe {
        if agc_compress_create(&params, c_samples.as_ptr(), c_samples.len()) {
            Ok(())
        } else {
            Err(anyhow::anyhow!("C++ AGC compression failed"))
        }
    }
}

/// Compress FASTA files using C++ AGC with Rust splitter detection
///
/// This uses Rust's splitter detection algorithm (verified to match C++ AGC exactly),
/// then delegates the rest of compression to C++ AGC.
///
/// # Purpose
/// First step in progressive Rust replacement - splitters computed in Rust,
/// everything else in C++. Produces byte-identical archives to native C++ AGC.
pub fn compress_with_rust_splitters(
    out_archive: impl AsRef<Path>,
    sample_files: &[(String, String)], // (sample_name, file_path)
    kmer_length: u32,
    segment_size: u32,
    min_match_length: u32,
    pack_cardinality: u32,
    concatenated_genomes: bool,
    adaptive_compression: bool,
    verbosity: u32,
    no_threads: u32,
    fallback_frac: f64,
) -> anyhow::Result<()> {
    if sample_files.is_empty() {
        return Err(anyhow::anyhow!("No sample files provided"));
    }

    // Convert strings to CStrings (must keep alive during FFI call)
    let c_out_archive = CString::new(out_archive.as_ref().to_string_lossy().as_ref())?;
    let c_ref_file = CString::new(sample_files[0].1.as_str())?; // First file is reference

    let c_sample_names: Vec<CString> = sample_files
        .iter()
        .map(|(name, _)| CString::new(name.as_str()).unwrap())
        .collect();

    let c_file_paths: Vec<CString> = sample_files
        .iter()
        .map(|(_, path)| CString::new(path.as_str()).unwrap())
        .collect();

    let c_samples: Vec<AgcSampleFile> = c_sample_names
        .iter()
        .zip(c_file_paths.iter())
        .map(|(name, path)| AgcSampleFile {
            sample_name: name.as_ptr(),
            file_path: path.as_ptr(),
        })
        .collect();

    let params = AgcCompressParams {
        out_archive_name: c_out_archive.as_ptr(),
        pack_cardinality,
        kmer_length,
        reference_file: c_ref_file.as_ptr(),
        segment_size,
        min_match_length,
        concatenated_genomes,
        adaptive_compression,
        verbosity,
        no_threads,
        fallback_frac,
    };

    unsafe {
        if agc_compress_create_with_rust_splitters(&params, c_samples.as_ptr(), c_samples.len()) {
            Ok(())
        } else {
            Err(anyhow::anyhow!("[RAGC FORK] C++ AGC compression with Rust splitters failed"))
        }
    }
}
