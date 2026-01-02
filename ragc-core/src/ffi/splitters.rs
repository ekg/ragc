// Rust FFI for splitter detection (called from C++ AGC)

use crate::{determine_splitters, GenomeIO};
use std::ffi::CStr;
use std::fs::File;
use std::os::raw::c_char;

/// C-compatible result structure for splitter detection
#[repr(C)]
pub struct SplitterResult {
    pub splitters: *mut u64,
    pub n_splitters: usize,
    pub singletons: *mut u64,
    pub n_singletons: usize,
    pub duplicates: *mut u64,
    pub n_duplicates: usize,
}

/// Free memory allocated by splitter detection
#[no_mangle]
pub extern "C" fn ragc_free_splitters(result: *mut SplitterResult) {
    if result.is_null() {
        return;
    }

    unsafe {
        let r = Box::from_raw(result);

        if !r.splitters.is_null() {
            drop(Vec::from_raw_parts(
                r.splitters,
                r.n_splitters,
                r.n_splitters,
            ));
        }
        if !r.singletons.is_null() {
            drop(Vec::from_raw_parts(
                r.singletons,
                r.n_singletons,
                r.n_singletons,
            ));
        }
        if !r.duplicates.is_null() {
            drop(Vec::from_raw_parts(
                r.duplicates,
                r.n_duplicates,
                r.n_duplicates,
            ));
        }
    }
}

/// Determine splitters from a FASTA file using RAGC's algorithm
///
/// This function is called from C++ AGC to replace its native splitter detection.
/// It produces identical results but uses Rust implementation.
///
/// # Safety
/// - `file_path` must be a valid null-terminated C string
/// - Caller must call `ragc_free_splitters()` to free the result
///
/// # Returns
/// Pointer to SplitterResult (NULL on error)
#[no_mangle]
pub extern "C" fn ragc_determine_splitters(
    file_path: *const c_char,
    k: u32,
    segment_size: u32,
    verbosity: u32,
) -> *mut SplitterResult {
    // Convert C string to Rust
    let path_str = unsafe {
        if file_path.is_null() {
            eprintln!("ragc_determine_splitters: NULL file path");
            return std::ptr::null_mut();
        }
        match CStr::from_ptr(file_path).to_str() {
            Ok(s) => s,
            Err(e) => {
                eprintln!("ragc_determine_splitters: Invalid UTF-8 in path: {}", e);
                return std::ptr::null_mut();
            }
        }
    };

    if verbosity > 0 {
        eprintln!("[RAGC] Gathering reference k-mers from: {}", path_str);
    }

    // Read contigs from file
    let file = match File::open(path_str) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("ragc_determine_splitters: Cannot open file: {}", e);
            return std::ptr::null_mut();
        }
    };

    let mut gio = GenomeIO::new(file);
    let mut contigs = Vec::new();

    while let Ok(Some((_name, contig))) = gio.read_contig_converted() {
        contigs.push(contig);
    }

    if verbosity > 0 {
        eprintln!("[RAGC] Loaded {} contigs", contigs.len());
    }

    // Run splitter detection
    if verbosity > 0 {
        eprintln!("[RAGC] Determination of splitters");
    }

    let (splitters_set, singletons_set, duplicates_set) =
        determine_splitters(&contigs, k as usize, segment_size as usize);

    if verbosity > 1 {
        eprintln!("[RAGC] No. of singletons: {}", singletons_set.len());
        eprintln!("[RAGC] No. of duplicates: {}", duplicates_set.len());
        eprintln!("[RAGC] No. of splitters: {}", splitters_set.len());
    }

    // Convert to C-compatible arrays
    let mut splitters: Vec<u64> = splitters_set.into_iter().collect();
    let mut singletons: Vec<u64> = singletons_set.into_iter().collect();
    let mut duplicates: Vec<u64> = duplicates_set.into_iter().collect();

    // Sort for deterministic output (C++ AGC expects sorted)
    splitters.sort_unstable();
    singletons.sort_unstable();
    duplicates.sort_unstable();

    let n_splitters = splitters.len();
    let n_singletons = singletons.len();
    let n_duplicates = duplicates.len();

    // Transfer ownership to C
    let splitters_ptr = splitters.as_mut_ptr();
    let singletons_ptr = singletons.as_mut_ptr();
    let duplicates_ptr = duplicates.as_mut_ptr();

    std::mem::forget(splitters);
    std::mem::forget(singletons);
    std::mem::forget(duplicates);

    let result = Box::new(SplitterResult {
        splitters: splitters_ptr,
        n_splitters,
        singletons: singletons_ptr,
        n_singletons,
        duplicates: duplicates_ptr,
        n_duplicates,
    });

    Box::into_raw(result)
}
