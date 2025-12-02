// Test that Rust splitter FFI works correctly
use std::ffi::CString;

extern "C" {
    fn ragc_determine_splitters(
        file_path: *const std::os::raw::c_char,
        k: u32,
        segment_size: u32,
        verbosity: u32,
    ) -> *mut SplitterResult;

    fn ragc_free_splitters(result: *mut SplitterResult);
}

#[repr(C)]
struct SplitterResult {
    splitters: *mut u64,
    n_splitters: usize,
    singletons: *mut u64,
    n_singletons: usize,
    duplicates: *mut u64,
    n_duplicates: usize,
}

fn main() {
    let test_file = "/tmp/yeast_test.fa";
    let k = 21;
    let segment_size = 10000;
    let verbosity = 2;

    println!("Testing Rust FFI for splitter detection");
    println!("File: {}", test_file);
    println!("k: {}, segment_size: {}\n", k, segment_size);

    let c_path = CString::new(test_file).unwrap();

    unsafe {
        let result = ragc_determine_splitters(c_path.as_ptr(), k, segment_size, verbosity);

        if result.is_null() {
            eprintln!("ERROR: FFI returned NULL");
            return;
        }

        let r = &*result;

        println!("\n=== FFI Results ===");
        println!("Singletons: {}", r.n_singletons);
        println!("Duplicates: {}", r.n_duplicates);
        println!("Splitters:  {}", r.n_splitters);

        // Show splitters
        if r.n_splitters > 0 {
            println!("\nSplitters:");
            let splitters = std::slice::from_raw_parts(r.splitters, r.n_splitters);
            for (i, &sp) in splitters.iter().enumerate() {
                println!("  {}: {:016x}", i, sp);
            }
        }

        // Clean up
        ragc_free_splitters(result);

        println!("\n✓ FFI test completed successfully!");
    }
}
