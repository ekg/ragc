// Simple test of Rust splitter FFI (called directly from Rust)
use std::ffi::CString;

fn main() {
    let test_file = "/tmp/yeast_test.fa";
    let k = 21;
    let segment_size = 10000;
    let verbosity = 2;

    println!("Testing Rust FFI for splitter detection");
    println!("File: {}", test_file);
    println!("k: {}, segment_size: {}\n", k, segment_size);

    let c_path = CString::new(test_file).unwrap();

    let result = unsafe {
        ragc_core::splitters_ffi::ragc_determine_splitters(
            c_path.as_ptr(),
            k,
            segment_size,
            verbosity,
        )
    };

    if result.is_null() {
        eprintln!("ERROR: FFI returned NULL");
        return;
    }

    unsafe {
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
        ragc_core::splitters_ffi::ragc_free_splitters(result);

        println!("\n✓ FFI test completed successfully!");
    }
}
