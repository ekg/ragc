use ragc_core::agc_compress_ffi::compress_with_cpp_agc;
fn main() -> anyhow::Result<()> {
    compress_with_cpp_agc(
        "/tmp/ffi_yeast.agc",
        &[("yeast_test".to_string(), "/tmp/yeast_test.fa".to_string())],
        21, 10000, 20, 50, false, false, 0, 1, 0.0
    )?;
    println!("Created /tmp/ffi_yeast.agc");
    Ok(())
}
