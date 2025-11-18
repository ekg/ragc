fn main() {
    #[cfg(feature = "ffi_cost")]
    {
        println!("cargo:rerun-if-changed=src/ffi/cost.cpp");
        let mut build = cc::Build::new();
        build.cpp(true).flag_if_supported("-std=c++17").file("src/ffi/cost.cpp");
        build.compile("agc_cost");
        // Link C++ stdlib on Linux
        #[cfg(target_os = "linux")]
        println!("cargo:rustc-link-lib=stdc++");
    }
}

