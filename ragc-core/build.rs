fn main() {
    // Only compile C++ FFI when cpp_agc feature is enabled
    #[cfg(feature = "cpp_agc")]
    {
        // Path to C++ AGC source directories (local fork in ragc repo)
        let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set");
        let agc_common = format!("{}/../agc/src/common", manifest_dir);
        let agc_core = format!("{}/../agc/src/core", manifest_dir);

        // Also need HOME for 3rd party dependencies
        let home = std::env::var("HOME").expect("HOME not set");

        // Compile C++ AGC complete compression FFI (for byte-identical archives)
        {
            println!("cargo:rerun-if-changed=src/ffi/agc_compress.cpp");
            println!("cargo:rerun-if-changed=src/ffi/agc_index.cpp");
            println!("cargo:rerun-if-changed=src/ffi/splitters.cpp");
            println!("cargo:rerun-if-changed=src/ffi/agc_compressor_rust.cpp");
            println!("cargo:rerun-if-changed=src/ffi/cost.cpp");
            println!("cargo:rerun-if-changed={}/agc_compressor.cpp", agc_core);

            let mut build = cc::Build::new();
            build
                .cpp(true)
                .flag_if_supported("-std=c++20")
                .flag_if_supported("-O3")
                .flag_if_supported("-msse4.2")
                .flag_if_supported("-mavx2")
                .include(format!("{}/agc/src", home))
                .include(&agc_common)
                .include(&agc_core)
                .include(format!("{}/agc/3rd_party", home))
                .include(format!("{}/agc/3rd_party/raduls-inplace/Raduls", home))
                // FFI wrappers
                .file("src/ffi/agc_compress.cpp")
                .file("src/ffi/agc_index.cpp")
                .file("src/ffi/splitters.cpp")
                .file("src/ffi/agc_compressor_rust.cpp")
                .file("src/ffi/cost.cpp")
                // C++ AGC core compression
                .file(format!("{}/agc_compressor.cpp", agc_core))
                .file(format!("{}/agc_decompressor.cpp", agc_core))
                .file(format!("{}/genome_io.cpp", agc_core))
                // C++ AGC common
                .file(format!("{}/agc_basic.cpp", agc_common))
                .file(format!("{}/agc_decompressor_lib.cpp", agc_common))
                .file(format!("{}/archive.cpp", agc_common))
                .file(format!("{}/collection.cpp", agc_common))
                .file(format!("{}/collection_v1.cpp", agc_common))
                .file(format!("{}/collection_v2.cpp", agc_common))
                .file(format!("{}/collection_v3.cpp", agc_common))
                .file(format!("{}/lz_diff.cpp", agc_common))
                .file(format!("{}/segment.cpp", agc_common))
                .file(format!("{}/utils.cpp", agc_common));

            build.compile("agc_full");

            // Link against C++ AGC dependencies
            println!(
                "cargo:rustc-link-search=native={}/agc/3rd_party/zstd/lib",
                home
            );
            println!(
                "cargo:rustc-link-search=native={}/agc/3rd_party/libdeflate/build",
                home
            );
            println!(
                "cargo:rustc-link-search=native={}/agc/3rd_party/raduls-inplace/Raduls",
                home
            );

            println!("cargo:rustc-link-lib=static=zstd");
            println!("cargo:rustc-link-lib=static=deflate");
            println!("cargo:rustc-link-lib=static=raduls");
            println!("cargo:rustc-link-lib=z"); // System zlib for inflate/deflate

            // Link C++ stdlib
            #[cfg(target_os = "linux")]
            println!("cargo:rustc-link-lib=stdc++");
        }
    }
}
