// Extended C++ AGC compressor with Rust integration
// Adds ability to use pre-computed splitters from Rust

#ifndef _AGC_COMPRESSOR_RUST_H
#define _AGC_COMPRESSOR_RUST_H

#include "../../../agc/src/core/agc_compressor.h"
#include <vector>
#include <cstdint>

// Rust splitter FFI
extern "C" {
    struct SplitterResult {
        uint64_t* splitters;
        size_t n_splitters;
        uint64_t* singletons;
        size_t n_singletons;
        uint64_t* duplicates;
        size_t n_duplicates;
    };

    SplitterResult* ragc_determine_splitters(
        const char* file_path,
        uint32_t k,
        uint32_t segment_size,
        uint32_t verbosity
    );

    void ragc_free_splitters(SplitterResult* result);
}

// Extended compressor with Rust splitter support
class CAGCCompressorRust : public CAGCCompressor {
public:
    // New Create() that uses Rust for splitter detection
    bool CreateWithRustSplitters(
        const std::string& _file_name,
        const uint32_t _pack_cardinality,
        const uint32_t _kmer_length,
        const std::string& reference_file_name,
        const uint32_t _segment_size,
        const uint32_t _min_match_len,
        const bool _concatenated_genomes,
        const bool _adaptive_compression,
        const uint32_t _verbosity,
        const uint32_t _no_threads,
        double _fallback_frac
    );

    // Set pre-computed splitters (called by CreateWithRustSplitters)
    void SetPrecomputedSplitters(
        const std::vector<uint64_t>& splitters,
        const std::vector<uint64_t>& singletons,
        const std::vector<uint64_t>& duplicates
    );

private:
    bool splitters_precomputed = false;
};

#endif // _AGC_COMPRESSOR_RUST_H
