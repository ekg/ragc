// Extended C++ AGC compressor with Rust integration
// MINIMAL VERSION - just test that fork works

#include "agc_compressor_rust.h"
#include <iostream>

using namespace std;

void CAGCCompressorRust::SetPrecomputedSplitters(
    const std::vector<uint64_t>& splitters,
    const std::vector<uint64_t>& singletons,
    const std::vector<uint64_t>& duplicates
) {
    // Clear existing data
    v_candidate_kmers.clear();
    v_duplicated_kmers.clear();
    v_candidate_kmers_offset = 0;
    hs_splitters.clear();

    // Set singletons (for adaptive mode and splitter finding if needed)
    v_candidate_kmers = singletons;
    v_duplicated_kmers = duplicates;

    // Populate hs_splitters with actual splitters
    for (const auto& splitter : splitters) {
        hs_splitters.insert(splitter);
    }

    // Populate bloom filter for fast splitter lookup
    // Using same sizing as C++ AGC (line 554 of agc_compressor.cpp)
    bloom_splitters.resize((uint64_t)(hs_splitters.size() / 0.25));
    bloom_splitters.insert(hs_splitters.begin(), hs_splitters.end());

    splitters_precomputed = true;

    cerr << "[RAGC FORK] Loaded " << singletons.size() << " singletons, "
         << duplicates.size() << " duplicates, "
         << splitters.size() << " splitters" << endl;
}

bool CAGCCompressorRust::CreateWithRustSplitters(
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
) {
    cerr << "[RAGC FORK] CreateWithRustSplitters - computing splitters with Rust FFI" << endl;

    // Get splitters from Rust
    SplitterResult* rust_result = ragc_determine_splitters(
        reference_file_name.c_str(),
        _kmer_length,
        _segment_size,
        _verbosity
    );

    if (!rust_result) {
        cerr << "[RAGC FORK] Error: Rust splitter detection failed" << endl;
        return false;
    }

    // Convert to C++ vectors
    std::vector<uint64_t> splitters(
        rust_result->splitters,
        rust_result->splitters + rust_result->n_splitters
    );

    std::vector<uint64_t> singletons(
        rust_result->singletons,
        rust_result->singletons + rust_result->n_singletons
    );

    std::vector<uint64_t> duplicates(
        rust_result->duplicates,
        rust_result->duplicates + rust_result->n_duplicates
    );

    ragc_free_splitters(rust_result);

    cerr << "[RAGC FORK] Rust found: " << singletons.size() << " singletons, "
         << splitters.size() << " splitters" << endl;

    // Set the splitters BEFORE calling Create()
    SetPrecomputedSplitters(splitters, singletons, duplicates);

    // Now call normal Create(), which will detect pre-set splitters and skip determine_splitters()
    cerr << "[RAGC FORK] Calling Create() with pre-computed splitters..." << endl;

    bool result = CAGCCompressor::Create(
        _file_name,
        _pack_cardinality,
        _kmer_length,
        reference_file_name,  // Still needed for file path but won't be read
        _segment_size,
        _min_match_len,
        _concatenated_genomes,
        _adaptive_compression,
        _verbosity,
        _no_threads,
        _fallback_frac
    );

    if (result) {
        cerr << "[RAGC FORK] Create() succeeded with Rust splitters!" << endl;
    } else {
        cerr << "[RAGC FORK] Create() failed" << endl;
    }

    return result;
}
