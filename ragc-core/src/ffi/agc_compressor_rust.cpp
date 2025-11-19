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
    // Clear existing splitters
    v_candidate_kmers.clear();
    v_duplicated_kmers.clear();
    v_candidate_kmers_offset = 0;

    // Set new splitters (these are the singleton k-mers)
    v_candidate_kmers = singletons;  // All singletons
    v_duplicated_kmers = duplicates; // Duplicates for adaptive mode

    splitters_precomputed = true;

    cerr << "[RAGC FORK] Loaded " << singletons.size() << " singletons, "
         << duplicates.size() << " duplicates" << endl;
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

    // Set the splitters
    SetPrecomputedSplitters(splitters, singletons, duplicates);

    // Now call normal Create, which will use the pre-set splitters
    // TODO: Need to modify CAGCCompressor::Create to skip determine_splitters if already set
    cerr << "[RAGC FORK] TODO: Need to implement rest of Create() logic" << endl;

    return false;  // Not fully implemented yet
}
