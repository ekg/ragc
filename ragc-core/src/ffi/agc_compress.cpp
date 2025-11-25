// FFI bridge to C++ AGC complete compression pipeline
// Provides a simple C interface to CAGCCompressor for byte-identical archives

#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <cstring>

// C++ AGC headers
#include "../../../agc/src/core/agc_compressor.h"
#include "../../../agc/src/common/lz_diff.h"
#include "agc_compressor_rust.h"

extern "C" {

// Compression parameters (matches CLI args)
struct AgcCompressParams {
    const char* out_archive_name;
    uint32_t pack_cardinality;      // batch size (default 50)
    uint32_t kmer_length;            // k-mer length
    const char* reference_file;      // first file for splitter detection
    uint32_t segment_size;           // segment size
    uint32_t min_match_length;       // LZ min match
    bool concatenated_genomes;       // PanSN mode
    bool adaptive_compression;       // adaptive mode
    uint32_t verbosity;              // verbosity level
    uint32_t no_threads;             // thread count
    double fallback_frac;            // fallback fraction (default 0.0)
};

// Sample file specification
struct AgcSampleFile {
    const char* sample_name;
    const char* file_path;
};

// Main compression function - wraps CAGCCompressor::Create/AddSampleFiles/Close
bool agc_compress_create(
    const AgcCompressParams* params,
    const AgcSampleFile* samples,
    size_t n_samples
) {
    if (!params || !samples || n_samples == 0) {
        return false;
    }

    try {
        CAGCCompressor compressor;

        // Step 1: Create compressor and detect splitters from reference
        bool ok = compressor.Create(
            std::string(params->out_archive_name),
            params->pack_cardinality,
            params->kmer_length,
            std::string(params->reference_file),
            params->segment_size,
            params->min_match_length,
            params->concatenated_genomes,
            params->adaptive_compression,
            params->verbosity,
            params->no_threads,
            params->fallback_frac
        );

        if (!ok) {
            return false;
        }

        // Step 2: Add all sample files
        std::vector<std::pair<std::string, std::string>> v_sample_files;
        v_sample_files.reserve(n_samples);

        for (size_t i = 0; i < n_samples; ++i) {
            v_sample_files.emplace_back(
                std::string(samples[i].sample_name),
                std::string(samples[i].file_path)
            );
        }

        ok = compressor.AddSampleFiles(v_sample_files, params->no_threads);
        if (!ok) {
            return false;
        }

        // Step 3: Finalize and write archive
        ok = compressor.Close(params->no_threads);
        return ok;

    } catch (const std::exception& e) {
        // Log error if verbosity is enabled
        if (params->verbosity > 0) {
            fprintf(stderr, "C++ AGC compression error: %s\n", e.what());
        }
        return false;
    } catch (...) {
        if (params->verbosity > 0) {
            fprintf(stderr, "C++ AGC compression error: unknown exception\n");
        }
        return false;
    }
}

// Compression function using Rust splitter detection
bool agc_compress_create_with_rust_splitters(
    const AgcCompressParams* params,
    const AgcSampleFile* samples,
    size_t n_samples
) {
    if (!params || !samples || n_samples == 0) {
        return false;
    }

    try {
        CAGCCompressorRust compressor;

        // Step 1: Create compressor using Rust splitter detection
        bool ok = compressor.CreateWithRustSplitters(
            std::string(params->out_archive_name),
            params->pack_cardinality,
            params->kmer_length,
            std::string(params->reference_file),
            params->segment_size,
            params->min_match_length,
            params->concatenated_genomes,
            params->adaptive_compression,
            params->verbosity,
            params->no_threads,
            params->fallback_frac
        );

        if (!ok) {
            return false;
        }

        // Step 2: Add all sample files
        std::vector<std::pair<std::string, std::string>> v_sample_files;
        v_sample_files.reserve(n_samples);

        for (size_t i = 0; i < n_samples; ++i) {
            v_sample_files.emplace_back(
                std::string(samples[i].sample_name),
                std::string(samples[i].file_path)
            );
        }

        ok = compressor.AddSampleFiles(v_sample_files, params->no_threads);
        if (!ok) {
            return false;
        }

        // Step 3: Finalize and write archive
        ok = compressor.Close(params->no_threads);
        return ok;

    } catch (const std::exception& e) {
        if (params->verbosity > 0) {
            fprintf(stderr, "[RAGC FORK] C++ AGC compression error: %s\n", e.what());
        }
        return false;
    } catch (...) {
        if (params->verbosity > 0) {
            fprintf(stderr, "[RAGC FORK] C++ AGC compression error: unknown exception\n");
        }
        return false;
    }
}

// Compute encoding cost estimate using actual CLZDiff_V2::Estimate
// This uses the EXACT same algorithm as C++ AGC's find_cand_segment_with_one_splitter
uint32_t agc_lzdiff_v2_estimate(
    const uint8_t* ref, size_t ref_len,
    const uint8_t* text, size_t text_len,
    uint32_t min_match_len,
    uint32_t bound
) {
    if (!ref || !text || ref_len == 0 || text_len == 0) {
        return UINT32_MAX;
    }

    // Create CLZDiff_V2 with min_match_len
    CLZDiff_V2 lz(min_match_len);

    // Build reference contig
    contig_t reference_contig(ref, ref + ref_len);

    // Prepare() sets reference and prepares for encoding
    lz.Prepare(reference_contig);

    // AssureIndex() ensures the hash index is built (called before Estimate)
    lz.AssureIndex();

    // Build text contig
    contig_t text_contig(text, text + text_len);

    // Call the actual Estimate function
    size_t result = lz.Estimate(text_contig, bound);

    return (uint32_t)result;
}

} // extern "C"
