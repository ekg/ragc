// C++ wrapper for Rust splitter detection
// Allows C++ AGC to use RAGC's splitter algorithm

#include <string>
#include <vector>
#include <cstdint>

extern "C" {

// Rust FFI result structure
struct SplitterResult {
    uint64_t* splitters;
    size_t n_splitters;
    uint64_t* singletons;
    size_t n_singletons;
    uint64_t* duplicates;
    size_t n_duplicates;
};

// Rust FFI functions
SplitterResult* ragc_determine_splitters(
    const char* file_path,
    uint32_t k,
    uint32_t segment_size,
    uint32_t verbosity
);

void ragc_free_splitters(SplitterResult* result);

} // extern "C"

// C++ wrapper for easier usage
class RustSplitterDetector {
public:
    RustSplitterDetector(
        const std::string& file_path,
        uint32_t k,
        uint32_t segment_size,
        uint32_t verbosity
    ) : result_(nullptr) {
        result_ = ragc_determine_splitters(
            file_path.c_str(),
            k,
            segment_size,
            verbosity
        );
    }

    ~RustSplitterDetector() {
        if (result_) {
            ragc_free_splitters(result_);
            result_ = nullptr;
        }
    }

    // No copy
    RustSplitterDetector(const RustSplitterDetector&) = delete;
    RustSplitterDetector& operator=(const RustSplitterDetector&) = delete;

    bool is_valid() const { return result_ != nullptr; }

    std::vector<uint64_t> get_splitters() const {
        if (!result_) return {};
        return std::vector<uint64_t>(
            result_->splitters,
            result_->splitters + result_->n_splitters
        );
    }

    std::vector<uint64_t> get_singletons() const {
        if (!result_) return {};
        return std::vector<uint64_t>(
            result_->singletons,
            result_->singletons + result_->n_singletons
        );
    }

    std::vector<uint64_t> get_duplicates() const {
        if (!result_) return {};
        return std::vector<uint64_t>(
            result_->duplicates,
            result_->duplicates + result_->n_duplicates
        );
    }

    size_t n_splitters() const { return result_ ? result_->n_splitters : 0; }
    size_t n_singletons() const { return result_ ? result_->n_singletons : 0; }
    size_t n_duplicates() const { return result_ ? result_->n_duplicates : 0; }

private:
    SplitterResult* result_;
};
