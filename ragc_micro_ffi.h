// Rust micro-function FFI declarations
// All functions match C++ AGC logic exactly

#ifndef RAGC_MICRO_FFI_H
#define RAGC_MICRO_FFI_H

#include <cstdint>
#include <cstddef>

extern "C" {

// 1. Contig slicing (get_part)
struct ContigSlice {
    uint8_t* data;
    size_t len;
};
ContigSlice ragc_get_contig_part(const uint8_t* contig_data, size_t contig_len, uint64_t pos, uint64_t len);
void ragc_free_contig_slice(ContigSlice slice);

// 2. K-mer extraction
struct KmerArray {
    uint64_t* data;
    size_t len;
};
KmerArray ragc_extract_canonical_kmers(const uint8_t* contig_data, size_t contig_len, uint32_t k);
uint64_t ragc_extract_kmer_at_position(const uint8_t* contig_data, size_t contig_len, uint32_t k, size_t position);
void ragc_free_kmer_array(KmerArray array);

// 3. Splitter checking
bool ragc_is_splitter(uint64_t kmer_value, const uint64_t* splitters_ptr, size_t splitters_len);
struct SplitterCheckResult {
    bool* results;
    size_t len;
};
SplitterCheckResult ragc_check_splitters_batch(const uint64_t* kmer_values, size_t kmer_count, const uint64_t* splitters_ptr, size_t splitters_len);
void ragc_free_splitter_results(SplitterCheckResult result);

// 4. Segment boundary calculation
uint64_t ragc_calculate_split_position(uint64_t pos, uint32_t kmer_length);
uint64_t ragc_calculate_segment_end(uint64_t current_pos);
uint64_t ragc_calculate_segment_length(uint64_t segment_end, uint64_t split_pos);
struct SegmentBoundary {
    uint64_t segment_end;
    uint64_t new_split_pos;
    uint64_t segment_length;
};
SegmentBoundary ragc_calculate_segment_boundary(uint64_t current_pos, uint64_t current_split_pos, uint32_t kmer_length);

// 5. Base validation
bool ragc_is_valid_base(uint8_t base);
bool ragc_should_reset_kmer(uint8_t base);
struct BaseCounts {
    size_t n_valid;
    size_t n_invalid;
};
BaseCounts ragc_count_base_validity(const uint8_t* sequence, size_t length);

// 6. Reverse complement
uint8_t ragc_complement_base(uint8_t base);
void ragc_reverse_complement_inplace(uint8_t* sequence, size_t length);
struct Sequence {
    uint8_t* data;
    size_t len;
};
Sequence ragc_reverse_complement_copy(const uint8_t* src, size_t src_len);
void ragc_free_sequence(Sequence seq);

// 7. Segment split calculation
uint32_t ragc_calculate_seg2_start_pos(uint32_t left_size, uint32_t kmer_length);
uint32_t ragc_calculate_segment1_size(uint32_t seg2_start_pos, uint32_t kmer_length);
struct SegmentSplitInfo {
    uint32_t seg2_start_pos;
    uint32_t segment1_new_size;
};
SegmentSplitInfo ragc_calculate_segment_split(uint32_t left_size, uint32_t kmer_length);

// 8. K-mer pair operations
struct KmerPair {
    uint64_t first;
    uint64_t second;
};
KmerPair ragc_create_kmer_pair(uint64_t kmer1, uint64_t kmer2);
bool ragc_kmer_pair_equals(uint64_t p1_first, uint64_t p1_second, uint64_t p2_first, uint64_t p2_second);
bool ragc_is_empty_kmer(uint64_t kmer);
bool ragc_is_valid_kmer_pair(uint64_t first, uint64_t second);
KmerPair ragc_create_empty_kmer_pair();

} // extern "C"

#endif // RAGC_MICRO_FFI_H
