// Integration test - using multiple Rust micro-functions together
// Simulates the core loop of compress_contig()

#include <iostream>
#include <vector>
#include <cstdint>
#include "ragc_micro_ffi.h"

using namespace std;

int main() {
    cout << "=== Micro-function Integration Test ===\n\n";

    // Test sequence: ACGTACGTACGTACGTACGT (20 bases)
    vector<uint8_t> contig = {
        0, 1, 2, 3,  // ACGT
        0, 1, 2, 3,  // ACGT
        0, 1, 2, 3,  // ACGT
        0, 1, 2, 3,  // ACGT
        0, 1, 2, 3   // ACGT
    };

    cout << "Test contig: 20 bases (ACGTACGTACGTACGTACGT)\n";
    cout << "K-mer length: 21\n\n";

    // Test 1: Base validation
    cout << "Test 1: Base validation\n";
    BaseCounts counts = ragc_count_base_validity(contig.data(), contig.size());
    cout << "  Valid bases: " << counts.n_valid << " (expect 20)\n";
    cout << "  Invalid bases: " << counts.n_invalid << " (expect 0)\n";
    if (counts.n_valid == 20 && counts.n_invalid == 0) {
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    cout << "\n";

    // Test 2: Extract k-mers
    cout << "Test 2: K-mer extraction (k=4)\n";
    KmerArray kmers = ragc_extract_canonical_kmers(contig.data(), contig.size(), 4);
    cout << "  Extracted " << kmers.len << " k-mers (expect 17)\n";
    if (kmers.len == 17) {  // 20 - 4 + 1 = 17
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    ragc_free_kmer_array(kmers);
    cout << "\n";

    // Test 3: Contig slicing (get_part)
    cout << "Test 3: Contig slicing\n";
    ContigSlice slice = ragc_get_contig_part(contig.data(), contig.size(), 5, 10);
    cout << "  Slice from pos 5, length 10: got " << slice.len << " bases\n";
    if (slice.len == 10) {
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    ragc_free_contig_slice(slice);
    cout << "\n";

    // Test 4: Segment boundary calculation
    cout << "Test 4: Segment boundary calculation\n";
    uint64_t current_pos = 100;
    uint64_t current_split_pos = 50;
    uint32_t kmer_length = 21;
    SegmentBoundary boundary = ragc_calculate_segment_boundary(current_pos, current_split_pos, kmer_length);
    cout << "  Current pos: " << current_pos << "\n";
    cout << "  Current split pos: " << current_split_pos << "\n";
    cout << "  Segment end: " << boundary.segment_end << " (expect 101)\n";
    cout << "  New split pos: " << boundary.new_split_pos << " (expect 80)\n";
    cout << "  Segment length: " << boundary.segment_length << " (expect 51)\n";
    if (boundary.segment_end == 101 && boundary.new_split_pos == 80 && boundary.segment_length == 51) {
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    cout << "\n";

    // Test 5: Reverse complement
    cout << "Test 5: Reverse complement\n";
    vector<uint8_t> test_seq = {0, 1, 2};  // ACG
    Sequence rc = ragc_reverse_complement_copy(test_seq.data(), test_seq.size());
    cout << "  Input: ACG (0,1,2)\n";
    cout << "  Output: ";
    for (size_t i = 0; i < rc.len; i++) {
        cout << (int)rc.data[i];
        if (i + 1 < rc.len) cout << ",";
    }
    cout << " (expect CGT: 1,2,3)\n";
    bool rc_match = (rc.len == 3 && rc.data[0] == 1 && rc.data[1] == 2 && rc.data[2] == 3);
    if (rc_match) {
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    ragc_free_sequence(rc);
    cout << "\n";

    // Test 6: K-mer pair operations
    cout << "Test 6: K-mer pair operations\n";
    KmerPair pair1 = ragc_create_kmer_pair(200, 100);  // Should swap
    cout << "  Input: (200, 100)\n";
    cout << "  Output: (" << pair1.first << ", " << pair1.second << ")\n";
    cout << "  Expected: (100, 200)\n";
    if (pair1.first == 100 && pair1.second == 200) {
        cout << "  ✅ PASS - correctly ordered\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    cout << "\n";

    // Test 7: Segment split calculation
    cout << "Test 7: Segment split calculation\n";
    SegmentSplitInfo split = ragc_calculate_segment_split(250, 21);
    cout << "  Input: left_size=250, k=21\n";
    cout << "  seg2_start_pos: " << split.seg2_start_pos << " (expect 240)\n";
    cout << "  segment1_new_size: " << split.segment1_new_size << " (expect 261)\n";
    if (split.seg2_start_pos == 240 && split.segment1_new_size == 261) {
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    cout << "\n";

    // Test 8: Splitter checking
    cout << "Test 8: Splitter checking\n";
    vector<uint64_t> splitters = {100, 200, 300, 400};
    bool is_splitter_200 = ragc_is_splitter(200, splitters.data(), splitters.size());
    bool is_splitter_150 = ragc_is_splitter(150, splitters.data(), splitters.size());
    cout << "  is_splitter(200): " << is_splitter_200 << " (expect true)\n";
    cout << "  is_splitter(150): " << is_splitter_150 << " (expect false)\n";
    if (is_splitter_200 == true && is_splitter_150 == false) {
        cout << "  ✅ PASS\n";
    } else {
        cout << "  ❌ FAIL\n";
        return 1;
    }
    cout << "\n";

    cout << "=== All integration tests passed! ===\n";
    cout << "✅ All 8 micro-functions work correctly together\n";

    return 0;
}
