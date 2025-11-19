// Test that Rust segment boundary calculation matches C++ AGC exactly
// Fourth micro-replacement: segment boundary arithmetic

#include <iostream>
#include <cstdint>
#include <algorithm>

using namespace std;

// Rust FFI
extern "C" {
    uint64_t ragc_calculate_split_position(uint64_t pos, uint32_t kmer_length);
    uint64_t ragc_calculate_segment_end(uint64_t pos);
    uint64_t ragc_calculate_segment_length(uint64_t split_pos, uint64_t segment_end);

    struct SegmentBoundary {
        uint64_t segment_end;
        uint64_t new_split_pos;
        uint64_t segment_length;
    };

    SegmentBoundary ragc_calculate_segment_boundary(
        uint64_t current_pos,
        uint64_t current_split_pos,
        uint32_t kmer_length
    );
}

// C++ AGC versions (from compress_contig logic)
uint64_t cpp_calculate_split_position(uint64_t pos, uint32_t kmer_length) {
    // Line 2044: split_pos = pos + 1 - kmer_length;
    if (pos + 1 < kmer_length) {
        return 0; // saturating subtraction
    }
    return pos + 1 - kmer_length;
}

uint64_t cpp_calculate_segment_end(uint64_t pos) {
    // Segment ends at pos + 1 (just after the splitter k-mer)
    return pos + 1;
}

uint64_t cpp_calculate_segment_length(uint64_t split_pos, uint64_t segment_end) {
    // Length is segment_end - split_pos
    if (segment_end < split_pos) {
        return 0; // saturating
    }
    return segment_end - split_pos;
}

void test_split_position(uint64_t pos, uint32_t k, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  pos=" << pos << " k=" << k << endl;

    uint64_t cpp_result = cpp_calculate_split_position(pos, k);
    uint64_t rust_result = ragc_calculate_split_position(pos, k);

    cout << "  C++:  " << cpp_result << endl;
    cout << "  Rust: " << rust_result << endl;

    if (cpp_result == rust_result) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    cout << endl;
}

void test_segment_end(uint64_t pos, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  pos=" << pos << endl;

    uint64_t cpp_result = cpp_calculate_segment_end(pos);
    uint64_t rust_result = ragc_calculate_segment_end(pos);

    cout << "  C++:  " << cpp_result << endl;
    cout << "  Rust: " << rust_result << endl;

    if (cpp_result == rust_result) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    cout << endl;
}

void test_segment_length(uint64_t split_pos, uint64_t segment_end, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  split_pos=" << split_pos << " segment_end=" << segment_end << endl;

    uint64_t cpp_result = cpp_calculate_segment_length(split_pos, segment_end);
    uint64_t rust_result = ragc_calculate_segment_length(split_pos, segment_end);

    cout << "  C++:  " << cpp_result << endl;
    cout << "  Rust: " << rust_result << endl;

    if (cpp_result == rust_result) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    cout << endl;
}

void test_complete_boundary(uint64_t current_pos, uint64_t current_split_pos, uint32_t k, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  current_pos=" << current_pos << " current_split_pos=" << current_split_pos << " k=" << k << endl;

    // Calculate with C++ functions
    uint64_t cpp_seg_end = cpp_calculate_segment_end(current_pos);
    uint64_t cpp_seg_len = cpp_calculate_segment_length(current_split_pos, cpp_seg_end);
    uint64_t cpp_new_split = cpp_calculate_split_position(current_pos, k);

    // Calculate with Rust
    SegmentBoundary rust_result = ragc_calculate_segment_boundary(current_pos, current_split_pos, k);

    cout << "  C++ - segment_end=" << cpp_seg_end << " new_split_pos=" << cpp_new_split << " segment_length=" << cpp_seg_len << endl;
    cout << "  Rust - segment_end=" << rust_result.segment_end << " new_split_pos=" << rust_result.new_split_pos << " segment_length=" << rust_result.segment_length << endl;

    bool match = (cpp_seg_end == rust_result.segment_end &&
                  cpp_new_split == rust_result.new_split_pos &&
                  cpp_seg_len == rust_result.segment_length);

    if (match) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    cout << endl;
}

int main() {
    cout << "=== Testing Rust segment boundary calculations vs C++ AGC ===" << endl << endl;

    cout << "=== Split position calculation ===" << endl << endl;

    test_split_position(100, 21, "Standard case (pos=100, k=21)");
    test_split_position(50, 3, "Small k-mer (pos=50, k=3)");
    test_split_position(10, 21, "pos < k (saturating subtraction)");
    test_split_position(0, 21, "pos=0");
    test_split_position(1000, 21, "Large position");

    cout << "=== Segment end calculation ===" << endl << endl;

    test_segment_end(100, "Standard case (pos=100)");
    test_segment_end(0, "pos=0");
    test_segment_end(999, "Large position");

    cout << "=== Segment length calculation ===" << endl << endl;

    test_segment_length(0, 101, "From start (split_pos=0, end=101)");
    test_segment_length(80, 101, "Mid-segment (split_pos=80, end=101)");
    test_segment_length(100, 50, "Invalid (end < start, saturating)");

    cout << "=== Complete boundary calculation ===" << endl << endl;

    test_complete_boundary(100, 0, 21, "First segment (splitter at 100)");
    test_complete_boundary(250, 80, 21, "Second segment (splitter at 250)");
    test_complete_boundary(500, 230, 21, "Third segment (splitter at 500)");
    test_complete_boundary(100, 0, 3, "Small k-mer (k=3)");

    cout << "=== Realistic scenario ===" << endl << endl;

    // Simulate compress_contig() logic for a few splitters
    cout << "Scenario: Contig with splitters at positions 100, 250, 500 (k=21)" << endl;
    cout << "Starting split_pos=0" << endl << endl;

    uint64_t split_pos = 0;
    uint32_t k = 21;

    for (uint64_t splitter_pos : {100, 250, 500}) {
        SegmentBoundary boundary = ragc_calculate_segment_boundary(splitter_pos, split_pos, k);

        cout << "  Splitter at pos=" << splitter_pos << endl;
        cout << "    Segment: [" << split_pos << ", " << boundary.segment_end << ") length=" << boundary.segment_length << endl;
        cout << "    Next split_pos: " << boundary.new_split_pos << endl;

        // Update for next iteration
        split_pos = boundary.new_split_pos;
        cout << endl;
    }

    cout << "=== All tests passed! ===" << endl;
    cout << "✅ Rust segment boundary calculations match C++ AGC exactly" << endl;
    cout << "Ready to integrate all micro-functions" << endl;

    return 0;
}
