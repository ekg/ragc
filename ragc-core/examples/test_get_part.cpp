// Test that Rust ragc_get_contig_part() produces identical results to C++ get_part()
// First micro-replacement: simplest segmentation helper function

#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <algorithm>

using namespace std;

// Rust FFI
extern "C" {
    struct ContigSlice {
        uint8_t* data;
        size_t len;
    };

    ContigSlice ragc_get_contig_part(
        const uint8_t* contig_data,
        size_t contig_len,
        uint64_t pos,
        uint64_t len
    );

    void ragc_free_contig_slice(ContigSlice slice);
}

// C++ AGC version (from agc_compressor.cpp:2101-2107)
vector<uint8_t> cpp_get_part(const vector<uint8_t>& contig, uint64_t pos, uint64_t len) {
    if (pos + len < contig.size())
        return vector<uint8_t>(contig.begin() + pos, contig.begin() + pos + len);
    else
        return vector<uint8_t>(contig.begin() + pos, contig.end());
}

void test_case(const vector<uint8_t>& contig, uint64_t pos, uint64_t len, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  Input: contig_len=" << contig.size() << " pos=" << pos << " len=" << len << endl;

    // Call C++ version
    auto cpp_result = cpp_get_part(contig, pos, len);

    // Call Rust version
    ContigSlice rust_result = ragc_get_contig_part(contig.data(), contig.size(), pos, len);

    // Compare
    bool match = (cpp_result.size() == rust_result.len);
    if (match) {
        match = (memcmp(cpp_result.data(), rust_result.data, rust_result.len) == 0);
    }

    cout << "  C++ result: len=" << cpp_result.size() << " data=[";
    for (size_t i = 0; i < min(cpp_result.size(), size_t(10)); i++) {
        cout << (int)cpp_result[i];
        if (i + 1 < min(cpp_result.size(), size_t(10))) cout << ",";
    }
    if (cpp_result.size() > 10) cout << "...";
    cout << "]" << endl;

    cout << "  Rust result: len=" << rust_result.len << " data=[";
    for (size_t i = 0; i < min(rust_result.len, size_t(10)); i++) {
        cout << (int)rust_result.data[i];
        if (i + 1 < min(rust_result.len, size_t(10))) cout << ",";
    }
    if (rust_result.len > 10) cout << "...";
    cout << "]" << endl;

    if (match) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    // Clean up Rust allocation
    ragc_free_contig_slice(rust_result);

    cout << endl;
}

int main() {
    cout << "=== Testing ragc_get_contig_part() vs C++ get_part() ===" << endl << endl;

    // Test 1: Extract middle portion (full length available)
    vector<uint8_t> contig1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    test_case(contig1, 2, 5, "Extract middle (full length)");

    // Test 2: Extract near end (exceeds contig length)
    vector<uint8_t> contig2 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    test_case(contig2, 7, 10, "Extract near end (exceeds length)");

    // Test 3: Extract exactly to end
    vector<uint8_t> contig3 = {0, 1, 2, 3, 4};
    test_case(contig3, 2, 3, "Extract exactly to end");

    // Test 4: Extract from start
    vector<uint8_t> contig4 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    test_case(contig4, 0, 5, "Extract from start");

    // Test 5: Extract single base
    vector<uint8_t> contig5 = {0, 1, 2, 3, 4};
    test_case(contig5, 2, 1, "Extract single base");

    // Test 6: Larger contig (more realistic)
    vector<uint8_t> contig6(1000);
    for (size_t i = 0; i < contig6.size(); i++) {
        contig6[i] = i % 256;
    }
    test_case(contig6, 100, 50, "Larger contig (middle)");
    test_case(contig6, 950, 100, "Larger contig (exceeds end)");

    cout << "=== All tests passed! ===" << endl;
    cout << "✅ Rust ragc_get_contig_part() produces identical output to C++ get_part()" << endl;
    cout << "Ready to integrate into compress_contig()" << endl;

    return 0;
}
