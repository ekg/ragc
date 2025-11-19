// Test base validation - Fifth micro-replacement

#include <iostream>
#include <vector>
#include <cstdint>

using namespace std;

extern "C" {
    bool ragc_is_valid_base(uint8_t base);
    bool ragc_should_reset_kmer(uint8_t base);

    struct BaseCounts {
        size_t n_valid;
        size_t n_invalid;
    };

    BaseCounts ragc_count_base_validity(const uint8_t* sequence, size_t length);

    struct PositionArray {
        size_t* data;
        size_t len;
    };

    PositionArray ragc_find_invalid_base_positions(const uint8_t* sequence, size_t length);
    void ragc_free_position_array(PositionArray array);
}

// C++ AGC version (from compress_contig line 2025)
bool cpp_is_valid_base(uint8_t x) {
    return !(x >> 2);  // Equivalent to x <= 3
}

bool cpp_should_reset_kmer(uint8_t x) {
    return x >> 2;  // Equivalent to x > 3
}

int main() {
    cout << "=== Testing Rust base validation vs C++ AGC ===" << endl << endl;

    // Test individual bases
    cout << "Single base validation:" << endl;
    for (uint8_t base = 0; base < 10; base++) {
        bool cpp_valid = cpp_is_valid_base(base);
        bool rust_valid = ragc_is_valid_base(base);
        bool cpp_reset = cpp_should_reset_kmer(base);
        bool rust_reset = ragc_should_reset_kmer(base);

        cout << "  base=" << (int)base
             << " C++ valid=" << cpp_valid << " Rust valid=" << rust_valid
             << " C++ reset=" << cpp_reset << " Rust reset=" << rust_reset;

        if (cpp_valid == rust_valid && cpp_reset == rust_reset) {
            cout << " ✅" << endl;
        } else {
            cout << " ❌ MISMATCH" << endl;
            return 1;
        }
    }

    cout << endl;

    // Test sequence counting
    vector<uint8_t> seq1 = {0, 1, 2, 3, 4, 0, 1, 5}; // ACGTNACX
    BaseCounts counts = ragc_count_base_validity(seq1.data(), seq1.size());

    // Count manually with C++ logic
    size_t cpp_valid = 0, cpp_invalid = 0;
    for (auto base : seq1) {
        if (cpp_is_valid_base(base)) cpp_valid++;
        else cpp_invalid++;
    }

    cout << "Sequence counting (ACGTNACX):" << endl;
    cout << "  C++ - valid=" << cpp_valid << " invalid=" << cpp_invalid << endl;
    cout << "  Rust - valid=" << counts.n_valid << " invalid=" << counts.n_invalid << endl;

    if (cpp_valid == counts.n_valid && cpp_invalid == counts.n_invalid) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        return 1;
    }

    cout << endl;

    // Test position finding
    vector<uint8_t> seq2 = {0, 1, 4, 2, 3, 5, 0}; // ACNGTXA
    PositionArray positions = ragc_find_invalid_base_positions(seq2.data(), seq2.size());

    cout << "Invalid base positions (ACNGTXA):" << endl;
    cout << "  Expected: 2, 5" << endl;
    cout << "  Rust found: ";
    for (size_t i = 0; i < positions.len; i++) {
        cout << positions.data[i];
        if (i + 1 < positions.len) cout << ", ";
    }
    cout << endl;

    bool pos_match = (positions.len == 2 && positions.data[0] == 2 && positions.data[1] == 5);
    if (pos_match) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        return 1;
    }

    ragc_free_position_array(positions);

    cout << endl << "=== All tests passed! ===" << endl;
    cout << "✅ Rust base validation matches C++ AGC exactly" << endl;

    return 0;
}
