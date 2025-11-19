// Test reverse complement - Sixth micro-replacement

#include <iostream>
#include <vector>
#include <cstdint>

using namespace std;

extern "C" {
    uint8_t ragc_complement_base(uint8_t base);

    void ragc_reverse_complement_inplace(uint8_t* sequence, size_t length);

    struct Sequence {
        uint8_t* data;
        size_t len;
    };

    Sequence ragc_reverse_complement_copy(const uint8_t* src, size_t src_len);
    void ragc_free_sequence(Sequence seq);
}

// C++ AGC version (from agc_basic.cpp:257)
void cpp_reverse_complement_inplace(vector<uint8_t>& seq) {
    size_t n = seq.size();
    if (n == 0) return;

    size_t i = 0;
    size_t j = n - 1;

    while (i < j) {
        uint8_t x = (seq[j] < 4) ? (3 - seq[j]) : seq[j];
        uint8_t y = (seq[i] < 4) ? (3 - seq[i]) : seq[i];

        seq[i] = x;
        seq[j] = y;

        i++;
        j--;
    }

    // Handle middle element if odd length
    if (i == j) {
        seq[i] = (seq[i] < 4) ? (3 - seq[i]) : seq[i];
    }
}

// C++ AGC complement_base (from agc_basic.cpp)
uint8_t cpp_complement_base(uint8_t base) {
    return (base < 4) ? (3 - base) : base;
}

int main() {
    cout << "=== Testing Rust reverse complement vs C++ AGC ===\n\n";

    // Test 1: complement_base
    cout << "Test 1: complement_base\n";
    for (uint8_t base = 0; base < 10; base++) {
        uint8_t cpp_result = cpp_complement_base(base);
        uint8_t rust_result = ragc_complement_base(base);

        cout << "  base=" << (int)base
             << " C++=" << (int)cpp_result
             << " Rust=" << (int)rust_result;

        if (cpp_result == rust_result) {
            cout << " ✅\n";
        } else {
            cout << " ❌ MISMATCH\n";
            return 1;
        }
    }
    cout << "\n";

    // Test 2: reverse complement in-place (ACG -> CGT)
    cout << "Test 2: reverse_complement_inplace (ACG -> CGT)\n";
    vector<uint8_t> cpp_seq1 = {0, 1, 2}; // ACG
    vector<uint8_t> rust_seq1 = {0, 1, 2}; // ACG

    cpp_reverse_complement_inplace(cpp_seq1);
    ragc_reverse_complement_inplace(rust_seq1.data(), rust_seq1.size());

    cout << "  C++ result: ";
    for (auto b : cpp_seq1) cout << (int)b << " ";
    cout << "\n";

    cout << "  Rust result: ";
    for (auto b : rust_seq1) cout << (int)b << " ";
    cout << "\n";

    if (cpp_seq1 == rust_seq1) {
        cout << "  ✅ MATCH (both produce CGT: 1,2,3)\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 3: reverse complement in-place with N (ACNG -> CNGT)
    cout << "Test 3: reverse_complement_inplace with N (ACNG -> CNGT)\n";
    vector<uint8_t> cpp_seq2 = {0, 1, 4, 2}; // ACNG
    vector<uint8_t> rust_seq2 = {0, 1, 4, 2}; // ACNG

    cpp_reverse_complement_inplace(cpp_seq2);
    ragc_reverse_complement_inplace(rust_seq2.data(), rust_seq2.size());

    cout << "  C++ result: ";
    for (auto b : cpp_seq2) cout << (int)b << " ";
    cout << "\n";

    cout << "  Rust result: ";
    for (auto b : rust_seq2) cout << (int)b << " ";
    cout << "\n";

    if (cpp_seq2 == rust_seq2) {
        cout << "  ✅ MATCH (both produce CNGT: 1,4,2,3)\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 4: reverse complement copy (ACG -> CGT)
    cout << "Test 4: reverse_complement_copy (ACG -> CGT)\n";
    vector<uint8_t> src = {0, 1, 2}; // ACG
    vector<uint8_t> cpp_dest = src;
    cpp_reverse_complement_inplace(cpp_dest);

    Sequence rust_result = ragc_reverse_complement_copy(src.data(), src.size());

    cout << "  C++ result: ";
    for (auto b : cpp_dest) cout << (int)b << " ";
    cout << "\n";

    cout << "  Rust result: ";
    for (size_t i = 0; i < rust_result.len; i++) {
        cout << (int)rust_result.data[i] << " ";
    }
    cout << "\n";

    bool match = (rust_result.len == cpp_dest.size());
    if (match) {
        for (size_t i = 0; i < rust_result.len; i++) {
            if (rust_result.data[i] != cpp_dest[i]) {
                match = false;
                break;
            }
        }
    }

    if (match) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        ragc_free_sequence(rust_result);
        return 1;
    }

    ragc_free_sequence(rust_result);
    cout << "\n";

    // Test 5: even length palindrome (ACGT -> ACGT)
    cout << "Test 5: even length palindrome (ACGT -> ACGT)\n";
    vector<uint8_t> cpp_seq3 = {0, 1, 2, 3}; // ACGT
    vector<uint8_t> rust_seq3 = {0, 1, 2, 3}; // ACGT

    cpp_reverse_complement_inplace(cpp_seq3);
    ragc_reverse_complement_inplace(rust_seq3.data(), rust_seq3.size());

    cout << "  C++ result: ";
    for (auto b : cpp_seq3) cout << (int)b << " ";
    cout << "\n";

    cout << "  Rust result: ";
    for (auto b : rust_seq3) cout << (int)b << " ";
    cout << "\n";

    if (cpp_seq3 == rust_seq3) {
        cout << "  ✅ MATCH (both produce ACGT: 0,1,2,3)\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 6: odd length (single base A -> T)
    cout << "Test 6: single base (A -> T)\n";
    vector<uint8_t> cpp_seq4 = {0}; // A
    vector<uint8_t> rust_seq4 = {0}; // A

    cpp_reverse_complement_inplace(cpp_seq4);
    ragc_reverse_complement_inplace(rust_seq4.data(), rust_seq4.size());

    cout << "  C++ result: " << (int)cpp_seq4[0] << "\n";
    cout << "  Rust result: " << (int)rust_seq4[0] << "\n";

    if (cpp_seq4 == rust_seq4) {
        cout << "  ✅ MATCH (both produce T: 3)\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 7: empty sequence
    cout << "Test 7: empty sequence\n";
    Sequence empty_result = ragc_reverse_complement_copy(nullptr, 0);

    if (empty_result.len == 0 && empty_result.data == nullptr) {
        cout << "  ✅ Empty handled correctly\n";
    } else {
        cout << "  ❌ Empty not handled correctly\n";
        ragc_free_sequence(empty_result);
        return 1;
    }
    cout << "\n";

    cout << "=== All tests passed! ===\n";
    cout << "✅ Rust reverse complement matches C++ AGC exactly\n";

    return 0;
}
