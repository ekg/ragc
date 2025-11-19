// Test that Rust k-mer extraction produces identical results to C++ CKmer
// Second micro-replacement: k-mer canonicalization

#include <iostream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <algorithm>

// Need to include C++ AGC's CKmer class
#include "kmer.h"

using namespace std;

// Rust FFI
extern "C" {
    struct KmerArray {
        uint64_t* data;
        size_t len;
    };

    KmerArray ragc_extract_canonical_kmers(
        const uint8_t* contig_data,
        size_t contig_len,
        uint32_t k
    );

    void ragc_free_kmer_array(KmerArray array);

    uint64_t ragc_extract_kmer_at_position(
        const uint8_t* contig_data,
        size_t contig_len,
        uint32_t k,
        size_t pos
    );
}

// C++ AGC version (matches compress_contig loop logic)
vector<uint64_t> cpp_extract_canonical_kmers(const vector<uint8_t>& contig, uint32_t k) {
    vector<uint64_t> kmers;
    CKmer kmer(k, kmer_mode_t::canonical);

    for (auto x : contig) {
        if (x > 3) {
            // Non-ACGT base, reset k-mer
            kmer.Reset();
        } else {
            kmer.insert_canonical(x);

            if (kmer.is_full()) {
                kmers.push_back(kmer.data_canonical());
            }
        }
    }

    return kmers;
}

void test_kmer_extraction(const vector<uint8_t>& contig, uint32_t k, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  Input: contig_len=" << contig.size() << " k=" << k << endl;

    // Call C++ version
    auto cpp_kmers = cpp_extract_canonical_kmers(contig, k);

    // Call Rust version
    KmerArray rust_kmers = ragc_extract_canonical_kmers(contig.data(), contig.size(), k);

    // Compare
    bool match = (cpp_kmers.size() == rust_kmers.len);
    if (match) {
        for (size_t i = 0; i < cpp_kmers.size(); i++) {
            if (cpp_kmers[i] != rust_kmers.data[i]) {
                match = false;
                cout << "  ❌ Mismatch at position " << i << ":" << endl;
                cout << "     C++:  " << cpp_kmers[i] << endl;
                cout << "     Rust: " << rust_kmers.data[i] << endl;
                break;
            }
        }
    }

    cout << "  C++ k-mers: " << cpp_kmers.size() << " [";
    for (size_t i = 0; i < min(cpp_kmers.size(), size_t(5)); i++) {
        cout << cpp_kmers[i];
        if (i + 1 < min(cpp_kmers.size(), size_t(5))) cout << ",";
    }
    if (cpp_kmers.size() > 5) cout << "...";
    cout << "]" << endl;

    cout << "  Rust k-mers: " << rust_kmers.len << " [";
    for (size_t i = 0; i < min(rust_kmers.len, size_t(5)); i++) {
        cout << rust_kmers.data[i];
        if (i + 1 < min(rust_kmers.len, size_t(5))) cout << ",";
    }
    if (rust_kmers.len > 5) cout << "...";
    cout << "]" << endl;

    if (match) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    // Clean up Rust allocation
    ragc_free_kmer_array(rust_kmers);

    cout << endl;
}

void test_kmer_at_position(const vector<uint8_t>& contig, uint32_t k, size_t pos, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  Input: contig_len=" << contig.size() << " k=" << k << " pos=" << pos << endl;

    // Build C++ k-mer at position
    uint64_t cpp_kmer = ~0ull;
    if (pos + k <= contig.size()) {
        bool valid = true;
        CKmer kmer(k, kmer_mode_t::canonical);
        for (size_t i = 0; i < k; i++) {
            if (contig[pos + i] > 3) {
                valid = false;
                break;
            }
            kmer.insert_canonical(contig[pos + i]);
        }
        if (valid && kmer.is_full()) {
            cpp_kmer = kmer.data_canonical();
        }
    }

    // Call Rust version
    uint64_t rust_kmer = ragc_extract_kmer_at_position(contig.data(), contig.size(), k, pos);

    cout << "  C++:  " << cpp_kmer << endl;
    cout << "  Rust: " << rust_kmer << endl;

    if (cpp_kmer == rust_kmer) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ MISMATCH" << endl;
        exit(1);
    }

    cout << endl;
}

int main() {
    cout << "=== Testing Rust k-mer extraction vs C++ CKmer ===" << endl << endl;

    // Test 1: Simple ACGT sequence
    vector<uint8_t> seq1 = {0, 1, 2, 3, 0, 1, 2, 3}; // ACGTACGT
    test_kmer_extraction(seq1, 3, "Simple ACGT sequence (k=3)");
    test_kmer_extraction(seq1, 4, "Simple ACGT sequence (k=4)");
    test_kmer_extraction(seq1, 5, "Simple ACGT sequence (k=5)");

    // Test 2: Sequence with non-ACGT base (should reset)
    vector<uint8_t> seq2 = {0, 1, 2, 4, 0, 1, 2, 3}; // ACGTNACGT (N=4)
    test_kmer_extraction(seq2, 3, "Sequence with N (resets k-mer)");

    // Test 3: Multiple N bases
    vector<uint8_t> seq3 = {0, 1, 4, 2, 3, 4, 0, 1}; // ACNGTN...
    test_kmer_extraction(seq3, 3, "Multiple N bases");

    // Test 4: Short sequence (too short for k-mer)
    vector<uint8_t> seq4 = {0, 1}; // AC
    test_kmer_extraction(seq4, 5, "Too short for k-mer");

    // Test 5: Exact k-length
    vector<uint8_t> seq5 = {0, 1, 2}; // ACG
    test_kmer_extraction(seq5, 3, "Exact k-length (1 k-mer)");

    // Test 6: Longer realistic sequence
    vector<uint8_t> seq6(100);
    for (size_t i = 0; i < seq6.size(); i++) {
        seq6[i] = i % 4; // Repeating ACGT pattern
    }
    test_kmer_extraction(seq6, 21, "Longer sequence (k=21)");

    cout << "=== Testing k-mer extraction at specific positions ===" << endl << endl;

    // Test position-specific extraction
    vector<uint8_t> seq7 = {0, 1, 2, 3, 0, 1, 2, 3}; // ACGTACGT
    test_kmer_at_position(seq7, 3, 0, "Extract at position 0");
    test_kmer_at_position(seq7, 3, 3, "Extract at position 3");
    test_kmer_at_position(seq7, 3, 5, "Extract at position 5");

    // Test with N base
    vector<uint8_t> seq8 = {0, 4, 2, 3}; // ANGT
    test_kmer_at_position(seq8, 3, 0, "Extract at position with N");

    // Test out of bounds
    vector<uint8_t> seq9 = {0, 1, 2, 3};
    test_kmer_at_position(seq9, 3, 10, "Out of bounds");

    cout << "=== All tests passed! ===" << endl;
    cout << "✅ Rust k-mer extraction produces identical output to C++ CKmer" << endl;
    cout << "Ready for next micro-replacement" << endl;

    return 0;
}
