// Test that Rust splitter checking produces identical results to C++ AGC
// Third micro-replacement: splitter lookup decision

#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <unordered_set>

using namespace std;

// Rust FFI
extern "C" {
    bool ragc_is_splitter(
        uint64_t kmer_value,
        const uint64_t* splitters_ptr,
        size_t splitters_len
    );

    struct SplitterChecker {
        uint64_t* splitters;
        size_t len;
    };

    SplitterChecker ragc_create_splitter_checker(
        const uint64_t* splitters_ptr,
        size_t splitters_len
    );

    void ragc_free_splitter_checker(SplitterChecker checker);

    struct BoolArray {
        bool* data;
        size_t len;
    };

    BoolArray ragc_check_splitters_batch(
        const uint64_t* kmers_ptr,
        size_t kmers_len,
        const uint64_t* splitters_ptr,
        size_t splitters_len
    );

    void ragc_free_bool_array(BoolArray array);
}

// C++ AGC-style splitter checking (hash set)
bool cpp_is_splitter(uint64_t kmer_value, const unordered_set<uint64_t>& splitters) {
    return splitters.find(kmer_value) != splitters.end();
}

void test_single_checks(const vector<uint64_t>& splitters, const vector<uint64_t>& test_kmers, const char* desc) {
    cout << "Test: " << desc << endl;
    cout << "  Splitters: " << splitters.size() << " k-mers" << endl;
    cout << "  Testing: " << test_kmers.size() << " queries" << endl;

    // Build C++ hash set (matches C++ AGC's hs_splitters)
    unordered_set<uint64_t> cpp_splitters(splitters.begin(), splitters.end());

    // Test each k-mer
    size_t matches = 0;
    size_t mismatches = 0;

    for (uint64_t kmer : test_kmers) {
        bool cpp_result = cpp_is_splitter(kmer, cpp_splitters);
        bool rust_result = ragc_is_splitter(kmer, splitters.data(), splitters.size());

        if (cpp_result != rust_result) {
            cout << "  ❌ Mismatch for k-mer " << kmer << ":" << endl;
            cout << "     C++:  " << (cpp_result ? "true" : "false") << endl;
            cout << "     Rust: " << (rust_result ? "true" : "false") << endl;
            mismatches++;
        } else {
            matches++;
        }
    }

    cout << "  Results: " << matches << " matches, " << mismatches << " mismatches" << endl;

    if (mismatches == 0) {
        cout << "  ✅ MATCH" << endl;
    } else {
        cout << "  ❌ FAILED" << endl;
        exit(1);
    }

    cout << endl;
}

void test_batch_check(const vector<uint64_t>& splitters, const vector<uint64_t>& test_kmers, const char* desc) {
    cout << "Test: " << desc << endl;

    // Build C++ hash set
    unordered_set<uint64_t> cpp_splitters(splitters.begin(), splitters.end());

    // Get C++ results
    vector<bool> cpp_results;
    for (uint64_t kmer : test_kmers) {
        cpp_results.push_back(cpp_is_splitter(kmer, cpp_splitters));
    }

    // Get Rust results
    BoolArray rust_results = ragc_check_splitters_batch(
        test_kmers.data(), test_kmers.size(),
        splitters.data(), splitters.size()
    );

    // Compare
    bool match = true;
    for (size_t i = 0; i < test_kmers.size(); i++) {
        if (cpp_results[i] != rust_results.data[i]) {
            cout << "  ❌ Mismatch at index " << i << " (k-mer " << test_kmers[i] << ")" << endl;
            cout << "     C++:  " << (cpp_results[i] ? "true" : "false") << endl;
            cout << "     Rust: " << (rust_results.data[i] ? "true" : "false") << endl;
            match = false;
            break;
        }
    }

    if (match) {
        cout << "  ✅ MATCH (all " << test_kmers.size() << " checks)" << endl;
    } else {
        cout << "  ❌ FAILED" << endl;
        exit(1);
    }

    ragc_free_bool_array(rust_results);

    cout << endl;
}

int main() {
    cout << "=== Testing Rust splitter checking vs C++ AGC ===" << endl << endl;

    // Test 1: Small splitter set
    vector<uint64_t> splitters1 = {100, 200, 300, 400, 500};
    vector<uint64_t> queries1 = {50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550};
    test_single_checks(splitters1, queries1, "Small splitter set");

    // Test 2: Empty splitter set (no splitters)
    vector<uint64_t> splitters2 = {};
    vector<uint64_t> queries2 = {100, 200, 300};
    test_single_checks(splitters2, queries2, "Empty splitter set");

    // Test 3: Single splitter
    vector<uint64_t> splitters3 = {12345};
    vector<uint64_t> queries3 = {12345, 12346, 12344};
    test_single_checks(splitters3, queries3, "Single splitter");

    // Test 4: Larger realistic set (like actual splitters)
    vector<uint64_t> splitters4;
    for (uint64_t i = 0; i < 1000; i += 10) {
        splitters4.push_back(i * 1000000);
    }
    vector<uint64_t> queries4;
    for (uint64_t i = 0; i < 100; i++) {
        queries4.push_back(i * 1000000); // Mix of splitters and non-splitters
    }
    test_single_checks(splitters4, queries4, "Larger splitter set (100 splitters)");

    cout << "=== Testing batch checking ===" << endl << endl;

    // Test 5: Batch check small set
    test_batch_check(splitters1, queries1, "Batch check - small set");

    // Test 6: Batch check larger set
    test_batch_check(splitters4, queries4, "Batch check - large set");

    cout << "=== Testing splitter checker creation ===" << endl << endl;

    // Test 7: Verify sorted splitter array
    vector<uint64_t> unsorted = {500, 100, 300, 200, 400};
    SplitterChecker checker = ragc_create_splitter_checker(unsorted.data(), unsorted.size());

    cout << "Test: Splitter checker creation" << endl;
    cout << "  Input (unsorted): ";
    for (auto x : unsorted) cout << x << " ";
    cout << endl;

    cout << "  Output (should be sorted): ";
    for (size_t i = 0; i < checker.len; i++) {
        cout << checker.splitters[i] << " ";
    }
    cout << endl;

    // Verify it's sorted
    bool sorted = true;
    for (size_t i = 1; i < checker.len; i++) {
        if (checker.splitters[i-1] > checker.splitters[i]) {
            sorted = false;
            break;
        }
    }

    if (sorted) {
        cout << "  ✅ Correctly sorted" << endl;
    } else {
        cout << "  ❌ Not sorted" << endl;
        exit(1);
    }

    ragc_free_splitter_checker(checker);

    cout << endl;
    cout << "=== All tests passed! ===" << endl;
    cout << "✅ Rust splitter checking produces identical output to C++ AGC" << endl;
    cout << "Ready for next micro-replacement" << endl;

    return 0;
}
