// Test k-mer pair operations - Eighth micro-replacement

#include <iostream>
#include <cstdint>
#include <algorithm>

using namespace std;

extern "C" {
    struct KmerPair {
        uint64_t first;
        uint64_t second;
    };

    KmerPair ragc_create_kmer_pair(uint64_t kmer1, uint64_t kmer2);
    bool ragc_kmer_pair_equals(uint64_t p1_first, uint64_t p1_second, uint64_t p2_first, uint64_t p2_second);
    bool ragc_is_empty_kmer(uint64_t kmer);
    bool ragc_is_valid_kmer_pair(uint64_t first, uint64_t second);
    KmerPair ragc_create_empty_kmer_pair();
}

// C++ AGC version - using std::minmax
pair<uint64_t, uint64_t> cpp_create_kmer_pair(uint64_t kmer1, uint64_t kmer2) {
    return minmax(kmer1, kmer2);
}

bool cpp_is_empty_kmer(uint64_t kmer) {
    return kmer == ~0ull;
}

bool cpp_is_valid_kmer_pair(uint64_t first, uint64_t second) {
    return first != ~0ull && second != ~0ull;
}

int main() {
    cout << "=== Testing Rust k-mer pair vs C++ AGC ===\n\n";

    // Test 1: Create k-mer pair (forward order)
    cout << "Test 1: Create k-mer pair (100, 200)\n";
    auto cpp_pair = cpp_create_kmer_pair(100, 200);
    KmerPair rust_pair = ragc_create_kmer_pair(100, 200);

    cout << "  C++ pair: (" << cpp_pair.first << ", " << cpp_pair.second << ")\n";
    cout << "  Rust pair: (" << rust_pair.first << ", " << rust_pair.second << ")\n";

    if (cpp_pair.first == rust_pair.first && cpp_pair.second == rust_pair.second) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 2: Create k-mer pair (reverse order - should swap)
    cout << "Test 2: Create k-mer pair (200, 100) - should swap\n";
    cpp_pair = cpp_create_kmer_pair(200, 100);
    rust_pair = ragc_create_kmer_pair(200, 100);

    cout << "  C++ pair: (" << cpp_pair.first << ", " << cpp_pair.second << ")\n";
    cout << "  Rust pair: (" << rust_pair.first << ", " << rust_pair.second << ")\n";

    if (cpp_pair.first == rust_pair.first && cpp_pair.second == rust_pair.second &&
        rust_pair.first == 100 && rust_pair.second == 200) {
        cout << "  ✅ MATCH (correctly swapped to (100, 200))\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 3: Equal values
    cout << "Test 3: Create k-mer pair (150, 150)\n";
    cpp_pair = cpp_create_kmer_pair(150, 150);
    rust_pair = ragc_create_kmer_pair(150, 150);

    cout << "  C++ pair: (" << cpp_pair.first << ", " << cpp_pair.second << ")\n";
    cout << "  Rust pair: (" << rust_pair.first << ", " << rust_pair.second << ")\n";

    if (cpp_pair.first == rust_pair.first && cpp_pair.second == rust_pair.second) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 4: Empty k-mer check
    cout << "Test 4: Empty k-mer check (~0ull)\n";
    bool cpp_empty = cpp_is_empty_kmer(~0ull);
    bool rust_empty = ragc_is_empty_kmer(~0ull);

    cout << "  C++ is_empty(~0ull): " << cpp_empty << "\n";
    cout << "  Rust is_empty(~0ull): " << rust_empty << "\n";

    bool cpp_not_empty = cpp_is_empty_kmer(100);
    bool rust_not_empty = ragc_is_empty_kmer(100);

    cout << "  C++ is_empty(100): " << cpp_not_empty << "\n";
    cout << "  Rust is_empty(100): " << rust_not_empty << "\n";

    if (cpp_empty == rust_empty && cpp_not_empty == rust_not_empty) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 5: Valid k-mer pair check
    cout << "Test 5: Valid k-mer pair check\n";
    bool cpp_valid = cpp_is_valid_kmer_pair(100, 200);
    bool rust_valid = ragc_is_valid_kmer_pair(100, 200);

    cout << "  C++ is_valid(100, 200): " << cpp_valid << "\n";
    cout << "  Rust is_valid(100, 200): " << rust_valid << "\n";

    bool cpp_invalid1 = cpp_is_valid_kmer_pair(~0ull, 200);
    bool rust_invalid1 = ragc_is_valid_kmer_pair(~0ull, 200);

    cout << "  C++ is_valid(~0ull, 200): " << cpp_invalid1 << "\n";
    cout << "  Rust is_valid(~0ull, 200): " << rust_invalid1 << "\n";

    bool cpp_invalid2 = cpp_is_valid_kmer_pair(100, ~0ull);
    bool rust_invalid2 = ragc_is_valid_kmer_pair(100, ~0ull);

    cout << "  C++ is_valid(100, ~0ull): " << cpp_invalid2 << "\n";
    cout << "  Rust is_valid(100, ~0ull): " << rust_invalid2 << "\n";

    if (cpp_valid == rust_valid && cpp_invalid1 == rust_invalid1 && cpp_invalid2 == rust_invalid2) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 6: Empty k-mer pair
    cout << "Test 6: Empty k-mer pair (~0ull, ~0ull)\n";
    KmerPair rust_empty_pair = ragc_create_empty_kmer_pair();

    cout << "  Rust empty pair: (" << rust_empty_pair.first << ", " << rust_empty_pair.second << ")\n";
    cout << "  Expected: (" << ~0ull << ", " << ~0ull << ")\n";

    if (rust_empty_pair.first == ~0ull && rust_empty_pair.second == ~0ull) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 7: Pair equality
    cout << "Test 7: K-mer pair equality\n";
    bool equal1 = ragc_kmer_pair_equals(100, 200, 100, 200);
    bool equal2 = ragc_kmer_pair_equals(100, 200, 100, 201);
    bool equal3 = ragc_kmer_pair_equals(100, 200, 200, 100);

    cout << "  (100,200) == (100,200): " << equal1 << " (expect true)\n";
    cout << "  (100,200) == (100,201): " << equal2 << " (expect false)\n";
    cout << "  (100,200) == (200,100): " << equal3 << " (expect false)\n";

    if (equal1 == true && equal2 == false && equal3 == false) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    // Test 8: Realistic k-mer values (21-mers)
    cout << "Test 8: Realistic 21-mer values\n";
    uint64_t kmer_a = 0x123456789ABULL;  // Example 21-mer
    uint64_t kmer_b = 0x987654321DCULL;

    cpp_pair = cpp_create_kmer_pair(kmer_a, kmer_b);
    rust_pair = ragc_create_kmer_pair(kmer_a, kmer_b);

    cout << "  C++ pair: (" << hex << cpp_pair.first << ", " << cpp_pair.second << dec << ")\n";
    cout << "  Rust pair: (" << hex << rust_pair.first << ", " << rust_pair.second << dec << ")\n";

    if (cpp_pair.first == rust_pair.first && cpp_pair.second == rust_pair.second) {
        cout << "  ✅ MATCH\n";
    } else {
        cout << "  ❌ MISMATCH\n";
        return 1;
    }
    cout << "\n";

    cout << "=== All tests passed! ===\n";
    cout << "✅ Rust k-mer pair operations match C++ AGC exactly\n";

    return 0;
}
