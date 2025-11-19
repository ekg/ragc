#!/bin/bash
# Test the eighth micro-replacement: k-mer pair operations
# Verifies that Rust k-mer pair logic matches C++ AGC exactly

set -e

echo "=== Building Rust library ==="
cargo build --release

echo ""
echo "=== Compiling C++ test ==="
g++ -std=c++17 -O2 \
    -o /tmp/test_kmer_pair \
    ragc-core/examples/test_kmer_pair.cpp \
    -L./target/release \
    -lragc_core \
    -lpthread -ldl

echo ""
echo "=== Running test ==="
LD_LIBRARY_PATH=./target/release /tmp/test_kmer_pair

echo ""
echo "Test binary: /tmp/test_kmer_pair"
