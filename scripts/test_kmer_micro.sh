#!/bin/bash
# Test the second micro-replacement: k-mer extraction and canonicalization
# Verifies that Rust k-mer operations produce identical output to C++ CKmer

set -e

echo "=== Building Rust library ==="
cargo build --release

echo ""
echo "=== Compiling C++ test ==="
HOME_DIR=$(cd ~ && pwd)

g++ -std=c++17 -O2 \
    -o /tmp/test_kmer_extraction \
    -I./agc/src/common \
    -I./agc/src/core \
    -I$HOME_DIR/agc/3rd_party/refresh/compression/lib \
    ragc-core/examples/test_kmer_extraction.cpp \
    -L./target/release \
    -lragc_core \
    -lpthread -ldl

echo ""
echo "=== Running test ==="
LD_LIBRARY_PATH=./target/release /tmp/test_kmer_extraction

echo ""
echo "Test binary: /tmp/test_kmer_extraction"
