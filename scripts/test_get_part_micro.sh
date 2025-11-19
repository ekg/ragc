#!/bin/bash
# Test the first micro-replacement: get_part() function
# Verifies that Rust ragc_get_contig_part() produces identical output to C++ get_part()

set -e

echo "=== Building Rust library ==="
cargo build --release

echo ""
echo "=== Compiling C++ test ==="
g++ -std=c++17 -O2 \
    -o /tmp/test_get_part \
    ragc-core/examples/test_get_part.cpp \
    -L./target/release \
    -lragc_core \
    -lpthread -ldl

echo ""
echo "=== Running test ==="
LD_LIBRARY_PATH=./target/release /tmp/test_get_part

echo ""
echo "Test binary: /tmp/test_get_part"
