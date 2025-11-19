#!/bin/bash
# Test the sixth micro-replacement: reverse complement
# Verifies that Rust reverse complement matches C++ AGC exactly

set -e

echo "=== Building Rust library ==="
cargo build --release

echo ""
echo "=== Compiling C++ test ==="
g++ -std=c++17 -O2 \
    -o /tmp/test_reverse_complement \
    ragc-core/examples/test_reverse_complement.cpp \
    -L./target/release \
    -lragc_core \
    -lpthread -ldl

echo ""
echo "=== Running test ==="
LD_LIBRARY_PATH=./target/release /tmp/test_reverse_complement

echo ""
echo "Test binary: /tmp/test_reverse_complement"
