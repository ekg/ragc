#!/bin/bash
# Test the third micro-replacement: splitter checking
# Verifies that Rust splitter lookup produces identical results to C++ AGC

set -e

echo "=== Building Rust library ==="
cargo build --release

echo ""
echo "=== Compiling C++ test ==="
g++ -std=c++17 -O2 \
    -o /tmp/test_splitter_check \
    ragc-core/examples/test_splitter_check.cpp \
    -L./target/release \
    -lragc_core \
    -lpthread -ldl

echo ""
echo "=== Running test ==="
LD_LIBRARY_PATH=./target/release /tmp/test_splitter_check

echo ""
echo "Test binary: /tmp/test_splitter_check"
