#!/bin/bash
# Test the seventh micro-replacement: segment splitting
# Verifies that Rust segment split arithmetic matches C++ AGC exactly

set -e

echo "=== Building Rust library ==="
cargo build --release

echo ""
echo "=== Compiling C++ test ==="
g++ -std=c++17 -O2 \
    -o /tmp/test_segment_split \
    ragc-core/examples/test_segment_split.cpp \
    -L./target/release \
    -lragc_core \
    -lpthread -ldl

echo ""
echo "=== Running test ==="
LD_LIBRARY_PATH=./target/release /tmp/test_segment_split

echo ""
echo "Test binary: /tmp/test_segment_split"
