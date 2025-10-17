#!/bin/bash
set -e

RUST_AGC="./target/release/ragc"
CPP_AGC="/home/erik/agc/bin/agc"
TEST_DATA="/home/erik/wfmash.3/data/scerevisiae8.fa.gz"
TEST_DIR="/tmp/agc_compat"

echo "=== AGC Cross-Compatibility Test ==="
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"

# Get original hash
echo "1. Extracting original data..."
zcat "$TEST_DATA" | grep -v "^>" | tr -d '\n' > "$TEST_DIR/original.txt"
ORIG_SHA=$(sha256sum "$TEST_DIR/original.txt" | awk '{print $1}')
echo "   Original SHA-256: $ORIG_SHA"

# Test 1: Rust → C++ 
echo ""
echo "2. Rust creates → C++ extracts"
$RUST_AGC create -o "$TEST_DIR/rust.agc" "$TEST_DATA" >/dev/null 2>&1
RUST_SIZE=$(stat -c%s "$TEST_DIR/rust.agc")
echo "   Rust archive: $((RUST_SIZE/1024/1024))MB"

$CPP_AGC getset "$TEST_DIR/rust.agc" scerevisiae8 2>/dev/null | grep -v "^>" | tr -d '\n' > "$TEST_DIR/cpp_out.txt"
CPP_SHA=$(sha256sum "$TEST_DIR/cpp_out.txt" | awk '{print $1}')

if [ "$ORIG_SHA" = "$CPP_SHA" ]; then
    echo "   ✓ C++ AGC extracted Rust archive correctly"
else
    echo "   ✗ FAIL: Mismatch!"
    exit 1
fi

# Test 2: C++ → Rust
echo ""
echo "3. C++ creates → Rust extracts"
$CPP_AGC create "$TEST_DIR/cpp.agc" "$TEST_DATA" >/dev/null 2>&1
CPP_SIZE=$(stat -c%s "$TEST_DIR/cpp.agc")
echo "   C++ archive: $((CPP_SIZE/1024/1024))MB"

$RUST_AGC getset "$TEST_DIR/cpp.agc" scerevisiae8.fa 2>/dev/null | grep -v "^>" | tr -d '\n' > "$TEST_DIR/rust_out.txt"
RUST_SHA=$(sha256sum "$TEST_DIR/rust_out.txt" | awk '{print $1}')

if [ "$ORIG_SHA" = "$RUST_SHA" ]; then
    echo "   ✓ Rust AGC extracted C++ archive correctly"
else
    echo "   ✗ FAIL: Mismatch!"
    exit 1
fi

echo ""
echo "=== ✓ All tests passed! ==="
echo "Compression: Rust $((RUST_SIZE/1024/1024))MB vs C++ $((CPP_SIZE/1024/1024))MB"

rm -rf "$TEST_DIR"
