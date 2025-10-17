#!/bin/bash
# Cross-compatibility test: Verify Rust AGC ↔ C++ AGC compatibility

set -e

RUST_AGC="./target/release/ragc"
CPP_AGC="/home/erik/agc/bin/agc"
TEST_DATA="/home/erik/wfmash.3/data/scerevisiae8.fa.gz"
TEST_DIR="/tmp/agc_compat_test"

echo "=== AGC Cross-Compatibility Test ==="
echo ""

# Clean up and create test directory
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"

# Extract original data for comparison
echo "1. Extracting original data for comparison..."
zcat "$TEST_DATA" | grep -v "^>" | tr -d '\n' > "$TEST_DIR/original.txt"
ORIGINAL_SHA=$(sha256sum "$TEST_DIR/original.txt" | cut -d' ' -f1)
echo "   Original SHA-256: $ORIGINAL_SHA"
echo ""

# Test 1: Rust creates, C++ extracts
echo "2. Test: Rust AGC creates → C++ AGC extracts"
$RUST_AGC create -o "$TEST_DIR/rust_created.agc" "$TEST_DATA" 2>&1 | grep -E "(Created|complete|error)" || true
RUST_SIZE=$(stat -c%s "$TEST_DIR/rust_created.agc")
echo "   Rust archive size: $(($RUST_SIZE / 1024 / 1024))MB"

$CPP_AGC getset "$TEST_DIR/rust_created.agc" scerevisiae8 | grep -v "^>" | tr -d '\n' > "$TEST_DIR/cpp_extracted.txt"
CPP_EXTRACTED_SHA=$(sha256sum "$TEST_DIR/cpp_extracted.txt" | cut -d' ' -f1)

if [ "$ORIGINAL_SHA" = "$CPP_EXTRACTED_SHA" ]; then
    echo "   ✓ C++ AGC successfully extracted Rust archive"
else
    echo "   ✗ FAIL: SHA-256 mismatch!"
    echo "     Expected: $ORIGINAL_SHA"
    echo "     Got:      $CPP_EXTRACTED_SHA"
    exit 1
fi
echo ""

# Test 2: C++ creates, Rust extracts
echo "3. Test: C++ AGC creates → Rust AGC extracts"
$CPP_AGC create "$TEST_DIR/cpp_created.agc" "$TEST_DATA" 2>&1 | grep -E "(complete|error)" || true
CPP_SIZE=$(stat -c%s "$TEST_DIR/cpp_created.agc")
echo "   C++ archive size: $(($CPP_SIZE / 1024 / 1024))MB"

$RUST_AGC getset "$TEST_DIR/cpp_created.agc" scerevisiae8.fa 2>/dev/null | grep -v "^>" | tr -d '\n' > "$TEST_DIR/rust_extracted.txt"
RUST_EXTRACTED_SHA=$(sha256sum "$TEST_DIR/rust_extracted.txt" | cut -d' ' -f1)

if [ "$ORIGINAL_SHA" = "$RUST_EXTRACTED_SHA" ]; then
    echo "   ✓ Rust AGC successfully extracted C++ archive"
else
    echo "   ✗ FAIL: SHA-256 mismatch!"
    echo "     Expected: $ORIGINAL_SHA"
    echo "     Got:      $RUST_EXTRACTED_SHA"
    exit 1
fi
echo ""

# Summary
echo "=== Summary ==="
echo "✓ All cross-compatibility tests passed!"
echo ""
echo "Archive sizes:"
echo "  Rust AGC:  $(($RUST_SIZE / 1024 / 1024))MB"
echo "  C++ AGC:   $(($CPP_SIZE / 1024 / 1024))MB"
echo "  Ratio:     $(echo "scale=2; $RUST_SIZE / $CPP_SIZE" | bc)x"
echo ""
echo "Both implementations can read each other's archives correctly."

# Cleanup
rm -rf "$TEST_DIR"
