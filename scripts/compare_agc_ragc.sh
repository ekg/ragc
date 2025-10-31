#!/bin/bash
# Standard comparison test: RAGC vs C++ AGC
# Tests both performance and correctness

set -e

# Parse arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_files...>"
    echo "Example: $0 genome1.fa genome2.fa genome3.fa"
    echo "Example: $0 *.fa"
    exit 1
fi

INPUT_FILES="$@"
RAGC_ARCHIVE="/tmp/test_ragc.agc"
AGC_ARCHIVE="/tmp/test_agc.agc"

# Determine RAGC and AGC binaries
RAGC_BIN="${RAGC_BIN:-$(pwd)/target/release/ragc}"
AGC_BIN="${AGC_BIN:-/home/erik/agc/bin/agc}"

if [ ! -f "$RAGC_BIN" ]; then
    echo "Error: RAGC binary not found at $RAGC_BIN"
    echo "Set RAGC_BIN environment variable or build with: cargo build --release"
    exit 1
fi

if [ ! -f "$AGC_BIN" ]; then
    echo "Error: C++ AGC binary not found at $AGC_BIN"
    echo "Set AGC_BIN environment variable"
    exit 1
fi

# Standard parameters
KMER_LENGTH=21
SEGMENT_SIZE=10000
MIN_MATCH=20

echo "==============================================="
echo "  AGC vs RAGC Comparison Test"
echo "==============================================="
echo "Input files: $INPUT_FILES"
echo "RAGC binary: $RAGC_BIN"
echo "AGC binary:  $AGC_BIN"
echo "Parameters:  k=$KMER_LENGTH, s=$SEGMENT_SIZE, m=$MIN_MATCH"
echo ""

# Clean up old archives
rm -f "$RAGC_ARCHIVE" "$AGC_ARCHIVE"

# Test RAGC
echo "=== Creating RAGC archive ==="
/usr/bin/time -v "$RAGC_BIN" create \
    -o "$RAGC_ARCHIVE" \
    -k "$KMER_LENGTH" \
    -s "$SEGMENT_SIZE" \
    -m "$MIN_MATCH" \
    -v 0 \
    $INPUT_FILES \
    2>&1 | grep -E "(Elapsed|Maximum resident|Exit status)" | sed 's/^/RAGC: /'

RAGC_SIZE=$(stat -c %s "$RAGC_ARCHIVE")
echo "RAGC archive size: $RAGC_SIZE bytes"
echo ""

# Test C++ AGC
echo "=== Creating C++ AGC archive ==="
/usr/bin/time -v "$AGC_BIN" create \
    -k "$KMER_LENGTH" \
    -s "$SEGMENT_SIZE" \
    -l "$MIN_MATCH" \
    -o "$AGC_ARCHIVE" \
    $INPUT_FILES \
    2>&1 | grep -E "(Elapsed|Maximum resident|Exit status)" | sed 's/^/AGC:  /'

AGC_SIZE=$(stat -c %s "$AGC_ARCHIVE")
echo "AGC archive size: $AGC_SIZE bytes"
echo ""

# Compare sizes
echo "=== Size Comparison ==="
python3 <<EOF
ragc = $RAGC_SIZE
agc = $AGC_SIZE
diff = ragc - agc
pct = (diff / agc) * 100
print(f"RAGC:    {ragc:,} bytes")
print(f"C++ AGC: {agc:,} bytes")
print(f"Diff:    {diff:+,} bytes ({pct:+.2f}%)")
if abs(pct) < 2:
    print("✅ Size difference < 2%")
elif diff < 0:
    print("🎉 RAGC produces smaller archives!")
else:
    print("⚠️  RAGC produces larger archives")
EOF
echo ""

# Test correctness (first sample only)
echo "=== Testing Correctness (first sample) ==="
FIRST_FILE=$(echo $INPUT_FILES | awk '{print $1}')
SAMPLE_NAME=$("$RAGC_BIN" listset "$RAGC_ARCHIVE" 2>/dev/null | head -1)

if [ -z "$SAMPLE_NAME" ]; then
    echo "⚠️  Cannot determine sample name for correctness test"
else
    echo "Testing sample: $SAMPLE_NAME"

    # Test RAGC
    "$RAGC_BIN" getset -o /tmp/ragc_extract.fa "$RAGC_ARCHIVE" "$SAMPLE_NAME" 2>/dev/null
    RAGC_EXTRACT_SIZE=$(stat -c %s /tmp/ragc_extract.fa)

    # Test C++ AGC
    "$AGC_BIN" getset "$AGC_ARCHIVE" "$SAMPLE_NAME" > /tmp/agc_extract.fa 2>/dev/null
    AGC_EXTRACT_SIZE=$(stat -c %s /tmp/agc_extract.fa)

    ORIG_SIZE=$(stat -c %s "$FIRST_FILE")

    echo "Original:    $ORIG_SIZE bytes"
    echo "RAGC extract: $RAGC_EXTRACT_SIZE bytes"
    echo "AGC extract:  $AGC_EXTRACT_SIZE bytes"

    if [ "$RAGC_EXTRACT_SIZE" -eq "$ORIG_SIZE" ] && [ "$AGC_EXTRACT_SIZE" -eq "$ORIG_SIZE" ]; then
        if diff -q "$FIRST_FILE" /tmp/ragc_extract.fa >/dev/null 2>&1; then
            echo "✅ RAGC: Byte-for-byte identical"
        else
            echo "⚠️  RAGC: Size matches but content differs"
        fi

        if diff -q "$FIRST_FILE" /tmp/agc_extract.fa >/dev/null 2>&1; then
            echo "✅ C++ AGC: Byte-for-byte identical"
        else
            echo "⚠️  C++ AGC: Size matches but content differs"
        fi
    else
        if [ "$RAGC_EXTRACT_SIZE" -ne "$ORIG_SIZE" ]; then
            RAGC_DIFF=$((RAGC_EXTRACT_SIZE - ORIG_SIZE))
            echo "❌ RAGC: Size mismatch (${RAGC_DIFF:+} bytes)"
        fi
        if [ "$AGC_EXTRACT_SIZE" -ne "$ORIG_SIZE" ]; then
            AGC_DIFF=$((AGC_EXTRACT_SIZE - ORIG_SIZE))
            echo "❌ C++ AGC: Size mismatch (${AGC_DIFF:+} bytes)"
        fi
    fi
fi

echo ""
echo "==============================================="
echo "  Test Complete"
echo "==============================================="

# Cleanup
rm -f /tmp/ragc_extract.fa /tmp/agc_extract.fa
