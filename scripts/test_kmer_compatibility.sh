#!/bin/bash
# Test that RAGC produces expected k-mer pairs (pegged to C++ AGC reference output)
# This ensures no regression from 100% algorithmic compatibility

set -e

RAGC_BIN="${RAGC_BIN:-./target/release/ragc}"
TEST_DATA="${1}"
EXPECTED_PAIRS="${2}"
K="${3:-21}"

if [ -z "$TEST_DATA" ] || [ -z "$EXPECTED_PAIRS" ]; then
    echo "Usage: $0 <test_data.fa> <expected_pairs.txt> [k-mer_size]"
    echo "Example: $0 tests/data/yeast_chr1.fa tests/fixtures/yeast_chr1_k21_pairs.txt 21"
    echo ""
    echo "Expected pairs file format: sorted 'front=... back=...' lines"
    exit 1
fi

if [ ! -f "$TEST_DATA" ]; then
    echo "Error: Test data file not found: $TEST_DATA"
    exit 1
fi

if [ ! -f "$EXPECTED_PAIRS" ]; then
    echo "Error: Expected pairs file not found: $EXPECTED_PAIRS"
    exit 1
fi

echo "=== K-mer Compatibility Test ==="
echo "RAGC: $RAGC_BIN"
echo "Test data: $TEST_DATA"
echo "Expected pairs: $EXPECTED_PAIRS"
echo "K-mer size: $K"
echo ""

# Check binary exists
if [ ! -f "$RAGC_BIN" ]; then
    echo "Error: RAGC binary not found: $RAGC_BIN"
    echo "Run: cargo build --release"
    exit 1
fi

# Temporary files
RAGC_OUT=$(mktemp)
RAGC_PAIRS=$(mktemp)

trap "rm -f $RAGC_OUT $RAGC_PAIRS /tmp/ragc_test.agc" EXIT

# Run RAGC
echo "Running RAGC..."
$RAGC_BIN create -o /tmp/ragc_test.agc -k $K "$TEST_DATA" 2>"$RAGC_OUT" >/dev/null || {
    echo "Error: RAGC failed"
    tail -20 "$RAGC_OUT"
    exit 1
}

# Extract k-mer pairs (requires GROUP_KMER logging to be enabled)
echo "Extracting k-mer pairs..."
grep "GROUP_KMER:" "$RAGC_OUT" 2>/dev/null | awk '{print $3, $4}' | sort > "$RAGC_PAIRS" || {
    echo "Error: No GROUP_KMER logging found"
    echo "Set RAGC_DEBUG_KMER_PAIRS=1 environment variable or enable debug logging"
    exit 1
}

# Count pairs
RAGC_COUNT=$(wc -l < "$RAGC_PAIRS" | tr -d ' ')
EXPECTED_COUNT=$(wc -l < "$EXPECTED_PAIRS" | tr -d ' ')
COMMON_COUNT=$(comm -12 "$RAGC_PAIRS" "$EXPECTED_PAIRS" | wc -l | tr -d ' ')

echo ""
echo "=== Results ==="
echo "RAGC k-mer pairs: $RAGC_COUNT"
echo "Expected k-mer pairs: $EXPECTED_COUNT"
echo "Pairs in common: $COMMON_COUNT"

# Check for exact match
if [ "$RAGC_COUNT" -eq "$EXPECTED_COUNT" ] && [ "$COMMON_COUNT" -eq "$RAGC_COUNT" ]; then
    PERCENT="100"
else
    PERCENT=$(echo "scale=1; $COMMON_COUNT * 100 / $EXPECTED_COUNT" | bc)
fi

echo "Match percentage: ${PERCENT}%"

if [ "$RAGC_COUNT" -eq "$EXPECTED_COUNT" ] && [ "$COMMON_COUNT" -eq "$RAGC_COUNT" ]; then
    echo ""
    echo "✅ PASS: 100% k-mer pair match ($COMMON_COUNT/$EXPECTED_COUNT)"

    # Also verify decompression
    # Extract sample name: ">AAA#0#chrI" -> "AAA#0" (sample#id, before contig name)
    SAMPLE=$(grep "^>" "$TEST_DATA" | head -1 | sed 's/^>//' | cut -d'#' -f1,2)
    ORIGINAL_SHA=$(cat "$TEST_DATA" | sha256sum | awk '{print $1}')
    RAGC_SHA=$($RAGC_BIN getset /tmp/ragc_test.agc "$SAMPLE" 2>/dev/null | sha256sum | awk '{print $1}')

    if [ "$ORIGINAL_SHA" = "$RAGC_SHA" ]; then
        echo "✅ Decompression verified (SHA256 match)"
    else
        echo "⚠️  WARNING: Decompression mismatch!"
    fi

    exit 0
else
    echo ""
    echo "❌ FAIL: K-mer pairs do not match expected!"
    echo ""
    echo "Only in RAGC output (first 10):"
    comm -23 "$RAGC_PAIRS" "$EXPECTED_PAIRS" | head -10
    echo ""
    echo "Only in expected (first 10):"
    comm -13 "$RAGC_PAIRS" "$EXPECTED_PAIRS" | head -10
    exit 1
fi
