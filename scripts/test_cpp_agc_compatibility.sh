#!/bin/bash
set -euo pipefail

# test_cpp_agc_compatibility.sh - Test RAGC vs C++ AGC compatibility
#
# Tests:
# 1. Creates archives with both implementations
# 2. Verifies byte-identical output
# 3. Tests cross-compatibility (RAGC reads C++, C++ reads RAGC)
# 4. Verifies extraction produces identical sequences

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Paths
RAGC="./target/release/ragc"
CPP_AGC="/home/erik/agc/bin/agc"
VERIFY_SCRIPT="./scripts/verify_byte_identical.sh"

# Default test data
DEFAULT_TEST_DATA="/home/erik/scrapy/cpp_AAA.fa"
TEST_DATA="${1:-$DEFAULT_TEST_DATA}"

# Check if test data exists, if not create minimal test
if [ ! -f "$TEST_DATA" ]; then
    echo -e "${YELLOW}Test data not found: $TEST_DATA${NC}"
    echo "Creating minimal test data..."
    TEST_DATA="/tmp/minimal_test.fa"

    # Generate 10KB test sequence (random DNA)
    cat > "$TEST_DATA" << 'EOF'
>test_sample#1#chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG
EOF

    # Repeat to reach ~10KB
    for i in {1..40}; do
        tail -n +2 "$TEST_DATA" >> "$TEST_DATA.tmp"
    done
    mv "$TEST_DATA.tmp" "$TEST_DATA"

    echo -e "${GREEN}Created test data: $TEST_DATA ($(stat -c%s "$TEST_DATA") bytes)${NC}"
    echo ""
fi

echo "=========================================="
echo "C++ AGC Compatibility Test"
echo "=========================================="
echo ""
echo "Test data: $TEST_DATA"
echo "Test size: $(stat -c%s "$TEST_DATA" | numfmt --to=iec-i --suffix=B)"
echo ""

# Check binaries exist
if [ ! -f "$RAGC" ]; then
    echo -e "${RED}ERROR: RAGC binary not found: $RAGC${NC}"
    echo "Build it with: cargo build --release"
    exit 1
fi

if [ ! -f "$CPP_AGC" ]; then
    echo -e "${RED}ERROR: C++ AGC binary not found: $CPP_AGC${NC}"
    echo "Expected at: $CPP_AGC"
    exit 1
fi

# Create temp directory
TMPDIR=$(mktemp -d)
# Only clean up on success - preserve on failure for debugging
cleanup() {
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        rm -rf "$TMPDIR"
    else
        echo ""
        echo "Temp directory preserved for debugging: $TMPDIR"
    fi
}
trap cleanup EXIT

RAGC_ARCHIVE="$TMPDIR/ragc_test.agc"
CPP_ARCHIVE="$TMPDIR/cpp_test.agc"

# Standard parameters (matching between implementations)
# RAGC: -k (kmer) -s (step) -m (min-seg)
# C++ AGC: -k (kmer) -s (step) -l (min-len, equivalent to min-seg)
KMER=21
STEP=10000
MINSEG=20
THREADS=1  # Single-threaded for determinism

echo "=========================================="
echo "Step 1: Create archive with RAGC"
echo "=========================================="
echo "Command: $RAGC create -o $RAGC_ARCHIVE -k $KMER -s $STEP -m $MINSEG -t $THREADS $TEST_DATA"
echo ""

/usr/bin/time -v "$RAGC" create -o "$RAGC_ARCHIVE" -k "$KMER" -s "$STEP" -m "$MINSEG" -t "$THREADS" "$TEST_DATA" 2>&1 | grep -E "(Maximum resident|User time|System time)" || true
echo ""

if [ ! -f "$RAGC_ARCHIVE" ]; then
    echo -e "${RED}ERROR: RAGC failed to create archive${NC}"
    exit 1
fi

RAGC_SIZE=$(stat -c%s "$RAGC_ARCHIVE")
echo -e "${GREEN}✅ RAGC archive created: $(numfmt --to=iec-i --suffix=B $RAGC_SIZE)${NC}"
echo ""

echo "=========================================="
echo "Step 2: Create archive with C++ AGC"
echo "=========================================="
echo "Command: $CPP_AGC create -o $CPP_ARCHIVE -k $KMER -s $STEP -l $MINSEG -t $THREADS $TEST_DATA"
echo ""

/usr/bin/time -v "$CPP_AGC" create -o "$CPP_ARCHIVE" -k "$KMER" -s "$STEP" -l "$MINSEG" -t "$THREADS" "$TEST_DATA" 2>&1 | grep -E "(Maximum resident|User time|System time)" || true
echo ""

if [ ! -f "$CPP_ARCHIVE" ]; then
    echo -e "${RED}ERROR: C++ AGC failed to create archive${NC}"
    exit 1
fi

CPP_SIZE=$(stat -c%s "$CPP_ARCHIVE")
echo -e "${GREEN}✅ C++ AGC archive created: $(numfmt --to=iec-i --suffix=B $CPP_SIZE)${NC}"
echo ""

# Calculate size difference
SIZE_DIFF=$((RAGC_SIZE - CPP_SIZE))
SIZE_PERCENT=$(awk "BEGIN {printf \"%.2f\", ($SIZE_DIFF / $CPP_SIZE) * 100}")

echo "Archive size comparison:"
echo "  RAGC:   $(numfmt --to=iec-i --suffix=B $RAGC_SIZE) ($RAGC_SIZE bytes)"
echo "  C++ AGC: $(numfmt --to=iec-i --suffix=B $CPP_SIZE) ($CPP_SIZE bytes)"
echo "  Difference: $SIZE_DIFF bytes ($SIZE_PERCENT%)"
echo ""

echo "=========================================="
echo "Step 3: Compare archives byte-by-byte"
echo "=========================================="

if [ -f "$VERIFY_SCRIPT" ]; then
    "$VERIFY_SCRIPT" "$RAGC_ARCHIVE" "$CPP_ARCHIVE"
    BYTE_IDENTICAL=$?
else
    echo "Verification script not found, doing basic comparison..."
    RAGC_SHA=$(sha256sum "$RAGC_ARCHIVE" | awk '{print $1}')
    CPP_SHA=$(sha256sum "$CPP_ARCHIVE" | awk '{print $1}')

    echo "RAGC SHA256:   $RAGC_SHA"
    echo "C++ AGC SHA256: $CPP_SHA"
    echo ""

    if [ "$RAGC_SHA" = "$CPP_SHA" ]; then
        echo -e "${GREEN}✅ Archives are BYTE-IDENTICAL!${NC}"
        BYTE_IDENTICAL=0
    else
        echo -e "${YELLOW}⚠️  Archives differ (not byte-identical)${NC}"
        BYTE_IDENTICAL=1
    fi
fi
echo ""

echo "=========================================="
echo "Step 4: Test cross-compatibility"
echo "=========================================="

# Get sample names from each archive (they may differ!)
RAGC_SAMPLE=$("$RAGC" listset "$RAGC_ARCHIVE" | head -1)
CPP_SAMPLE=$("$CPP_AGC" listset "$CPP_ARCHIVE" | head -1)

echo "Sample name in RAGC archive: $RAGC_SAMPLE"
echo "Sample name in C++ AGC archive: $CPP_SAMPLE"
echo ""

# Test 1: RAGC reads C++ archive
echo "Test 4a: RAGC reads C++ AGC archive..."
"$RAGC" listset "$CPP_ARCHIVE" > "$TMPDIR/ragc_reads_cpp_samples.txt" 2>&1 || {
    echo -e "${RED}ERROR: RAGC cannot list samples from C++ archive${NC}"
    exit 1
}
echo -e "${GREEN}✅ RAGC can read C++ AGC archive${NC}"

"$RAGC" getset "$CPP_ARCHIVE" "$CPP_SAMPLE" > "$TMPDIR/ragc_from_cpp.fa" 2>&1 || {
    echo -e "${RED}ERROR: RAGC cannot extract from C++ archive${NC}"
    exit 1
}
echo -e "${GREEN}✅ RAGC can extract sequences from C++ AGC archive${NC}"
echo ""

# Test 2: C++ AGC reads RAGC archive
echo "Test 4b: C++ AGC reads RAGC archive..."
"$CPP_AGC" listset "$RAGC_ARCHIVE" > "$TMPDIR/cpp_reads_ragc_samples.txt" 2>&1 || {
    echo -e "${RED}ERROR: C++ AGC cannot list samples from RAGC archive${NC}"
    exit 1
}
echo -e "${GREEN}✅ C++ AGC can read RAGC archive${NC}"

"$CPP_AGC" getset "$RAGC_ARCHIVE" "$RAGC_SAMPLE" > "$TMPDIR/cpp_from_ragc.fa" 2>&1 || {
    echo -e "${RED}ERROR: C++ AGC cannot extract from RAGC archive${NC}"
    exit 1
}
echo -e "${GREEN}✅ C++ AGC can extract sequences from RAGC archive${NC}"
echo ""

echo "=========================================="
echo "Step 5: Verify extracted sequences"
echo "=========================================="

# Compare extracted sequences
echo "Comparing extracted sequences..."

# Extract from RAGC archive using RAGC
"$RAGC" getset "$RAGC_ARCHIVE" "$RAGC_SAMPLE" > "$TMPDIR/ragc_from_ragc.fa" 2>&1

# Extract from C++ archive using C++ AGC
"$CPP_AGC" getset "$CPP_ARCHIVE" "$CPP_SAMPLE" > "$TMPDIR/cpp_from_cpp.fa" 2>&1

# Normalize FASTA (remove line breaks in sequences)
normalize_fasta() {
    awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {if (seq) print seq}' "$1"
}

normalize_fasta "$TMPDIR/ragc_from_ragc.fa" > "$TMPDIR/ragc_from_ragc_norm.fa"
normalize_fasta "$TMPDIR/ragc_from_cpp.fa" > "$TMPDIR/ragc_from_cpp_norm.fa"
normalize_fasta "$TMPDIR/cpp_from_ragc.fa" > "$TMPDIR/cpp_from_ragc_norm.fa"
normalize_fasta "$TMPDIR/cpp_from_cpp.fa" > "$TMPDIR/cpp_from_cpp_norm.fa"

# Compare all extractions
ALL_MATCH=true

echo "Checking: RAGC extracts from RAGC archive vs C++ extracts from C++ archive..."
if diff -q "$TMPDIR/ragc_from_ragc_norm.fa" "$TMPDIR/cpp_from_cpp_norm.fa" > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Sequences match${NC}"
else
    echo -e "${RED}❌ Sequences differ!${NC}"
    ALL_MATCH=false
fi

echo "Checking: RAGC extracts from C++ archive vs original..."
if diff -q "$TMPDIR/ragc_from_cpp_norm.fa" "$TMPDIR/cpp_from_cpp_norm.fa" > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Cross-extraction matches${NC}"
else
    echo -e "${RED}❌ Cross-extraction differs!${NC}"
    ALL_MATCH=false
fi

echo "Checking: C++ extracts from RAGC archive vs original..."
if diff -q "$TMPDIR/cpp_from_ragc_norm.fa" "$TMPDIR/ragc_from_ragc_norm.fa" > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Cross-extraction matches${NC}"
else
    echo -e "${RED}❌ Cross-extraction differs!${NC}"
    ALL_MATCH=false
fi
echo ""

echo "=========================================="
echo "Final Summary"
echo "=========================================="
echo ""

if [ $BYTE_IDENTICAL -eq 0 ]; then
    echo -e "${GREEN}🎉 PERFECT: Archives are BYTE-IDENTICAL!${NC}"
    echo ""
    echo "RAGC and C++ AGC produce exactly the same output."
    exit 0
elif $ALL_MATCH; then
    echo -e "${YELLOW}⚠️  COMPATIBLE but not byte-identical:${NC}"
    echo ""
    echo "✅ Both implementations can read each other's archives"
    echo "✅ Extracted sequences are identical"
    echo "✅ Format compatibility verified"
    echo "⚠️  Archive files differ by $SIZE_DIFF bytes ($SIZE_PERCENT%)"
    echo ""
    echo "This indicates algorithmic differences but functional compatibility."
    echo "For research/production use, aim for byte-identical archives."
    exit 1
else
    echo -e "${RED}❌ FAILED: Compatibility issues detected${NC}"
    echo ""
    echo "Archives have differences that affect extracted sequences."
    echo "Check debug output above for details."
    exit 2
fi
