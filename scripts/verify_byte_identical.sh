#!/bin/bash
set -euo pipefail

# verify_byte_identical.sh - Compare two AGC archives byte-for-byte
#
# Usage: ./verify_byte_identical.sh <archive1.agc> <archive2.agc>
#
# Checks:
# 1. SHA256 checksums (byte-identical test)
# 2. If different, shows byte offset where they diverge
# 3. Compares segment layouts using ragc inspect
# 4. Compares sample lists

if [ $# -ne 2 ]; then
    echo "Usage: $0 <archive1.agc> <archive2.agc>"
    echo ""
    echo "Compares two AGC archives for byte-identical equality."
    echo "If archives differ, shows WHERE they differ and extracts debug info."
    exit 1
fi

ARCHIVE1="$1"
ARCHIVE2="$2"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "=========================================="
echo "Archive Comparison Tool"
echo "=========================================="
echo ""
echo "Archive 1: $ARCHIVE1"
echo "Archive 2: $ARCHIVE2"
echo ""

# Check files exist
if [ ! -f "$ARCHIVE1" ]; then
    echo -e "${RED}ERROR: Archive 1 not found: $ARCHIVE1${NC}"
    exit 1
fi

if [ ! -f "$ARCHIVE2" ]; then
    echo -e "${RED}ERROR: Archive 2 not found: $ARCHIVE2${NC}"
    exit 1
fi

# Get file sizes
SIZE1=$(stat -c%s "$ARCHIVE1")
SIZE2=$(stat -c%s "$ARCHIVE2")
echo "Archive 1 size: $(numfmt --to=iec-i --suffix=B $SIZE1) ($SIZE1 bytes)"
echo "Archive 2 size: $(numfmt --to=iec-i --suffix=B $SIZE2) ($SIZE2 bytes)"

if [ $SIZE1 -ne $SIZE2 ]; then
    DIFF=$((SIZE1 - SIZE2))
    PERCENT=$(awk "BEGIN {printf \"%.2f\", ($DIFF / $SIZE2) * 100}")
    echo -e "${YELLOW}Size difference: $DIFF bytes ($PERCENT%)${NC}"
else
    echo -e "${GREEN}Sizes match!${NC}"
fi
echo ""

# Compare SHA256 checksums
echo "Computing SHA256 checksums..."
SHA1=$(sha256sum "$ARCHIVE1" | awk '{print $1}')
SHA2=$(sha256sum "$ARCHIVE2" | awk '{print $1}')

echo "Archive 1: $SHA1"
echo "Archive 2: $SHA2"
echo ""

if [ "$SHA1" = "$SHA2" ]; then
    echo -e "${GREEN}✅ SUCCESS: Archives are BYTE-IDENTICAL!${NC}"
    echo ""
    echo "SHA256 checksums match exactly. Archives are 100% identical."
    exit 0
fi

echo -e "${RED}❌ FAILED: Archives differ!${NC}"
echo ""

# Find byte offset where they differ
echo "Finding first byte difference..."
DIFF_OUTPUT=$(cmp -l "$ARCHIVE1" "$ARCHIVE2" 2>&1 | head -1 || true)

if [ -n "$DIFF_OUTPUT" ]; then
    OFFSET=$(echo "$DIFF_OUTPUT" | awk '{print $1}')
    BYTE1=$(echo "$DIFF_OUTPUT" | awk '{print $2}')
    BYTE2=$(echo "$DIFF_OUTPUT" | awk '{print $3}')
    echo -e "${YELLOW}First difference at byte offset: $OFFSET${NC}"
    echo "  Archive 1: byte value $BYTE1 (octal)"
    echo "  Archive 2: byte value $BYTE2 (octal)"
    echo ""
fi

# Create temp directory for comparison outputs
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "Extracting archive metadata for comparison..."
echo ""

# Compare sample lists
echo "=========================================="
echo "Sample List Comparison"
echo "=========================================="

if command -v ragc &> /dev/null; then
    ragc listset "$ARCHIVE1" > "$TMPDIR/samples1.txt" 2>&1 || echo "Failed to list samples from archive 1"
    ragc listset "$ARCHIVE2" > "$TMPDIR/samples2.txt" 2>&1 || echo "Failed to list samples from archive 2"

    if [ -f "$TMPDIR/samples1.txt" ] && [ -f "$TMPDIR/samples2.txt" ]; then
        echo "Archive 1 samples:"
        cat "$TMPDIR/samples1.txt" | head -10
        echo ""
        echo "Archive 2 samples:"
        cat "$TMPDIR/samples2.txt" | head -10
        echo ""

        if diff -q "$TMPDIR/samples1.txt" "$TMPDIR/samples2.txt" > /dev/null 2>&1; then
            echo -e "${GREEN}✅ Sample lists match${NC}"
        else
            echo -e "${RED}❌ Sample lists differ:${NC}"
            diff -u "$TMPDIR/samples1.txt" "$TMPDIR/samples2.txt" || true
        fi
        echo ""
    fi
fi

# Compare segment layouts if ragc inspect supports it
echo "=========================================="
echo "Segment Layout Comparison"
echo "=========================================="

if command -v ragc &> /dev/null; then
    # Try to export segment layouts
    if ragc inspect "$ARCHIVE1" --segment-layout > "$TMPDIR/layout1.csv" 2>&1; then
        if ragc inspect "$ARCHIVE2" --segment-layout > "$TMPDIR/layout2.csv" 2>&1; then
            echo "Comparing segment layouts..."

            LAYOUT1_LINES=$(wc -l < "$TMPDIR/layout1.csv")
            LAYOUT2_LINES=$(wc -l < "$TMPDIR/layout2.csv")

            echo "Archive 1: $LAYOUT1_LINES segments"
            echo "Archive 2: $LAYOUT2_LINES segments"
            echo ""

            if diff -q "$TMPDIR/layout1.csv" "$TMPDIR/layout2.csv" > /dev/null 2>&1; then
                echo -e "${GREEN}✅ Segment layouts match exactly${NC}"
            else
                echo -e "${RED}❌ Segment layouts differ:${NC}"
                echo ""
                echo "First 10 differences:"
                diff -u "$TMPDIR/layout1.csv" "$TMPDIR/layout2.csv" | head -30 || true
                echo ""
                echo -e "${BLUE}Full diff saved to: $TMPDIR/layout_diff.txt${NC}"
                diff -u "$TMPDIR/layout1.csv" "$TMPDIR/layout2.csv" > "$TMPDIR/layout_diff.txt" || true
            fi
        else
            echo "Note: ragc inspect --segment-layout not available for archive 2"
        fi
    else
        echo "Note: ragc inspect --segment-layout not available for archive 1"
    fi
fi

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo -e "${RED}Archives are NOT byte-identical${NC}"
echo ""
echo "Debug artifacts saved to: $TMPDIR"
echo "  - samples1.txt, samples2.txt: Sample lists"
echo "  - layout1.csv, layout2.csv: Segment layouts (if available)"
echo "  - layout_diff.txt: Detailed segment differences"
echo ""
echo "To investigate further:"
echo "  1. Check segment layouts for first divergence point"
echo "  2. Compare splitter detection logs"
echo "  3. Compare compression parameters"
echo "  4. Verify both archives can be read correctly"

exit 1
