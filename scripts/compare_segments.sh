#!/bin/bash
# Compare actual segment creation between C++ AGC and RAGC

set -e

echo "=== Comparing Actual Segments Created ==="
echo ""

DATA_DIR="/home/erik/scrapy/yeast_split_proper"
CPP_AGC="/home/erik/agc/bin/agc"
RAGC="/home/erik/ragc/target/release/ragc"

# Create archives with verbose logging
echo "1. Creating C++ AGC archive with logging..."
cd "$DATA_DIR"
$CPP_AGC create -o /tmp/cpp_segments.agc -v 2 AAA#0.fa AAB#0.fa 2>&1 | tee /tmp/cpp_segments.log | tail -5

echo ""
echo "2. Creating RAGC archive with logging..."
cd "$DATA_DIR"
$RAGC create -o /tmp/ragc_segments.agc -v 2 AAA#0.fa AAB#0.fa 2>&1 | tee /tmp/ragc_segments.log | tail -5

echo ""
echo "3. Comparing segment statistics..."
echo ""
echo "C++ AGC:"
grep -E "No\. segments|splitters:" /tmp/cpp_segments.log

echo ""
echo "RAGC:"
grep -E "segments|splitter" /tmp/ragc_segments.log | grep -v "DEBUG:" | grep -v "Pass"

echo ""
echo "4. Archive sizes:"
ls -lh /tmp/cpp_segments.agc /tmp/ragc_segments.agc | awk '{print "  " $9 ": " $5}'

echo ""
echo "5. Decompression test:"
$CPP_AGC getset /tmp/cpp_segments.agc "AAA#0" 2>/dev/null > /tmp/cpp_aaa.fa
$RAGC getset /tmp/ragc_segments.agc "AAA#0" 2>/dev/null > /tmp/ragc_aaa.fa

CPP_SHA=$(grep -v "^>" /tmp/cpp_aaa.fa | tr -d '\n' | sha256sum | awk '{print $1}')
RAGC_SHA=$(grep -v "^>" /tmp/ragc_aaa.fa | tr -d '\n' | sha256sum | awk '{print $1}')
ORIG_SHA=$(grep -v "^>" "$DATA_DIR/AAA#0.fa" | tr -d '\n' | sha256sum | awk '{print $1}')

echo "  Original:  $ORIG_SHA"
echo "  C++ AGC:   $CPP_SHA"
echo "  RAGC:      $RAGC_SHA"

if [ "$CPP_SHA" = "$ORIG_SHA" ] && [ "$RAGC_SHA" = "$ORIG_SHA" ]; then
    echo "  ✓ Both decompress correctly"
else
    echo "  ✗ Decompression error!"
fi

echo ""
echo "CRITICAL QUESTION: Do they create the same number of segments?"
echo "(This would prove same split pattern)"
