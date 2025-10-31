#!/bin/bash
# Compare segment group creation between C++ AGC and RAGC

set -e

echo "=== Dumping Segment Groups for Comparison ==="
echo ""

DATA_DIR="/home/erik/scrapy/yeast_split_proper"
CPP_AGC="/home/erik/agc/bin/agc"
RAGC="/home/erik/ragc/target/release/ragc"

# Create archives with verbose output
echo "1. Creating C++ AGC archive..."
cd "$DATA_DIR"
$CPP_AGC create -o /tmp/cpp_dump.agc AAA#0.fa AAB#0.fa 2>&1 > /tmp/cpp_create.log

# Count segments mentioned
echo "   Analyzing C++ AGC output..."
grep -c "CPP_SPLIT" /tmp/cpp_create.log || echo "   Split attempts logged"

echo ""
echo "2. Creating RAGC archive with group logging..."
cd "$DATA_DIR"
$RAGC create -o /tmp/ragc_dump.agc -k 21 -s 10000 -m 20 -c 17 -t 1 -v 2 AAA#0.fa AAB#0.fa 2>&1 > /tmp/ragc_create.log

# Extract group counts
echo "   Analyzing RAGC output..."
grep "total unique groups" /tmp/ragc_create.log | tail -1

echo ""
echo "3. Comparing archive sizes..."
ls -lh /tmp/cpp_dump.agc /tmp/ragc_dump.agc | awk '{print "   " $9 ": " $5}'

echo ""
echo "4. Testing decompression..."
$CPP_AGC getset /tmp/cpp_dump.agc "AAB#0" > /tmp/cpp_aab.fa 2>/dev/null
$RAGC getset /tmp/ragc_dump.agc "AAB#0" > /tmp/ragc_aab.fa 2>/dev/null

CPP_SHA=$(grep -v "^>" /tmp/cpp_aab.fa | tr -d '\n' | sha256sum | awk '{print $1}')
RAGC_SHA=$(grep -v "^>" /tmp/ragc_aab.fa | tr -d '\n' | sha256sum | awk '{print $1}')
ORIG_SHA=$(grep -v "^>" "$DATA_DIR/AAB#0.fa" | tr -d '\n' | sha256sum | awk '{print $1}')

echo "   Original SHA256: $ORIG_SHA"
echo "   C++ AGC SHA256:  $CPP_SHA"
echo "   RAGC SHA256:     $RAGC_SHA"

if [ "$CPP_SHA" = "$ORIG_SHA" ] && [ "$RAGC_SHA" = "$ORIG_SHA" ]; then
    echo "   ✓ Both decompress correctly"
else
    echo "   ✗ Decompression mismatch!"
fi

echo ""
echo "Logs saved to:"
echo "  - /tmp/cpp_create.log"
echo "  - /tmp/ragc_create.log"
echo ""
echo "Next: Analyze the logs to find grouping differences"
