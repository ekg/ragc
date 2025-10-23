#!/bin/bash
#
# Quick test: Channel buffer reduction from 100 → 10
# Compare against baseline (from /tmp/allocator_profile_20251022_215517/default_time.txt)

set -e

TEST_DIR="/home/erik/scrapy/yeast10_test"
OUTPUT_DIR="/tmp/buffer_test_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$OUTPUT_DIR"

echo "=== Testing Channel Buffer Reduction (100 \u2192 10) ==="
echo ""

# Prepare sample files
SAMPLES_DIR="$OUTPUT_DIR/samples"
mkdir -p "$SAMPLES_DIR"
for f in "$TEST_DIR"/*.fa; do
    basename_f=$(basename "$f")
    if [ "$basename_f" != "yeast10_pansn.fa" ]; then
        cp "$f" "$SAMPLES_DIR/"
    fi
done

# Test current build (with buffer=10)
echo "Running compression with reduced buffers (10)..."
ARCHIVE="$OUTPUT_DIR/archive.agc"
/usr/bin/time -v /home/erik/ragc/target/release/ragc create \
    -o "$ARCHIVE" \
    -k 21 \
    -s 10000 \
    -m 20 \
    -v 0 \
    "$SAMPLES_DIR"/*.fa \
    2>&1 | tee "$OUTPUT_DIR/time.txt"

echo ""
echo "=== RESULTS ==="
echo ""

MEM=$(grep "Maximum resident set size" "$OUTPUT_DIR/time.txt" | awk '{print $6}')
TIME=$(grep "Elapsed (wall clock)" "$OUTPUT_DIR/time.txt" | awk '{print $8}')
SYS=$(grep "System time" "$OUTPUT_DIR/time.txt" | awk '{print $4}')
FAULTS=$(grep "Minor (reclaiming" "$OUTPUT_DIR/time.txt" | awk '{print $6}')
SIZE=$(stat -c%s "$ARCHIVE")

# Baseline from previous default run
BASELINE_MEM=1006208
BASELINE_SYS=79.20
BASELINE_FAULTS=54592515

MEM_DIFF=$(echo "scale=2; (($MEM - $BASELINE_MEM) * 100.0) / $BASELINE_MEM" | bc)
SYS_DIFF=$(echo "scale=2; (($SYS - $BASELINE_SYS) * 100.0) / $BASELINE_SYS" | bc)
FAULT_DIFF=$(echo "scale=2; (($FAULTS - $BASELINE_FAULTS) * 100.0) / $BASELINE_FAULTS" | bc)

printf "%-20s %15s %15s %12s\n" "Metric" "Baseline" "Optimized" "Change"
printf "%s\n" "$(printf '=%.0s' {1..68})"
printf "%-20s %13s KB %13s KB %11s%%\n" "Peak Memory" "$BASELINE_MEM" "$MEM" "$MEM_DIFF"
printf "%-20s %14ss %14ss %11s%%\n" "System Time" "$BASELINE_SYS" "$SYS" "$SYS_DIFF"
printf "%-20s %15s %15s %11s%%\n" "Page Faults" "$BASELINE_FAULTS" "$FAULTS" "$FAULT_DIFF"
printf "%-20s %15s %15s\n" "Wall Time" "~15.1s" "$TIME"
printf "%-20s %15s\n" "Archive Size" "$SIZE bytes"

echo ""
echo "Negative percentage = improvement"
echo ""
echo "Results saved to: $OUTPUT_DIR"
