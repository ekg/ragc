#!/bin/bash
#
# Test ZSTD buffer pooling optimization
# Compare against baseline to measure impact

set -e

TEST_DIR="/home/erik/scrapy/yeast10_test"
OUTPUT_DIR="/tmp/zstd_pool_test_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$OUTPUT_DIR"

echo "=== Testing ZSTD Buffer Pooling Optimization ==="
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

# Test with buffer pooling
echo "Running compression with ZSTD buffer pooling..."
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
USER=$(grep "User time" "$OUTPUT_DIR/time.txt" | awk '{print $4}')
FAULTS=$(grep "Minor (reclaiming" "$OUTPUT_DIR/time.txt" | awk '{print $6}')
SIZE=$(stat -c%s "$ARCHIVE")

# Baseline from previous run (before optimization)
BASELINE_MEM=1006208
BASELINE_SYS=79.20
BASELINE_FAULTS=54592515
BASELINE_TIME="15.12s"

MEM_DIFF=$(echo "scale=2; (($MEM - $BASELINE_MEM) * 100.0) / $BASELINE_MEM" | bc)
SYS_DIFF=$(echo "scale=2; (($SYS - $BASELINE_SYS) * 100.0) / $BASELINE_SYS" | bc)
FAULT_DIFF=$(echo "scale=2; (($FAULTS - $BASELINE_FAULTS) * 100.0) / $BASELINE_FAULTS" | bc)

printf "%-25s %15s %15s %12s\n" "Metric" "Baseline" "With Pooling" "Change"
printf "%s\n" "$(printf '=%.0s' {1..75})"
printf "%-25s %13s KB %13s KB %11s%%\n" "Peak Memory" "$BASELINE_MEM" "$MEM" "$MEM_DIFF"
printf "%-25s %14ss %14ss %11s%%\n" "System Time" "$BASELINE_SYS" "$SYS" "$SYS_DIFF"
printf "%-25s %14ss %14ss\n" "User Time" "16.67s" "$USER"
printf "%-25s %15s %15s %11s%%\n" "Page Faults" "$BASELINE_FAULTS" "$FAULTS" "$FAULT_DIFF"
printf "%-25s %15s %15s\n" "Wall Time" "$BASELINE_TIME" "$TIME"
printf "%-25s %15s\n" "Archive Size" "$SIZE bytes"

echo ""
echo "Negative percentage = improvement (lower is better)"
echo ""
echo "Expected improvements from buffer pooling:"
echo "  - System time: -20 to -30%"
echo "  - Page faults: -30 to -50%"
echo "  - Memory: -5 to -10%"
echo ""
echo "Results saved to: $OUTPUT_DIR"
