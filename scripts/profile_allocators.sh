#!/bin/bash
#
# Compare memory allocators: system default vs jemalloc
# Tests both on the same dataset to measure impact

set -e

TEST_DIR="/home/erik/scrapy/yeast10_test"
OUTPUT_DIR="/tmp/allocator_profile_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$OUTPUT_DIR"

echo "=== Allocator Comparison: Default vs Jemalloc ==="
echo "Test data: $TEST_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Prepare sample files once
SAMPLES_DIR="$OUTPUT_DIR/samples"
mkdir -p "$SAMPLES_DIR"
for f in "$TEST_DIR"/*.fa; do
    basename_f=$(basename "$f")
    if [ "$basename_f" != "yeast10_pansn.fa" ]; then
        cp "$f" "$SAMPLES_DIR/"
    fi
done

echo "=== Building RAGC (default allocator) ==="
cd /home/erik/ragc
cargo build --release --quiet
echo "✓ Default build complete"
echo ""

echo "=== Profiling with DEFAULT allocator ==="
DEFAULT_ARCHIVE="$OUTPUT_DIR/default_archive.agc"
/usr/bin/time -v /home/erik/ragc/target/release/ragc create \
    -o "$DEFAULT_ARCHIVE" \
    -k 21 \
    -s 10000 \
    -m 20 \
    -v 0 \
    "$SAMPLES_DIR"/*.fa \
    2>&1 | tee "$OUTPUT_DIR/default_time.txt"

DEFAULT_SIZE=$(stat -c%s "$DEFAULT_ARCHIVE")
DEFAULT_MEM=$(grep "Maximum resident set size" "$OUTPUT_DIR/default_time.txt" | awk '{print $6}')
DEFAULT_TIME=$(grep "Elapsed (wall clock)" "$OUTPUT_DIR/default_time.txt" | awk '{print $8}')
DEFAULT_SYS=$(grep "System time" "$OUTPUT_DIR/default_time.txt" | awk '{print $4}')
DEFAULT_FAULTS=$(grep "Minor (reclaiming" "$OUTPUT_DIR/default_time.txt" | awk '{print $6}')

echo ""
echo "=== Building RAGC with JEMALLOC ==="
cargo build --release --features jemalloc --quiet
echo "✓ Jemalloc build complete"
echo ""

echo "=== Profiling with JEMALLOC ==="
JEMALLOC_ARCHIVE="$OUTPUT_DIR/jemalloc_archive.agc"
/usr/bin/time -v /home/erik/ragc/target/release/ragc create \
    -o "$JEMALLOC_ARCHIVE" \
    -k 21 \
    -s 10000 \
    -m 20 \
    -v 0 \
    "$SAMPLES_DIR"/*.fa \
    2>&1 | tee "$OUTPUT_DIR/jemalloc_time.txt"

JEMALLOC_SIZE=$(stat -c%s "$JEMALLOC_ARCHIVE")
JEMALLOC_MEM=$(grep "Maximum resident set size" "$OUTPUT_DIR/jemalloc_time.txt" | awk '{print $6}')
JEMALLOC_TIME=$(grep "Elapsed (wall clock)" "$OUTPUT_DIR/jemalloc_time.txt" | awk '{print $8}')
JEMALLOC_SYS=$(grep "System time" "$OUTPUT_DIR/jemalloc_time.txt" | awk '{print $4}')
JEMALLOC_FAULTS=$(grep "Minor (reclaiming" "$OUTPUT_DIR/jemalloc_time.txt" | awk '{print $6}')

echo ""
echo "=== COMPARISON RESULTS ===" | tee "$OUTPUT_DIR/comparison.txt"
echo "" | tee -a "$OUTPUT_DIR/comparison.txt"

printf "%-20s %15s %15s %10s\n" "Metric" "Default" "Jemalloc" "Change" | tee -a "$OUTPUT_DIR/comparison.txt"
printf "%s\n" "$(printf '=%.0s' {1..65})" | tee -a "$OUTPUT_DIR/comparison.txt"

# Memory
MEM_DIFF=$(echo "scale=2; (($JEMALLOC_MEM - $DEFAULT_MEM) * 100.0) / $DEFAULT_MEM" | bc)
printf "%-20s %13s KB %13s KB %9s%%\n" "Peak Memory" "$DEFAULT_MEM" "$JEMALLOC_MEM" "$MEM_DIFF" | tee -a "$OUTPUT_DIR/comparison.txt"

# System Time
SYS_DIFF=$(echo "scale=2; (($JEMALLOC_SYS - $DEFAULT_SYS) * 100.0) / $DEFAULT_SYS" | bc)
printf "%-20s %14ss %14ss %9s%%\n" "System Time" "$DEFAULT_SYS" "$JEMALLOC_SYS" "$SYS_DIFF" | tee -a "$OUTPUT_DIR/comparison.txt"

# Wall Time
printf "%-20s %15s %15s\n" "Wall Time" "$DEFAULT_TIME" "$JEMALLOC_TIME" | tee -a "$OUTPUT_DIR/comparison.txt"

# Page Faults
FAULT_DIFF=$(echo "scale=2; (($JEMALLOC_FAULTS - $DEFAULT_FAULTS) * 100.0) / $DEFAULT_FAULTS" | bc)
printf "%-20s %15s %15s %9s%%\n" "Page Faults" "$DEFAULT_FAULTS" "$JEMALLOC_FAULTS" "$FAULT_DIFF" | tee -a "$OUTPUT_DIR/comparison.txt"

# Archive Size
SIZE_DIFF=$(echo "scale=2; (($JEMALLOC_SIZE - $DEFAULT_SIZE) * 100.0) / $DEFAULT_SIZE" | bc)
printf "%-20s %15s %15s %9s%%\n" "Archive Size" "$DEFAULT_SIZE" "$JEMALLOC_SIZE" "$SIZE_DIFF" | tee -a "$OUTPUT_DIR/comparison.txt"

echo "" | tee -a "$OUTPUT_DIR/comparison.txt"
echo "Negative percentage = improvement (lower is better)" | tee -a "$OUTPUT_DIR/comparison.txt"
echo "" | tee -a "$OUTPUT_DIR/comparison.txt"
echo "Full results saved to: $OUTPUT_DIR" | tee -a "$OUTPUT_DIR/comparison.txt"
