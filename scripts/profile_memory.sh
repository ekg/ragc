#!/bin/bash
#
# Memory profiling script comparing RAGC vs C++ AGC
# Uses /usr/bin/time -v to measure peak memory usage

set -e

TEST_DIR="/home/erik/scrapy/yeast10_test"
PANSN_FILE="$TEST_DIR/yeast10_pansn.fa"
OUTPUT_DIR="/tmp/memory_profile_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$OUTPUT_DIR"

echo "=== Memory Profiling: RAGC vs C++ AGC ==="
echo "Test data: $PANSN_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Build RAGC in release mode
echo "Building RAGC (release mode)..."
cd /home/erik/ragc
cargo build --release --quiet
echo "✓ RAGC built"
echo ""

# Build C++ AGC
echo "Building C++ AGC..."
cd /home/erik/agc
make clean >/dev/null 2>&1 || true
make -j$(nproc) >/dev/null 2>&1
echo "✓ C++ AGC built"
echo ""

# Profile RAGC compression
echo "=== Profiling RAGC Compression ==="
RAGC_ARCHIVE="$OUTPUT_DIR/ragc_archive.agc"

# RAGC needs individual sample files
echo "Preparing individual sample files for RAGC..."
RAGC_SAMPLES_DIR="$OUTPUT_DIR/ragc_samples"
mkdir -p "$RAGC_SAMPLES_DIR"

# Use the individual sample files from test directory
for f in "$TEST_DIR"/*.fa; do
    basename_f=$(basename "$f")
    if [ "$basename_f" != "yeast10_pansn.fa" ]; then
        cp "$f" "$RAGC_SAMPLES_DIR/"
    fi
done

/usr/bin/time -v /home/erik/ragc/target/release/ragc create \
    -o "$RAGC_ARCHIVE" \
    -k 21 \
    -s 10000 \
    -m 20 \
    -v 0 \
    "$RAGC_SAMPLES_DIR"/*.fa \
    2>&1 | tee "$OUTPUT_DIR/ragc_time.txt"

RAGC_SIZE=$(stat -f%z "$RAGC_ARCHIVE" 2>/dev/null || stat -c%s "$RAGC_ARCHIVE")
echo ""
echo "RAGC archive size: $RAGC_SIZE bytes"
echo ""

# Profile C++ AGC compression
echo "=== Profiling C++ AGC Compression ==="
CPP_ARCHIVE="$OUTPUT_DIR/cpp_archive.agc"

# C++ AGC uses the same sample files we prepared for RAGC
# C++ AGC syntax: agc create -o output.agc -k 21 -s 10000 -l 20 -t 1 input1.fa input2.fa...
/usr/bin/time -v /home/erik/agc/bin/agc create -o "$CPP_ARCHIVE" -k 21 -s 10000 -l 20 -t 1 -v 0 "$RAGC_SAMPLES_DIR"/*.fa \
    2>&1 | tee "$OUTPUT_DIR/cpp_time.txt"

CPP_SIZE=$(stat -f%z "$CPP_ARCHIVE" 2>/dev/null || stat -c%s "$CPP_ARCHIVE")
echo ""
echo "C++ AGC archive size: $CPP_SIZE bytes"
echo ""

# Extract memory statistics
echo "=== Memory Usage Summary ===" | tee "$OUTPUT_DIR/summary.txt"
echo "" | tee -a "$OUTPUT_DIR/summary.txt"

RAGC_MEM=$(grep "Maximum resident set size" "$OUTPUT_DIR/ragc_time.txt" | awk '{print $6}')
CPP_MEM=$(grep "Maximum resident set size" "$OUTPUT_DIR/cpp_time.txt" | awk '{print $6}')

RAGC_TIME=$(grep "Elapsed (wall clock)" "$OUTPUT_DIR/ragc_time.txt" | awk '{print $8}')
CPP_TIME=$(grep "Elapsed (wall clock)" "$OUTPUT_DIR/cpp_time.txt" | awk '{print $8}')

echo "RAGC:" | tee -a "$OUTPUT_DIR/summary.txt"
echo "  Peak Memory: ${RAGC_MEM} KB" | tee -a "$OUTPUT_DIR/summary.txt"
echo "  Wall Time: ${RAGC_TIME}" | tee -a "$OUTPUT_DIR/summary.txt"
echo "  Archive Size: ${RAGC_SIZE} bytes" | tee -a "$OUTPUT_DIR/summary.txt"
echo "" | tee -a "$OUTPUT_DIR/summary.txt"

echo "C++ AGC:" | tee -a "$OUTPUT_DIR/summary.txt"
echo "  Peak Memory: ${CPP_MEM} KB" | tee -a "$OUTPUT_DIR/summary.txt"
echo "  Wall Time: ${CPP_TIME}" | tee -a "$OUTPUT_DIR/summary.txt"
echo "  Archive Size: ${CPP_SIZE} bytes" | tee -a "$OUTPUT_DIR/summary.txt"
echo "" | tee -a "$OUTPUT_DIR/summary.txt"

# Calculate ratios
if [ "$CPP_MEM" -gt 0 ]; then
    MEM_RATIO=$(echo "scale=2; $RAGC_MEM / $CPP_MEM" | bc)
    echo "Memory Ratio (RAGC/C++): ${MEM_RATIO}x" | tee -a "$OUTPUT_DIR/summary.txt"
fi

if [ "$CPP_SIZE" -gt 0 ]; then
    SIZE_RATIO=$(echo "scale=2; $RAGC_SIZE / $CPP_SIZE" | bc)
    echo "Archive Size Ratio (RAGC/C++): ${SIZE_RATIO}x" | tee -a "$OUTPUT_DIR/summary.txt"
fi

echo "" | tee -a "$OUTPUT_DIR/summary.txt"
echo "Full results saved to: $OUTPUT_DIR" | tee -a "$OUTPUT_DIR/summary.txt"
