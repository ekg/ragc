#!/bin/bash
#
# Test different thread counts to find optimal parallelism
# Compare 1, 2, 4, and 6 threads
#
# Hypothesis: Single-thread may be FASTER than multi-thread due to:
# - No channel coordination overhead
# - No crossbeam internal allocations
# - No rayon work-stealing overhead
# - Better cache locality

set -e

TEST_DIR="/home/erik/scrapy/yeast10_test"
OUTPUT_DIR="/tmp/thread_test_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$OUTPUT_DIR"

echo "=== Thread Count Performance Testing ==="
echo "Testing with: 1, 2, 4, 6 threads"
echo "Dataset: yeast10 (~221 MB, 10 samples)"
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

# Function to run compression with specific thread count
run_test() {
    local threads=$1
    local output_file="$OUTPUT_DIR/threads_${threads}.txt"
    local archive_file="$OUTPUT_DIR/archive_${threads}.agc"

    echo "=== Testing with $threads thread(s) ==="

    /usr/bin/time -v /home/erik/ragc/target/release/ragc create \
        -o "$archive_file" \
        -k 21 \
        -s 10000 \
        -m 20 \
        -v 0 \
        -t "$threads" \
        "$SAMPLES_DIR"/*.fa \
        2>&1 | tee "$output_file"

    echo ""
}

# Run tests with different thread counts
run_test 1
run_test 2
run_test 4
run_test 6

# Extract and compare results
echo "=== RESULTS SUMMARY ==="
echo ""

printf "%-10s %15s %12s %12s %15s %15s\n" "Threads" "Memory (KB)" "Wall Time" "System (s)" "User (s)" "Page Faults"
printf "%s\n" "$(printf '=%.0s' {1..85})"

for threads in 1 2 4 6; do
    file="$OUTPUT_DIR/threads_${threads}.txt"

    mem=$(grep "Maximum resident set size" "$file" | awk '{print $6}')
    wall=$(grep "Elapsed (wall clock)" "$file" | awk '{print $8}')
    sys=$(grep "System time" "$file" | awk '{print $4}')
    user=$(grep "User time" "$file" | awk '{print $4}')
    faults=$(grep "Minor (reclaiming" "$file" | awk '{print $6}')

    printf "%-10s %15s %12s %12ss %15ss %15s\n" "$threads" "$mem" "$wall" "$sys" "$user" "$faults"
done

echo ""
echo "=== ANALYSIS ==="
echo ""

# Get baseline (6 threads)
baseline_mem=$(grep "Maximum resident set size" "$OUTPUT_DIR/threads_6.txt" | awk '{print $6}')
baseline_sys=$(grep "System time" "$OUTPUT_DIR/threads_6.txt" | awk '{print $4}')
baseline_faults=$(grep "Minor (reclaiming" "$OUTPUT_DIR/threads_6.txt" | awk '{print $6}')

# Compare single-thread
single_mem=$(grep "Maximum resident set size" "$OUTPUT_DIR/threads_1.txt" | awk '{print $6}')
single_sys=$(grep "System time" "$OUTPUT_DIR/threads_1.txt" | awk '{print $4}')
single_faults=$(grep "Minor (reclaiming" "$OUTPUT_DIR/threads_1.txt" | awk '{print $6}')

mem_diff=$(echo "scale=2; (($single_mem - $baseline_mem) * 100.0) / $baseline_mem" | bc)
sys_diff=$(echo "scale=2; (($single_sys - $baseline_sys) * 100.0) / $baseline_sys" | bc)
fault_diff=$(echo "scale=2; (($single_faults - $baseline_faults) * 100.0) / $baseline_faults" | bc)

echo "Single-thread (1) vs Baseline (6 threads):"
printf "  Memory:      %s KB vs %s KB (%+.1f%%)\n" "$single_mem" "$baseline_mem" "$mem_diff"
printf "  System time: %s s vs %s s (%+.1f%%)\n" "$single_sys" "$baseline_sys" "$sys_diff"
printf "  Page faults: %s vs %s (%+.1f%%)\n" "$single_faults" "$baseline_faults" "$fault_diff"
echo ""

# Compare C++ AGC (for reference)
echo "For reference, C++ AGC (1 thread):"
echo "  Memory:      205584 KB"
echo "  System time: 0.09 s"
echo "  Page faults: 83984"
echo "  Wall time:   ~3.0 s"
echo ""

echo "Key Insights:"
echo "  - Negative % = improvement (lower is better)"
echo "  - If single-thread is faster: parallelism overhead dominates"
echo "  - If multi-thread is faster: parallelism helps despite overhead"
echo ""
echo "Results saved to: $OUTPUT_DIR"
