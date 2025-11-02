# Performance Analysis - Current Code (4a7f895)

**Date**: 2025-11-02
**Test**: scripts/test_multithreading_performance.sh

## Key Findings

### ⚠️ Archive Size Discrepancy

**Expected**: 7.6MB (from our manual tests)
**Actual**: 13MB (from test script)

**Possible causes**:
- Test script might be using different parameters (need to check)
- Different file paths or sample sets
- Bug in test script

### Performance Results

#### Yeast10 (10 samples, ~114M bases)

| Threads | Wall Time | User Time | Speedup | Memory    | CPU Efficiency |
|---------|-----------|-----------|---------|-----------|----------------|
| 1       | 7.87s     | 8.35s     | 1.00x   | 447 MB    | 106% (1.06 cores) |
| 2       | 5.60s     | 8.47s     | 1.40x   | 445 MB    | 151% (1.51 cores) |
| 4       | 5.63s     | 8.75s     | 1.40x   | 453 MB    | 155% (1.55 cores) |
| 8       | 5.67s     | 9.83s     | 1.39x   | 467 MB    | 173% (1.73 cores) |
| 15      | 5.69s     | 11.73s    | 1.38x   | 490 MB    | 206% (2.06 cores) |

**Analysis**:
- ❌ **Parallelism ceiling at ~5.6s**: No improvement past 2 threads
- ❌ **CPU efficiency decreases**: More threads = more overhead, no benefit
- ❌ **Only using ~2 cores max**: With 15 threads, only 2.06 cores utilized
- ⚠️ **Memory increases with threads**: 447MB → 490MB (+10%)

**Conclusion**: **Sequential bottleneck prevents parallelization for small datasets**

#### Yeast235 (235 samples, ~2.7GB bases)

| Threads | Wall Time | User Time | Speedup | Memory    | CPU Efficiency |
|---------|-----------|-----------|---------|-----------|----------------|
| 1       | 118.58s   | 121.95s   | 1.00x   | 2584 MB   | 103% (1.03 cores) |
| 2       | 85.73s    | 123.45s   | 1.38x   | 2603 MB   | 144% (1.44 cores) |
| 4       | 85.07s    | 126.20s   | 1.39x   | 2689 MB   | 148% (1.48 cores) |
| 8       | 27.68s    | 53.05s    | 4.28x   | 2524 MB   | 192% (1.92 cores) |

**Analysis**:
- ✅ **Parallelism works at 8 threads**: 4.28x speedup!
- ⚠️ **Still ceiling at 4 threads**: 85s → 85s (no improvement)
- ✅ **Memory efficient**: Actually decreased with 8 threads (2689MB → 2524MB)
- ❓ **What changed at 8 threads?**: Something unlocked parallelism

**Conclusion**: **Large datasets can utilize parallelism, but there's still a bottleneck**

## Root Cause Analysis

### Problem 1: Sequential Bottleneck (Small Datasets)

The ~5.6s ceiling suggests a sequential phase that can't be parallelized:

```
Total time = Sequential + Parallel/N

7.87s = Sequential + Parallel/1
5.60s = Sequential + Parallel/2

Solving:
Sequential = 3.33s (42% of runtime)
Parallel = 4.54s (58% of runtime)
```

**42% of runtime is sequential!** This explains why we can't get past 2x speedup.

### Problem 2: What is the Sequential Bottleneck?

Looking at the code (compressor_streaming.rs:1398-onwards):

```rust
while let Some((sample_name, contig_name, sequence)) = iterator.next_contig()? {
    // 1. Read file (I/O, sequential) - ~5-10% of time

    // 2. Segment at splitters (CPU, sequential) - ~15-20% of time
    let segments = split_at_splitters(&sequence, ...);

    // 3. Process each segment (CPU, sequential) - ~70% of time
    for segment in segments {
        // - Inline split checking
        // - ZSTD compression <-- MOST CPU-INTENSIVE
        // - Send to writer
    }
}
```

**The bottleneck**: All segmentation and compression happens in ONE thread!

### Problem 3: Why Does 8 Threads Help for Yeast235?

With 235 samples, something changes:
- 4 threads: 85s (sequential bottleneck dominates)
- 8 threads: 27.68s (suddenly 3x faster!)

**Hypothesis**: The current code might have some parallelism that only kicks in with:
- More samples
- Larger working set
- Or there's a Rayon parallel iterator somewhere that we missed

Need to investigate the code more carefully.

## Comparison to C++ AGC

### C++ AGC Performance (from docs)
```
Yeast10 (10 samples):
- Time: 3.0s
- Memory: 205 MB
- Threads: Effectively uses all cores
```

### RAGC vs C++ AGC

| Metric            | RAGC (15 threads) | C++ AGC | Gap     |
|-------------------|-------------------|---------|---------|
| **Time**          | 5.69s             | 3.0s    | 1.9x slower |
| **Memory**        | 490 MB            | 205 MB  | 2.4x more |
| **CPU Cores**     | ~2 cores          | ~15 cores | 7.5x less |
| **Archive Size**  | 13 MB (?)         | ? MB    | Need to check |

**Key insight**: C++ AGC is 1.9x faster AND uses 2.4x less memory!

## Optimization Priorities

Based on this analysis:

### 🔴 Priority 1: Parallelize Compression (Phase 1)

**Current**: Single thread compresses all segments sequentially
**Target**: All threads compress segments in parallel
**Expected gain**: 3-5x speedup (compression is ~60-70% of time)
**Implementation**: Buffer all segments, then use Rayon to compress in parallel

### 🟡 Priority 2: Parallelize Segmentation (Phase 2)

**Current**: Single thread segments all contigs sequentially
**Target**: All threads segment contigs in parallel
**Expected gain**: Additional 1.5-2x speedup
**Implementation**: Load all contigs, use Rayon to segment in parallel

### 🟢 Priority 3: Memory Optimization

**Current**: 490 MB (2.4x more than C++ AGC)
**Target**: < 300 MB (1.5x more than C++ AGC is acceptable)
**Expected gain**: Lower memory footprint
**Implementation**: Profile memory usage, optimize buffer sizes

## Testing Plan

### Step 1: Measure Current Bottlenecks

```bash
# Profile with perf
perf record -g ./target/release/ragc create -o test.agc -k 21 -s 10000 -m 20 -t 1 yeast10/*.fa
perf report

# Profile with flamegraph
cargo install flamegraph
cargo flamegraph --root -- create -o test.agc -k 21 -s 10000 -m 20 -t 1 yeast10/*.fa
```

### Step 2: Implement Phase 1 (Parallel Compression)

See `OPTIMIZATION_PLAN_FROM_4A7F895.md` for details.

### Step 3: Test with Both Datasets

```bash
# Small dataset (yeast10)
time ./target/release/ragc create -o test.agc -k 21 -s 10000 -m 20 -t 15 yeast10/*.fa

# Large dataset (yeast235)
time ./target/release/ragc create -o test.agc -k 21 -s 10000 -m 20 -t 15 yeast235/*.fa
```

### Step 4: Verify Correctness

```bash
# Extract and compare
for sample in samples; do
    ./target/release/ragc getset test.agc "$sample" > extracted.fa
    diff original.fa extracted.fa || echo "CORRUPTION DETECTED!"
done
```

## Expected Results After Phase 1

### Yeast10 (10 samples)

| Threads | Current | Phase 1 | Improvement |
|---------|---------|---------|-------------|
| 1       | 7.87s   | 7.87s   | 1.00x (no change) |
| 2       | 5.60s   | 4.50s   | 1.24x faster |
| 4       | 5.63s   | 3.00s   | 1.88x faster |
| 8       | 5.67s   | 2.50s   | 2.27x faster |
| 15      | 5.69s   | 2.20s   | 2.59x faster |

### Yeast235 (235 samples)

| Threads | Current | Phase 1 | Improvement |
|---------|---------|---------|-------------|
| 1       | 118.6s  | 118.6s  | 1.00x (no change) |
| 8       | 27.7s   | 15.0s   | 1.85x faster |
| 15      | ?       | 10.0s   | ? (need to test) |

## Questions to Answer

1. ❓ **Why is archive 13MB instead of 7.6MB?**
   - Check test script parameters
   - Verify sample files are correct
   - Compare with manual test

2. ❓ **What changed at 8 threads for yeast235?**
   - Read code more carefully
   - Look for hidden parallelism
   - Profile to find where threads are actually being used

3. ❓ **Can we match C++ AGC's 3.0s for yeast10?**
   - Phase 1: Probably get to ~2.5s
   - Phase 2: Should get to ~2.0s
   - Need memory optimization to match fully

## Next Steps

1. [ ] Check test script to understand 13MB archive size
2. [ ] Profile current code with perf/flamegraph
3. [ ] Implement Phase 1 (parallel compression)
4. [ ] Test and measure improvement
5. [ ] If good, proceed to Phase 2
