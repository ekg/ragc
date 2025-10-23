# RAGC Memory Profiling Results

## Summary

Comparative memory profiling of RAGC vs C++ AGC on yeast10 dataset (10 samples, ~221MB input).

### Key Findings

| Metric | RAGC | C++ AGC | Ratio (RAGC/C++) |
|--------|------|---------|------------------|
| **Peak Memory** | 1,008,300 KB (~984 MB) | 209,584 KB (~205 MB) | **4.81x** |
| **Wall Time** | 15.20s | 2.98s | **5.10x** |
| **Archive Size** | 9,328,016 bytes (~8.9 MB) | 6,207,610 bytes (~5.9 MB) | **1.50x** |

## Test Configuration

- **Dataset**: yeast10 (10 yeast samples)
- **Input Size**: ~221 MB (10 FASTA files)
- **Parameters**:
  - K-mer length: 21
  - Segment size: 10,000
  - Min match length: 20
  - Threads: 1 (single-threaded for fair comparison)
  - Compression level: RAGC default (11), C++ AGC default

## Detailed Results

### RAGC Performance
```
Command: ragc create -o archive.agc -k 21 -s 10000 -m 20 -v 0 [10 input files]
User time: 16.67s
System time: 74.48s
CPU usage: 599%
Wall time: 15.20s
Peak memory: 1,008,300 KB
Page faults (minor): 54,592,079
Voluntary context switches: 14,825
Involuntary context switches: 4,990
Output size: 9,328,016 bytes
```

### C++ AGC Performance
```
Command: agc create -o archive.agc -k 21 -s 10000 -l 20 -t 1 -v 0 [10 input files]
User time: 2.95s
System time: 0.09s
CPU usage: 101%
Wall time: 2.98s
Peak memory: 209,584 KB
Page faults (minor): 83,984
Voluntary context switches: 29
Involuntary context switches: 58
Output size: 6,207,610 bytes
```

## Analysis

### Memory Usage (4.81x Gap)

RAGC uses **4.81x more memory** than C++ AGC. Potential causes:

1. **Parallel Processing Overhead**: RAGC uses 599% CPU (6 cores) vs C++ 101% (1 core)
   - Multiple thread-local buffers
   - Crossbeam channel buffering
   - Rayon work-stealing overhead

2. **High System Time**: 74.48s system time vs 0.09s for C++
   - Suggests excessive memory allocations/deallocations
   - 54M page faults vs 84K for C++

3. **Known Memory-Intensive Operations**:
   - K-mer collection: ~91 MB for 11.5M k-mers (documented in debug output)
   - Segment buffering across multiple threads
   - HashMap/BTreeMap usage for group tracking

### Performance (5.10x Slower)

RAGC takes **5.10x longer** despite using 6 cores:

1. **System Time Dominates**: 74.48s system / 16.67s user = 4.47:1 ratio
   - C++ AGC: 0.09s system / 2.95s user = 0.03:1 ratio
   - **High system time indicates memory management bottleneck**

2. **Context Switching**: 14,825 voluntary + 4,990 involuntary
   - Suggests thread contention or I/O waits

### Compression Ratio (1.50x Larger Archives)

RAGC produces **1.50x larger archives**:

1. **Possible Causes**:
   - Different default compression levels
   - Segmentation differences
   - Encoding inefficiencies
   - Metadata overhead

## Optimization Opportunities

### High Priority (Memory)

1. **Reduce Parallel Buffering**
   - Profile per-thread memory usage
   - Consider streaming/flushing strategies
   - Investigate group buffer sizes

2. **Memory Allocator**
   - Try jemalloc or mimalloc
   - Profile allocation patterns
   - Consider object pooling for hot paths

3. **K-mer Collection**
   - 91 MB for 11.5M k-mers is ~8 bytes/k-mer
   - Consider compact representations
   - Streaming k-mer processing

### Medium Priority (Performance)

1. **System Time Investigation**
   - Profile allocation hotspots
   - Reduce memory churn
   - Consider zero-copy strategies

2. **Thread Coordination**
   - Analyze crossbeam channel usage
   - Review Rayon thread pool configuration
   - Reduce lock contention

### Low Priority (Compression)

1. **Compression Level Tuning**
   - Test different ZSTD levels
   - Compare segmentation strategies
   - Validate encoding correctness

## Next Steps

1. ✅ **Baseline established**: RAGC uses 4.81x memory, 5.10x time, 1.50x archive size
2. 🔄 **Deep profiling**: Use valgrind massif for heap profiling
3. 🔄 **Allocator swap**: Test jemalloc/mimalloc impact
4. 🔄 **Algorithmic review**: Compare k-mer collection with C++ implementation
5. 🔄 **Parallel tuning**: Experiment with thread counts and buffer sizes

## Raw Data

Full profiling output saved to: `/tmp/memory_profile_20251022_212350/`

- `ragc_time.txt`: RAGC /usr/bin/time output
- `cpp_time.txt`: C++ AGC /usr/bin/time output
- `summary.txt`: Side-by-side comparison
- `ragc_archive.agc`: RAGC output archive
- `cpp_archive.agc`: C++ AGC output archive
