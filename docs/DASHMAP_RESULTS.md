# DashMap Parallelization Fix Results

**Date**: 2025-10-23
**Fix**: Replaced `Arc<Mutex<HashMap>>` with `DashMap` for concurrent group access
**Dataset**: yeast235 (37 files, 644 MB)

---

## Summary

**Progress**: DashMap **doubled CPU utilization** and gave **1.4x speedup**!

**Limitation**: Still **20x slower than C++ AGC** - bottleneck shifted from locking to system time.

---

## Benchmark Results

### RAGC Performance

| Configuration | Wall Time | CPU % | User | System | Memory | Speedup |
|---------------|-----------|-------|------|--------|--------|---------|
| **Old (Mutex) 1T** | 127s | 101% | 47.7s | 81.4s | 918 MB | 1.0x |
| **Old (Mutex) 6T** | 136s | 111% | 56.7s | 95.2s | 1003 MB | **0.93x** |
| **DashMap 6T** | **91s** | **217%** | 67.7s | 129.8s | 1168 MB | **1.4x** ✅ |
| **DashMap 12T** | 91s | **256%** | 80.0s | 153.6s | 1403 MB | **1.4x** |

### C++ AGC Performance (Reference)

| Threads | Wall Time | CPU % | User | System | Memory | Speedup |
|---------|-----------|-------|------|--------|--------|---------|
| 1 | 23s | 101% | 23.2s | 0.4s | 968 MB | 1.0x |
| 6 | 6.2s | 516% | 31.5s | 0.5s | 947 MB | **3.7x** |
| 12 | 4.8s | 917% | 43.3s | 0.8s | 928 MB | **4.8x** |

---

## Analysis

### What DashMap Fixed

**Before (Arc<Mutex<HashMap>>)**:
- Global exclusive lock on entire HashMap
- Every segment blocked all other threads
- CPU utilization: 111% (only ~1 core)
- Speedup: 0.93x (NEGATIVE - slower with more threads!)

**After (DashMap)**:
- Per-shard locking (fine-grained concurrency)
- Threads can access different groups simultaneously
- CPU utilization: 217-256% (~2-2.5 cores)
- Speedup: **1.4x** (actual parallelism!)

**Improvement**: **2x CPU utilization**, **1.5x speedup**

### Remaining Bottleneck: System Time

| Version | User Time | System Time | Ratio |
|---------|-----------|-------------|-------|
| C++ AGC (6T) | 31.5s | 0.5s | **63:1** (user dominates) |
| RAGC DashMap (6T) | 67.7s | 129.8s | **1:2** (system dominates!) |

**Problem**: RAGC spends **2x more time in system calls** than actual computation!

System time sources:
1. **Memory allocation churn**: 120M page faults vs C++ AGC's ~5M
2. **I/O overhead**: Reading FASTAs, writing archives
3. **Memory management**: Rust allocator overhead

---

## CPU Utilization Analysis

### DashMap RAGC (6 threads)
```
User + System = 67.7s + 129.8s = 197.5s total CPU
Wall time = 90.7s
CPU % = 217% ≈ 2.2 cores utilized

Expected with perfect parallelism: 600% (6 cores)
Actual: 217%
Efficiency: 36%
```

### C++ AGC (6 threads)
```
User + System = 31.5s + 0.5s = 32s total CPU
Wall time = 6.2s
CPU % = 516% ≈ 5.2 cores utilized

Expected: 600%
Actual: 516%
Efficiency: 86%
```

**C++ AGC achieves 86% parallel efficiency vs RAGC's 36%**

---

## Why Didn't We Get Full Parallelism?

DashMap reduced lock contention, but threads are still bottlenecked by:

### 1. **High System Time** (129s vs 0.5s)

Likely causes:
- **Memory allocations**: 49M Vec growth calls in LZDiff (from heaptrack)
- **Page faults**: 120M faults requiring kernel intervention
- **I/O blocking**: FASTA reading, archive writing

### 2. **I/O Serialization**

Worker threads may be blocking on:
- Reader thread (FASTA parsing)
- Writer thread (archive writes)
- Bounded channels (backpressure)

### 3. **Remaining Contention**

DashMap uses per-shard locking, but:
- Default shard count may be too low
- Hot groups may still contend
- Writer thread serializes all writes

---

## Comparison with C++ AGC

### Speed

| Metric | RAGC DashMap (6T) | C++ AGC (6T) | Ratio |
|--------|-------------------|--------------|-------|
| **Wall Time** | 91s | 6.2s | **15x slower** |
| **CPU Time** | 197.5s | 32s | **6x more CPU** |
| **System Time** | 129.8s | 0.5s | **260x more system calls!** |

### Memory

| Metric | RAGC DashMap (6T) | C++ AGC (6T) | Ratio |
|--------|-------------------|--------------|-------|
| **Peak RSS** | 1168 MB | 947 MB | 1.23x more |
| **Page Faults** | 120M | ~5M | **24x more faults** |

---

## What This Tells Us

### Good News ✅

1. **DashMap works**: Doubled CPU utilization, eliminated global lock
2. **Correctness maintained**: All tests pass, C++ AGC can read archives
3. **Memory reasonable**: 1168 MB vs 947 MB C++ AGC (23% more)

### Bad News ❌

1. **System time dominates**: 2x system vs user (should be <<1x)
2. **Still 15x slower than C++ AGC** in wall time
3. **Only 36% parallel efficiency** (vs C++ AGC's 86%)
4. **120M page faults** suggests memory management issues

---

## Next Steps

To achieve C++ AGC-level performance, need to address **system time bottleneck**:

### Option 1: Reduce Memory Allocation Churn

**Problem**: 49M Vec growth calls in LZDiff::prepare

**Solutions**:
- Pre-allocate Vec capacities based on expected size
- Pool buffers per thread
- Use arena allocator for related allocations

**Expected**: -50% system time → ~45s wall, 300% CPU

### Option 2: Try Alternative Allocator

**Problem**: Rust's default allocator may have higher overhead than C++'s

**Solutions**:
- jemalloc: Thread-local caching, better scalability
- mimalloc: Microsoft's fast allocator
- Add `#[global_allocator]` and benchmark

**Expected**: -20-30% system time → ~63s wall

### Option 3: Reduce Page Faults

**Problem**: 120M page faults vs C++ AGC's 5M

**Solutions**:
- Use hugepages for large allocations
- Reduce memory churn (see Option 1)
- Better memory locality

**Expected**: -20% system time → ~72s wall

### Option 4: Optimize I/O Pipeline

**Problem**: Serialized reading/writing may block workers

**Solutions**:
- Larger bounded queues (reduce backpressure stalls)
- Async I/O for FASTA reading
- Multiple writer threads

**Expected**: -15% system time → ~77s wall

### Option 5: Profile with perf

**Need**: Detailed breakdown of where system time is spent

**Command**:
```bash
perf record -g ./target/release/ragc create ...
perf report
```

Identify hotspots in system calls (mmap, brk, futex, etc.)

---

## Implementation Changes

**Files Modified**:
- `ragc-core/src/compressor_streaming.rs`:
  - Line 2365: Changed `Arc<Mutex<HashMap>>` to `Arc<DashMap>`
  - Lines 2545-2558: Use DashMap entry API instead of lock()
  - Lines 2640-2644: Use DashMap::into_iter() instead of into_inner()

**Dependencies**: DashMap already in `Cargo.toml` (workspace dependency)

**Code diff**: ~10 lines changed

---

## Conclusion

**DashMap successfully eliminated the global lock bottleneck**, doubling CPU utilization from 111% to 217-256%. This proves the original diagnosis was correct - `Arc<Mutex<HashMap>>` was serializing all segment operations.

However, **the performance gap to C++ AGC remains large** (15x slower) due to a new bottleneck: **excessive system time** (129s vs 0.5s). This suggests memory management and I/O issues that require deeper optimization.

**Recommended next steps**:
1. Profile with `perf` to identify system time sources
2. Reduce memory allocation churn in LZDiff
3. Try jemalloc or mimalloc allocator
4. Optimize I/O pipeline

DashMap was the right fix for the parallelization bug, but RAGC needs further work to match C++ AGC's efficiency.

---

## Appendix: Test Commands

```bash
# Build release
cargo build --release

# Test correctness
cargo test --release

# Benchmark 6 threads
/usr/bin/time -v ./target/release/ragc create \
  -o /tmp/test.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa

# Verify C++ AGC compatibility
agc listset /tmp/test.agc
agc getset /tmp/test.agc sample1 > /tmp/extracted.fa
```
