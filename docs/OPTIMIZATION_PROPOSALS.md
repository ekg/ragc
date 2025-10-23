# RAGC Optimization Proposals

## Profiling Summary

### Baseline Comparison (yeast10 dataset)
| Metric | RAGC | C++ AGC | Ratio |
|--------|------|---------|-------|
| Peak Memory | 984 MB | 205 MB | **4.81x** |
| Wall Time | 15.2s | 3.0s | **5.10x** |
| System Time | 74.5s | 0.09s | **828x** ⚠️ |
| Page Faults | 54.6M | 84K | **650x** ⚠️ |
| Archive Size | 8.9 MB | 5.9 MB | **1.50x** |

### Allocator Testing
- **Jemalloc**: WORSE than default (+7% memory, +16% system time)
- **Conclusion**: Problem is excessive allocations, not allocator efficiency

## Root Causes Identified

### 1. Excessive Memory Allocations (54.6M page faults)
**Evidence:**
- 54.6 million page faults vs 84K for C++ (650x more)
- System time 74.5s vs 0.09s (828x more)
- Indicates constant allocation/deallocation

**Location**: `compressor_streaming.rs` parallel pipeline

**Current Design:**
```rust
// Three bounded channels with 100-item buffers each
bounded(100) // contig_rx
bounded(100) // uncompressed_rx
bounded(100) // compressed_rx
```

**Problem**: Each buffer slot holds:
- UncompressedPack with Vec<u8> (LZ-encoded data, ~10KB+)
- CompressedPack with Vec<u8> (ZSTD data, varies)
- Multiple Vec allocations per pack

**Impact**: With 100 slots × 3 channels × multiple threads, significant buffer overhead

### 2. K-mer Collection (91 MB - acceptable)
**Evidence:**
- 11.5M k-mers × 8 bytes = 91 MB
- Single Vec allocation, then radix sort

**Assessment**: This matches C++ AGC design and is not the primary issue

### 3. Parallel Thread Overhead
**Evidence:**
- 599% CPU usage (6 threads) vs C++'s 101% (1 thread configured)
- Multiple thread-local allocations
- Arc/RwLock coordination overhead

### 4. Compression Ratio (1.50x larger archives)
**Possible Causes:**
- Different ZSTD compression levels
- Segmentation strategy differences
- Encoding inefficiencies
- Metadata overhead

## Proposed Optimizations

### Priority 1: Reduce Channel Buffering

**Proposal 1A: Smaller Channel Buffers**
```rust
// Current
bounded(100)

// Proposed
bounded(10)  // Reduce 10x
```

**Expected Impact:**
- 10x less memory in channel buffers
- May slightly reduce parallelism but reduce memory pressure
- Easy to test

**Implementation**: `compressor_streaming.rs:1018,1024,1027`

---

**Proposal 1B: Object Pooling for Packs**
```rust
// Reuse Vec<u8> buffers instead of allocating new ones
struct PackPool {
    buffers: Vec<Vec<u8>>,
}

impl PackPool {
    fn get_buffer(&mut self, capacity: usize) -> Vec<u8> {
        self.buffers.pop()
            .map(|mut buf| { buf.clear(); buf.reserve(capacity); buf })
            .unwrap_or_else(|| Vec::with_capacity(capacity))
    }

    fn return_buffer(&mut self, buf: Vec<u8>) {
        if buf.capacity() < 1_000_000 {  // Don't keep huge buffers
            self.buffers.push(buf);
        }
    }
}
```

**Expected Impact:**
- Reduce allocation count dramatically
- Lower page faults by reusing memory
- More complex implementation

---

### Priority 2: Optimize Parallel Strategy

**Proposal 2A: Reduce Thread Count**
```rust
// Current: Uses all available cores
num_threads: num_cpus::get()

// Proposed: Test with fewer threads
num_threads: (num_cpus::get() / 2).max(1)
```

**Rationale:**
- C++ uses 1 thread by default and achieves 3s
- We use 6 threads and achieve 15s
- Parallelism overhead may exceed benefits

**Expected Impact:**
- Lower memory pressure from thread-local storage
- Simpler coordination
- Potentially faster if overhead currently dominates

---

**Proposal 2B: Batch Processing**
```rust
// Instead of processing contigs one-by-one, batch them
const BATCH_SIZE: usize = 10;

// Process BATCH_SIZE contigs, flush results, repeat
```

**Expected Impact:**
- More predictable memory usage
- Better cache locality
- Easier to reason about peak memory

---

### Priority 3: Compression Tuning

**Proposal 3A: Match C++ Compression Level**
- Current RAGC default: 11
- Current C++ default: 17
- Test with same level to isolate compression vs encoding differences

**Proposal 3B: Investigate Segmentation**
- Compare actual segment boundaries between RAGC and C++
- Ensure splitter logic matches exactly
- Validate in_group_id encoding (already fixed)

---

### Priority 4: Memory Profiling Tools

**Proposal 4A: Add Instrumentation**
```rust
// Add periodic memory reporting
fn log_memory_usage() {
    if let Ok(usage) = sys_info::mem_info() {
        eprintln!("Memory: used={}MB available={}MB",
                  usage.total - usage.avail,
                  usage.avail);
    }
}
```

**Proposal 4B: Integration Tests**
Add memory regression tests:
```rust
#[test]
fn test_memory_usage_stays_under_threshold() {
    let peak_mb = run_compression_and_measure_peak();
    assert!(peak_mb < 500, "Memory usage {peak_mb}MB exceeds 500MB threshold");
}
```

---

## Implementation Roadmap

### Phase 1: Quick Wins (Easy, High Impact)
1. ✅ Test jemalloc (DONE - WORSE: +7% memory, +16% system time)
2. ✅ Reduce channel buffer size from 100 → 10 (DONE - minimal: -0.74% memory)
3. 🔄 Test with fewer threads (1, 2, 4)
4. 🔄 Match C++ compression level

**Results So Far:**
- Jemalloc: Not helpful (allocator not the problem)
- Channel buffers: Not the root cause of 54M page faults

**Updated Theory:**
Page faults likely from:
- Segment compression Vec allocations (hot path)
- Rayon parallel iterator overhead
- Per-segment metadata allocations

**Estimated Time**: 2-4 hours remaining
**Expected Memory Reduction**: 10-20% (revised down)

### Phase 2: Structural Changes (Medium, High Impact)
1. Implement object pooling for Vec<u8> buffers
2. Add batch processing mode
3. Optimize segment metadata structures

**Estimated Time**: 1-2 days
**Expected Memory Reduction**: 40-60%

### Phase 3: Deep Optimization (Hard, Medium Impact)
1. Compare with C++ implementation line-by-line
2. Optimize k-mer data structures (bloom filters?)
3. Investigate zero-copy strategies

**Estimated Time**: 3-5 days
**Expected Memory Reduction**: Additional 10-20%

---

## Success Metrics

### Target Goals
| Metric | Current | Target | Stretch Goal |
|--------|---------|--------|--------------|
| Peak Memory | 984 MB | < 400 MB | < 250 MB |
| Wall Time | 15.2s | < 8s | < 5s |
| System Time | 74.5s | < 10s | < 1s |
| Page Faults | 54.6M | < 5M | < 500K |
| Archive Size | 8.9 MB | < 6.5 MB | < 6.0 MB |

### Verification
After each optimization:
1. Run `/home/erik/ragc/scripts/profile_memory.sh`
2. Compare against baseline
3. Verify correctness with compatibility tests
4. Update `docs/MEMORY_PROFILING.md`

---

## Next Steps

1. Implement Proposal 1A (channel buffer reduction)
2. Implement Proposal 2A (thread count testing)
3. Measure results
4. Decide on Proposals 1B/2B based on impact
5. Iterate

**Start with**: Channel buffer reduction (1 line change, easy to test)
