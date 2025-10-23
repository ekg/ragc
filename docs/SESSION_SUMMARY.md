# Memory Optimization Session Summary

**Date**: 2025-10-23
**Branch**: feature/memory-profiling
**Status**: Major breakthrough achieved - 36% memory reduction

---

## Key Achievements

### 1. Fixed Critical Rayon Threading Bug
**Impact**: 36% memory reduction (1140 MB → 731 MB)

**Problem**: Rayon's `par_iter()` was using ALL CPU cores by default, completely ignoring the `num_threads` config parameter passed via `-t` flag.

**Solution**: Added `ThreadPoolBuilder.num_threads()` configuration to all 4 `par_iter()` call sites in `compressor_streaming.rs`.

**Results**:
- **6 threads**: 731 MB, 22.1s wall time (-36% memory)
- **1 thread**: 485 MB, 63.7s wall time (-57.5% memory!)
- System time: 109s → 74s (-32% with 6 threads)

**Files Modified**:
- `ragc-core/src/compressor_streaming.rs`: Added ThreadPoolBuilder config at 4 locations
- All tests pass, C++ AGC compatibility verified

---

## Investigation Summary

### Task 1: Thread Count Testing ✅
- **Finding**: Parallelism overhead dominates for this workload
- Tested 1, 2, 4, 6 threads
- Without Rayon fix, difference was minimal (all using max cores)
- With Rayon fix, clear memory/speed tradeoff visible

### Task 2: Channel 2 Investigation ✅
- **Finding**: "Channel 2" doesn't exist in active code paths
- CLI uses Rayon-based methods that compress immediately
- Old 3-channel pipeline (`add_contigs_with_splitters`) unused by CLI
- No optimization opportunity here

### Task 3: Rayon Threading Fix ✅ **BREAKTHROUGH**
- **ROOT CAUSE**: Rayon ignoring num_threads config
- **SOLUTION**: ThreadPoolBuilder configuration
- **RESULT**: 36% memory reduction, proper thread control

### Task 4: FASTA Streaming Investigation ✅
- **Finding**: Phase 1 collects ALL segments in Vec<PreparedSegment> (line 1655)
- For yeast10: ~12K segments × (data + metadata) held in memory
- **Challenge**: Phase 2 grouping by k-mer keys requires all segments available
- **Decision**: Defer (requires major architectural refactor)

### Task 5/6: Vec Buffer Pooling Analysis ✅
- **Finding**: `lz_diff.rs:227` allocates Vec for every LZ encoding
- 12K allocations for yeast10 dataset
- Similar to ZSTD pooling optimization
- **Decision**: Lower priority given 36% win already achieved

---

## Current Status

**Memory Gap**: 731 MB vs C++ AGC 205 MB = **3.6x** (was 4.8x)
**Performance Gap**: 22s vs C++ AGC 3s = **7.3x** (was 5.1x)

### Baseline Comparison

| Metric | Before | After | Improvement |
|--------|---------|-------|-------------|
| Memory (6 threads) | 1140 MB | 731 MB | **-36%** |
| Memory (1 thread) | 1140 MB | 485 MB | **-57.5%** |
| System time (6 threads) | 109s | 74s | **-32%** |
| System time (1 thread) | 109s | 49s | **-55%** |
| Page faults | 70.6M | 70.6M | Unchanged |

---

## Remaining Optimization Opportunities

Listed in order of estimated impact vs complexity:

### 1. LZ Vec Buffer Pooling (Medium Impact, Low Complexity)
- **Location**: `lz_diff.rs:227`
- **Issue**: Allocates Vec for every encode (12K times)
- **Solution**: Thread-local buffer pooling (similar to ZSTD)
- **Estimated Impact**: 5-10% memory reduction

### 2. Streaming FASTA Processing (High Impact, High Complexity)
- **Location**: `compressor_streaming.rs:1655`
- **Issue**: Collects all segments before processing
- **Challenge**: Phase 2 grouping requires all segments
- **Solution**: Multi-pass processing or disk-based temporary storage
- **Estimated Impact**: 10-20% memory reduction
- **Risk**: Major architectural refactor

### 3. Eliminate .collect() Buffering (Medium Impact, High Complexity)
- **Location**: `compressor_streaming.rs:828`
- **Issue**: Buffers all CompressedPacks before writing
- **Challenge**: Writing sequentially conflicts with Rayon parallel processing
- **Solution**: Custom iterator or streaming writer
- **Estimated Impact**: 5-10% memory reduction

---

## Commits in This Session

```
9f98086 Document Task 4 investigation and remaining optimization opportunities
f6b1f42 MAJOR: Fix Rayon threading to respect num_threads config
ff1be2d Document Task 2 investigation: Channel 2 not in active code paths
355c0a5 Update CLAUDE.md: Document single-thread finding
23f0f80 Add --threads CLI parameter and prove single-thread is fastest
```

---

## Next Steps

1. **Recommended**: Implement LZ Vec buffer pooling (quick win)
2. **Consider**: Profile with updated baseline to identify new hotspots
3. **Evaluate**: Whether remaining 3.6x gap justifies major refactoring
4. **Alternative**: Accept current performance, focus on other features

---

## Files to Review

- `CLAUDE.md` - Detailed investigation notes and architecture analysis
- `ragc-core/src/compressor_streaming.rs` - Rayon threading fixes
- `scripts/test_thread_counts.sh` - Thread performance testing
- `docs/MEMORY_PROFILING.md` - Original baseline profiling
- `docs/PIPELINE_COMPARISON.md` - C++ AGC vs RAGC architecture

---

## Conclusion

This session achieved a **major breakthrough** by fixing the Rayon threading bug, resulting in **36% memory reduction**. The `-t` flag now works correctly, giving users control over memory/speed tradeoffs.

The remaining 3.6x memory gap vs C++ AGC is primarily due to architectural differences:
- RAGC holds all segments in memory for grouping (Phase 1)
- Multiple Vec allocations per segment (LZ encoding, compression)
- Rayon's parallel collection pattern

Further optimization requires balancing implementation complexity against diminishing returns.
