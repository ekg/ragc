# RAGC Development Activity Log

## Current Mission: Memory Optimization (Feature Branch: memory-profiling)

**Goal**: Reduce RAGC memory usage from 984 MB to ~300 MB (matching C++ AGC performance)

**Status**: Phase 1 complete (profiling + analysis), Phase 2 in progress (implementation)

---

## Progress Summary

### ✅ Completed

1. **Baseline Profiling** (PR #2)
   - Full C++ AGC compatibility achieved
   - Memory profiling infrastructure established
   - Identified 4.81x memory gap (RAGC: 984 MB vs C++ AGC: 205 MB)
   - Identified 5.10x performance gap (RAGC: 15.2s vs C++ AGC: 3.0s)
   - **Critical**: 54.6M page faults vs 84K (650x more)

2. **Root Cause Analysis**
   - ❌ NOT allocator issue (jemalloc made it worse)
   - ❌ NOT channel buffer size (minimal impact)
   - ✅ **FOUND**: 3-stage pipeline creates excessive Vec allocations
   - ✅ **FOUND**: `zstd::encode_all()` was creating contexts per-segment
   - ✅ **FOUND**: C++ AGC uses priority queue + context reuse, RAGC uses channels

3. **C++ AGC Architecture Analysis**
   - Documented in `docs/PIPELINE_COMPARISON.md`
   - C++ uses single priority queue with memory-based capacity (2GB or 192MB/thread)
   - C++ reuses ZSTD contexts per-thread
   - C++ writes directly: segment → compress → write (no intermediate buffering)
   - C++ uses 1 thread efficiently vs RAGC's 6 threads inefficiently

4. **Initial Optimizations**
   - ZSTD buffer pooling: -20 MB (-2% memory)
   - Channel buffer reduction: -0.74% memory (minimal)
   - Learned: 54M page faults from pipeline Vec allocations, not ZSTD contexts

### 🔄 Current Phase: Systematic Pipeline Optimization

**Latest Results** (Task 1 Complete):

✅ **Single-thread is FASTEST!**
- 1 thread: 14.33s, 997 MB (BEST)
- 6 threads: 15.08s, 1028 MB (5% slower, 3% more memory)
- **Proof**: Parallelism overhead > parallelism benefit
- **BUT**: Still 4.8x slower than C++ AGC (14.33s vs 3.0s)
- **Conclusion**: Pipeline architecture is the bottleneck, not threads

**Next Steps** (in priority order):

1. ✅ **Test single-thread performance** - PROVEN: 1 thread is optimal
2. 🔄 **Remove intermediate buffering stage** - Eliminate Channel 2, compress immediately after LZ
3. **Stream FASTA reading** - Don't load all contigs into memory
4. **Implement Vec buffer pooling** - Reuse buffers for LZ/compression output
5. **Replace 3-channel pipeline** - Move to C++ AGC-style priority queue (if needed)

---

## Key Files & Documentation

### Documentation
- `docs/MEMORY_PROFILING.md` - Baseline profiling results and analysis
- `docs/OPTIMIZATION_PROPOSALS.md` - Comprehensive optimization roadmap
- `docs/PIPELINE_COMPARISON.md` - **CRITICAL**: C++ AGC vs RAGC architecture comparison

### Profiling Scripts
- `scripts/profile_memory.sh` - RAGC vs C++ AGC comparison
- `scripts/profile_allocators.sh` - Allocator testing
- `scripts/test_buffer_optimization.sh` - Channel buffer testing
- `scripts/test_zstd_pooling.sh` - ZSTD pooling verification

### Source Files (Modified)
- `ragc-core/src/zstd_pool.rs` - Thread-local ZSTD buffer pooling
- `ragc-core/src/segment_compression.rs` - Uses pooled compression
- `ragc-core/src/compressor_streaming.rs` - Main compression pipeline (NEEDS REFACTOR)

---

## Performance Targets

| Metric | Current | Target | Stretch |
|--------|---------|--------|---------|
| **Peak Memory** | 984 MB | < 400 MB | < 250 MB |
| **Wall Time** | 15.2s | < 8s | < 5s |
| **System Time** | 74.5s | < 10s | < 1s |
| **Page Faults** | 54.6M | < 5M | < 500K |
| **Archive Size** | 8.9 MB | < 6.5 MB | < 6.0 MB |

**Success Criteria**:
- Memory within 2x of C++ AGC (< 400 MB)
- Performance within 2x of C++ AGC (< 6s)
- All C++ compatibility tests passing

---

## Testing Strategy

After each optimization:
1. Run `scripts/profile_memory.sh` to measure impact
2. Verify correctness: `cargo test --release`
3. Verify C++ compatibility: CI tests
4. Update `docs/MEMORY_PROFILING.md` with results
5. Commit with clear performance metrics

---

## Architecture Notes

### Current RAGC Pipeline (INEFFICIENT)
```
FASTA → Load ALL → Channel 1 → Workers → LZ encode → Vec<u8>
                       ↓
                   Channel 2 → Compress → Vec<u8>
                       ↓
                   Channel 3 → Write
```

**Problem**: 3 buffering stages × Vec allocations per segment = 36K+ allocations for yeast10

### C++ AGC Pipeline (EFFICIENT)
```
FASTA → Stream → Priority Queue (memory-limited) → Worker Threads
                                                         ↓
                                                    Per-thread ZSTD_CCtx (reused)
                                                         ↓
                                                    Segment → compress → write
```

**Advantage**: Direct path, memory-based limits, context reuse

---

## Branch Strategy

- **Branch**: `feature/memory-profiling`
- **PR**: #2 (https://github.com/ekg/ragc/pull/2)
- **Base**: main
- **CI Status**: All passing
- **Strategy**: Systematic optimization commits with performance metrics

---

## Debug Tips

### Measuring Memory Impact
```bash
# Quick test on yeast10
/usr/bin/time -v ./target/release/ragc create -o test.agc -k 21 -s 10000 -m 20 -v 0 samples/*.fa

# Look for:
# - Maximum resident set size (kbytes)
# - System time (seconds)
# - Minor (reclaiming a frame) page faults
```

### Profiling Hotspots
```bash
# Build with debug symbols
cargo build --release --profile release-with-debug

# Profile (if needed later)
cargo flamegraph --root -- create -o test.agc ...
```

### Verify Compatibility
```bash
# Create with RAGC
./target/release/ragc create -o ragc.agc ...

# Read with C++ AGC
/home/erik/agc/bin/agc getset ragc.agc sample_name > output.fa

# Verify roundtrip
diff original.fa output.fa
```

---

## Context for Next Session

**Where we are**:
- Profiling complete, root cause identified
- ZSTD buffer pooling implemented (-2% memory)
- Ready to systematically optimize pipeline

**What's next**:
- Test single-thread vs multi-thread performance
- Remove intermediate buffering stages
- Work towards C++ AGC-style architecture

**Key insight**: RAGC's parallelism is SLOWER than C++ AGC's single thread due to coordination overhead. The 3-stage pipeline creates unnecessary allocations that dwarf any benefits from parallelism.

**Expected timeline**: With focused work, should achieve 2x reduction (to ~400 MB) within a few optimization iterations. Full C++ AGC parity (~250 MB) may require more extensive refactoring.
