# RAGC Development Activity Log

## 🔥 CORE DEVELOPMENT PRINCIPLES 🔥

**NO SIMPLIFICATIONS. NO COMPROMISES.**

- **TRUST in your future self** - Conversation compaction will keep you on track
- **TRUST in long refactorings** - Complex changes are fine, do them right
- **DO THE EXACT THING REQUESTED** - If user asks for Option 3, deliver Option 3
- **NEVER suggest "simpler alternatives"** - The user chose the approach for a reason
- **NEVER say "this is complex, let's do something easier"** - That's not your call

When user requests a specific approach (especially after analysis and discussion):
- ✅ DO: Implement exactly what was requested
- ❌ DON'T: Propose "faster" or "simpler" alternatives
- ❌ DON'T: Second-guess the user's technical decisions
- ❌ DON'T: Worry about refactoring scope

**This is a research codebase. Precision and correctness matter more than speed of implementation.**

---

## 🚨 CRITICAL LESSON: CORRECTNESS BEFORE PERFORMANCE 🚨

**Date**: 2025-10-30

### The Lesson

When the user insists on checking correctness before performance optimization - **LISTEN**.

### What Happened

1. **The Setup**: Achieved "100% C++ AGC compatibility" with adaptive mode
   - Archives could be read by C++ AGC ✓
   - Sample lists matched ✓
   - First sample extracted correctly ✓
   - Archive size was 8% larger (seemed minor)

2. **The Mistake**: Wanted to start performance optimizations immediately
   - Found 11x slowdown, wanted to fix parallelization
   - Dismissed 8% size difference as "acceptable for reimplementation"
   - Ready to move on to performance tuning

3. **The User's Intuition**: "I'm very worried that the archive size difference represents a hidden bug"
   - Insisted on systematic correctness verification
   - Wanted to check EVERY part of the algorithm
   - Refused to proceed with performance work

4. **The Discovery**: **CRITICAL DATA CORRUPTION BUG**
   - AAA#0 (first sample): Perfect, byte-for-byte identical ✓
   - AAB#0+ (all other samples): **31% of genome data missing** ✗
   - Archives were readable but contained WRONG DATA
   - The 8% size difference WAS the bug signature

### The Truth

**These systems must produce 100% IDENTICAL output.**
- Like a zlib reimplementation - not "close enough"
- Not just "compatible format"
- Not just "reads correctly"
- **BYTE-FOR-BYTE IDENTICAL compression and decompression**

### Protocol Going Forward

When implementing compression/serialization algorithms:

1. ✅ **Correctness First**: Verify EVERY sample decompresses identically
2. ✅ **Trust Size Differences**: Any size difference is a bug until proven otherwise
3. ✅ **Systematic Verification**: Check all outputs, not just the first one
4. ✅ **Listen to User Intuition**: User's experience trumps passing tests
5. ❌ **Never Skip Verification**: "It reads correctly" ≠ "It's correct"
6. ❌ **Never Rush to Performance**: Correctness cannot be optimized in later

### Remember

**"you wanted to start performance tweaks. no. we weren't done. remember this lesson in claude.md"**

The user was right. Performance means nothing if the output is wrong.

---

## Current Mission: Fix Critical Data Corruption Bug

**Status**: CORRECTNESS BUG FOUND - All samples after first are corrupted (31% data missing)

**Priority**: Fix correctness before ANY performance work

---

## Previous Mission: Inline Split Implementation (Matching C++ AGC Exactly)

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

**Latest Results**:

✅ **Task 1 Complete: Thread Count Testing**
- 1 thread: 14.33s, 997 MB (measured with test_thread_counts.sh)
- 6 threads: 15.08s, 1028 MB (5% slower, 3% more memory)
- **Finding**: Parallelism overhead dominates benefit for this workload

⚠️ **Task 2 Investigation: "Channel 2" Analysis**
- Discovered: CLI uses Rayon-based methods (add_fasta_files_with_splitters)
- These methods ALREADY compress immediately - no "Channel 2" exists in active code
- The old add_contigs_with_splitters (with 3-channel pipeline) is UNUSED by CLI
- **Real bottleneck found**: Line 820 `.collect()` buffers ALL CompressedPacks before writing
- Testing shows actual performance: ~20s wall time, ~1140 MB memory
- **Conclusion**: UncompressedPack optimization had minimal impact (~2% faster)

✅ **Task 3 Complete: Fixed Rayon Threading**
- **ROOT CAUSE FOUND**: Rayon's par_iter() was ignoring num_threads config!
- Default par_iter() uses ALL CPU cores regardless of config setting
- **Fix**: Added ThreadPoolBuilder.num_threads() to all par_iter() call sites
- **Results with fix**:
  - 1 thread: 485 MB, 63.7s (-57.5% memory vs 1140 MB baseline!)
  - 6 threads: 731 MB, 22.1s (-36% memory, comparable speed)
- **Memory gap**: Now 731 MB vs C++ AGC 205 MB = 3.6x (was 4.8x)
- **Key insight**: Proper thread control essential for memory management

**Current Status (After Task 3)**:
- Memory: 731 MB (vs C++ AGC 205 MB = **3.6x gap**)
- Performance: 22s (vs C++ AGC 3s = **7.3x gap**)

✅ **Root Cause Analysis Complete** (see `docs/ACTUAL_ARCHITECTURE.md`)

**RAGC vs C++ AGC - Core Architectural Difference**:

RAGC uses **Batch Mode**:
1. Load ALL segments (Phase 1): **~216 MB**
2. Group segments by k-mer keys
3. Process with Rayon par_iter()
4. Collect ALL packs: **~200 MB**
5. Write sequentially

C++ AGC uses **Streaming Mode**:
1. Stream contigs one-at-a-time into priority queue
2. Workers pull, process, **write immediately**
3. Never holds everything in memory
4. Queue limited by memory size (2GB or 192MB/thread)

**Memory Gap Breakdown**:
```
ALL segments in Phase 1:        +216 MB
Collect all packs in Phase 3:   +200 MB
Rayon threading overhead:        +80 MB
HashMap grouping overhead:       +50 MB
                                --------
Total extra:                    +546 MB
C++ AGC baseline:                205 MB
                                --------
RAGC total:                      751 MB (measured: 731 MB ✓)
```

**The 3.6x gap is architectural** - batch vs streaming.

**Path Forward** (willing to do major refactor):

1. **Option A: Incremental** (~300-400 MB target)
   - Eliminate .collect() in Phase 3: -200 MB
   - Stream Phase 1 (complex): -216 MB
   - Total: ~315 MB (vs 205 MB C++ = 1.5x)

2. **Option B: Full Redesign** (~235 MB target)
   - Replace 3-phase with C++ AGC-style streaming
   - Priority queue + immediate writes
   - Per-thread context reuse
   - Total: ~235 MB (vs 205 MB C++ = 1.1x)

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
