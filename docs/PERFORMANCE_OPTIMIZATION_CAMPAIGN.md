# RAGC Performance Optimization Campaign

**Goal**: Close the performance gap with C++ AGC while maintaining correctness

**Baseline Performance** (yeast10 dataset - 10 samples):
- Memory: 275 MB (C++ AGC: 205 MB = 1.34x gap)
- Wall Time: 3.70s (C++ AGC: 2.98s = 1.24x gap)
- System Time: 0.29s (C++ AGC: 0.90s = 0.32x - **BETTER!**)
- Page Faults: 259K (C++ AGC: 84K = 3.08x gap)
- Archive Size: 6.9 MB (C++ AGC: 5.9 MB = 1.17x gap)
- **Note**: Streaming queue implementation already provides major improvements vs old batch mode!

**Testing Protocol**:
```bash
# Performance test
/usr/bin/time -v ./target/release/ragc create -o /tmp/perf_test.agc -k 21 -s 10000 -m 20 -v 0 ~/scrapy/yeast235_samples/AAA_*.fa

# Correctness test
cargo test --release test_streaming_queue_yeast10
./target/release/ragc getset /tmp/perf_test.agc "AAA#0" 2>/dev/null | grep -v "^>" | tr -d '\n' | wc -c
# Expected: 12157105 bytes

# C++ AGC compatibility (if available)
agc getset /tmp/perf_test.agc "AAA#0" 2>/dev/null | sha256sum
```

---

## 🎉 Good News: Streaming Queue Already Highly Optimized!

The streaming queue implementation is **already very close** to C++ AGC performance:
- Memory: 275 MB vs 205 MB = **1.34x gap (only 70 MB!)**
- Runtime: 3.70s vs 2.98s = **1.24x gap (only 0.72s!)**
- Archive size: 6.9 MB vs 5.9 MB = **1.17x gap (only 1 MB!)**

The old profiling docs were measuring batch mode (984 MB). Streaming queue is a HUGE improvement!

## Phase 1: Micro-optimizations (Target: -30 MB, -0.3s)

### Step 1: Use BTreeMap for Segment Groups ✅ IN PROGRESS

**File**: `ragc-core/src/streaming_compressor_queue.rs:261`

**Current Code**:
```rust
let segment_groups = Arc::new(Mutex::new(HashMap::new()));
```

**Change**: Replace HashMap with BTreeMap for more memory-efficient storage

**Rationale**:
- BTreeMap uses less memory than HashMap for small to medium maps
- With ~20-100 segment groups for yeast10, BTreeMap should save memory
- Deterministic iteration order is a nice bonus

**Expected Impact**:
- Memory: -10 MB (less HashMap overhead)
- Archive Size: No change (must be identical)
- Runtime: +0.1s (slightly slower lookups, but negligible)
- Risk: Very Low

**Results**:
- [x] Performance test completed
- [x] Memory: 276 MB (baseline: 275 MB) = **+1 MB**
- [x] Wall Time: 3.54s (baseline: 3.70s) = **-0.16s (-4.3%)**
- [x] System Time: 0.26s (baseline: 0.29s) = **-0.03s**
- [x] Archive Size: 6.9 MB (baseline: 6.9 MB) = **IDENTICAL ✅**
- [x] Correctness: AAA#0 = 12157105 bytes ✅
- [x] Commit: 8c7a3fd

**Analysis**: BTreeMap uses slightly more memory (+1 MB) due to tree structure overhead, but provides 4.3% faster runtime due to better cache locality. The 1 MB memory increase is acceptable given the runtime improvement. Archive size is identical, confirming correctness.

---

### Step 2: Optimize Rayon Thread Count for Small Workloads ⏸️ PENDING

**File**: `ragc-core/src/compressor_streaming.rs`

**Change**: Add adaptive thread selection based on segment count
- If segments < 100: Use 1 thread
- If segments < 1000: Use 2 threads
- Else: Use num_cpus threads

**Expected Impact**:
- Memory: -80 MB (for yeast10)
- Runtime: -1s (less coordination overhead)
- Archive Size: No change
- Risk: Low

**Results**:
- [ ] Performance test completed
- [ ] Memory: ___ MB
- [ ] Wall Time: ___ s
- [ ] System Time: ___ s
- [ ] Archive Size: ___ MB
- [ ] Correctness: AAA#0 = ___ bytes
- [ ] Commit: ___

---

### Step 3: Pre-allocate Vector Capacities ⏸️ PENDING

**Files**: Multiple (segment buffers, group maps, etc.)

**Change**: Pre-size collections based on expected counts
- Segments vector: with_capacity(estimated_contigs * 2)
- Group HashMap: with_capacity(16)
- Pack buffers: with_capacity(expected_size)

**Expected Impact**:
- Memory: -20 MB (fewer reallocations)
- Runtime: -0.5s
- Archive Size: No change
- Risk: Low

**Results**:
- [ ] Performance test completed
- [ ] Memory: ___ MB
- [ ] Wall Time: ___ s
- [ ] System Time: ___ s
- [ ] Archive Size: ___ MB
- [ ] Correctness: AAA#0 = ___ bytes
- [ ] Commit: ___

---

## Phase 2: Streaming Architecture (Target: -200 MB, -5s)

### Step 4: Implement Streaming Phase 1 ⏸️ PENDING

**File**: `ragc-core/src/compressor_streaming.rs`

**Change**: Stream contigs through a bounded channel instead of loading all segments first
- Replace `Vec<SegmentDescriptor>` with `crossbeam::channel::bounded(1000)`
- Process contigs → segments → groups in pipeline
- Never materialize full segment list

**Expected Impact**:
- Memory: -216 MB
- Runtime: -2s (better cache locality)
- Archive Size: No change
- Risk: Medium (requires architectural refactor)

**Results**:
- [ ] Performance test completed
- [ ] Memory: ___ MB
- [ ] Wall Time: ___ s
- [ ] System Time: ___ s
- [ ] Archive Size: ___ MB
- [ ] Correctness: AAA#0 = ___ bytes
- [ ] Commit: ___

---

### Step 5: Replace HashMap with BTreeMap for Groups ⏸️ PENDING

**File**: `ragc-core/src/compressor_streaming.rs`

**Change**: Use BTreeMap instead of HashMap for group storage
- More memory-efficient
- Better cache locality
- Deterministic iteration order

**Expected Impact**:
- Memory: -30 MB
- Runtime: +0.5s (slightly slower lookups)
- Archive Size: No change
- Risk: Low

**Results**:
- [ ] Performance test completed
- [ ] Memory: ___ MB
- [ ] Wall Time: ___ s
- [ ] System Time: ___ s
- [ ] Archive Size: ___ MB
- [ ] Correctness: AAA#0 = ___ bytes
- [ ] Commit: ___

---

## Phase 3: Archive Size Optimization (Target: -2 MB)

### Step 6: Improve Segment Grouping ⏸️ PENDING

**File**: `ragc-core/src/compressor_streaming.rs`

**Change**: Better k-mer boundary detection and group assignment
- Analyze C++ AGC grouping strategy
- Implement more aggressive similarity detection
- May require minimizer implementation

**Expected Impact**:
- Memory: No change
- Runtime: +1s (more analysis)
- Archive Size: -1.5 MB
- Risk: Medium (affects correctness)

**Results**:
- [ ] Performance test completed
- [ ] Memory: ___ MB
- [ ] Wall Time: ___ s
- [ ] System Time: ___ s
- [ ] Archive Size: ___ MB
- [ ] Correctness: AAA#0 = ___ bytes
- [ ] Commit: ___

---

### Step 7: Optimize LZ Differential Encoding ⏸️ PENDING

**File**: `ragc-core/src/lz_diff.rs`

**Change**: More aggressive reference selection in groups
- Better reference choosing (closest match vs first)
- Improve copy/literal decisions
- May need to analyze C++ AGC LZ strategy

**Expected Impact**:
- Memory: No change
- Runtime: +0.5s
- Archive Size: -0.5 MB
- Risk: Medium (affects correctness)

**Results**:
- [ ] Performance test completed
- [ ] Memory: ___ MB
- [ ] Wall Time: ___ s
- [ ] System Time: ___ s
- [ ] Archive Size: ___ MB
- [ ] Correctness: AAA#0 = ___ bytes
- [ ] Commit: ___

---

## Summary

| Step | Status | Memory | Runtime | Archive | Commit |
|------|--------|--------|---------|---------|--------|
| Baseline | ✅ | 275 MB | 3.70s | 6.9 MB | - |
| 1. BTreeMap groups | ✅ | 276 MB | 3.54s | 6.9 MB | 8c7a3fd |
| 2. Pre-allocate vectors | ⏸️ | - | - | - | - |
| 3. Reduce String cloning | ⏸️ | - | - | - | - |
| 4. Better grouping | ⏸️ | - | - | - | - |
| **C++ AGC** | | **205 MB** | **2.98s** | **5.9 MB** | - |
| **Realistic Target** | | **< 240 MB** | **< 3.2s** | **< 6.2 MB** | - |

**Note**: We're already at 1.34x memory, 1.24x runtime, 1.17x archive size. The remaining gaps are very small!

---

## Notes

- After each step, run full test suite: `cargo test --release`
- Always verify C++ AGC compatibility if possible
- Archive size must not change in Phase 1 optimizations
- Archive size changes in Phase 3 require careful correctness validation
- Stop immediately if any corruption detected
