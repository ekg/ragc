# RAGC Parallelization Bug Analysis

**Date**: 2025-10-23
**Issue**: RAGC does not parallelize despite creating multiple threads
**Impact**: 5.5x slower than C++ AGC, no speedup from 1 → 6 → 12 threads

---

## Benchmark Results (yeast235: 37 files, 644 MB)

### RAGC Performance

| Threads | Wall Time | User | System | CPU % | Memory | vs 1-thread |
|---------|-----------|------|--------|-------|---------|-------------|
| 1 | 127s | 47.7s | 81.4s | **101%** | 918 MB | 1.0x |
| 6 | 136s | 56.7s | 95.2s | **111%** | 1003 MB | **0.93x** (SLOWER!) |
| 12 | 136s | 56.7s | 95.2s | **111%** | 1041 MB | **0.93x** (SLOWER!) |

**Problem**: Only ~1 CPU core utilized (111%) despite 6 or 12 threads!

### C++ AGC Performance (Perfect Scaling)

| Threads | Wall Time | User | System | CPU % | Memory | Speedup |
|---------|-----------|------|--------|-------|---------|---------|
| 1 | 23s | 23.2s | 0.4s | 101% | 968 MB | 1.0x |
| 6 | 6.2s | 31.5s | 0.5s | **516%** | 947 MB | **3.7x** ✅ |
| 12 | 4.8s | 43.3s | 0.8s | **917%** | 928 MB | **4.8x** ✅ |

---

## Root Cause: Global Lock Contention

### The Bottleneck (ragc-core/src/compressor_streaming.rs:2544-2562)

```rust
// EVERY segment locks the global HashMap!
let pack_opt = {
    let mut groups_lock = groups.lock()  // ← ALL THREADS BLOCK HERE!
        .map_err(|e| anyhow::anyhow!("Failed to lock groups: {}", e))?;

    let (_group_id, group_writer) = groups_lock
        .entry(seg_info.key.clone())
        .or_insert_with(|| {
            let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
            (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
        });

    // Add segment - still holding global lock!
    group_writer.add_segment(seg_info.segment, &config)?
};  // Lock released
```

### Why This Kills Parallelism

1. **Global exclusive lock** on `Arc<Mutex<HashMap>>` for every segment
2. With 37 files × thousands of segments = **tens of thousands of lock acquisitions**
3. **Only 1 thread can process a segment at a time** - completely serialized
4. Other threads block waiting for the lock
5. Result: **111% CPU** (only ~1 core active)

---

## C++ AGC's Solution

### Fine-Grained Locking Architecture

From `/home/erik/agc/src/core/agc_compressor.h`:

```cpp
// Line 608: Shared mutex allows multiple readers!
shared_mutex seg_map_mtx;

// Line 628: Map protected by SHARED mutex
unordered_map<pair<uint64_t, uint64_t>, int32_t, MurMurPair64Hash> map_segments;
// Comment: "shared_mutex (seg_map_mtx)"

// Line 630: Each segment has internal mutex
vector<shared_ptr<CSegment>> v_segments;
// Comment: "shared_mutex to vector (seg_vec_mtx) + internal mutexes in stored objects"
```

### Key Differences

| Operation | C++ AGC | RAGC |
|-----------|---------|------|
| **Find group** | Shared lock (readers) | Exclusive lock |
| **Create group** | Exclusive lock (brief) | Exclusive lock |
| **Add to group** | Per-group mutex | Still global lock! |
| **Parallelism** | High (multiple readers + per-group locks) | None (global exclusive) |

**C++ AGC's approach**:
1. **shared_mutex** for `map_segments` allows concurrent lookups
2. **Per-segment internal mutexes** for adding data to groups
3. Only group *creation* requires brief exclusive lock
4. Result: **Nearly perfect scaling** (4.8x with 12 threads)

---

## Performance Impact

### Speed Comparison

| Benchmark | RAGC | C++ AGC | Ratio |
|-----------|------|---------|-------|
| **1 thread** | 127s | 23s | **5.5x slower** |
| **6 threads** | 136s | 6.2s | **22x slower** |
| **12 threads** | 136s | 4.8s | **28x slower** |

### Scaling Efficiency

| Threads | RAGC Speedup | C++ AGC Speedup | RAGC Efficiency |
|---------|--------------|-----------------|-----------------|
| 1 | 1.0x | 1.0x | 100% |
| 6 | **0.93x** | 3.7x | **-7%** (NEGATIVE!) |
| 12 | **0.93x** | 4.8x | **-7%** (NEGATIVE!) |

**RAGC gets SLOWER with more threads** due to lock contention overhead!

---

## Solution Options

### Option A: RwLock (Quick Fix)

Replace `Arc<Mutex<HashMap>>` with `Arc<RwLock<HashMap>>`:

```rust
let groups = Arc::new(RwLock::new(HashMap::new()));

// Workers:
let pack_opt = {
    let groups_read = groups.read()?;  // Shared lock for lookup

    if let Some(group) = groups_read.get(&seg_info.key) {
        // Add to existing group - but still need to lock HashMap for writes!
        // This doesn't fully solve it...
    } else {
        drop(groups_read);
        let mut groups_write = groups.write()?;  // Exclusive for insert
        // ...
    }
};
```

**Problem**: Still requires write lock for `group.add_segment()` mutations.

### Option B: DashMap (Concurrent HashMap) ✅ RECOMMENDED

Use `dashmap::DashMap` - a concurrent HashMap with per-shard locking:

```rust
use dashmap::DashMap;

// Concurrent HashMap - no Arc<Mutex> needed!
let groups = DashMap::<SegmentGroupKey, (u32, GroupWriter)>::new();

// Workers can access concurrently:
let pack_opt = {
    let mut group_entry = groups.entry(seg_info.key.clone())
        .or_insert_with(|| {
            let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
            (gid, GroupWriter::new(...))
        });

    // Add segment - only locks this specific entry!
    group_entry.1.add_segment(seg_info.segment, &config)?
};
```

**Advantages**:
- Per-key locking (like C++ AGC's per-segment mutexes)
- Multiple threads can access different groups concurrently
- Minimal API changes
- Battle-tested crate (used in production)

### Option C: Parking Lot RwLock (Middle Ground)

Use `parking_lot::RwLock` which is faster than std::sync::RwLock:

```rust
use parking_lot::RwLock;

let groups = Arc::new(RwLock::new(HashMap::new()));
```

**Advantages**:
- Faster than std RwLock
- Still allows concurrent reads

**Disadvantages**:
- Still requires exclusive lock for add_segment mutations
- Not as fine-grained as DashMap

---

## Recommended Solution: DashMap

### Implementation Steps

1. **Add DashMap dependency** to `ragc-core/Cargo.toml`:
   ```toml
   dashmap = "6.0"
   ```

2. **Change groups type** (line 2364):
   ```rust
   // OLD:
   let groups = Arc::new(Mutex::new(HashMap::new()));

   // NEW:
   let groups = Arc::new(DashMap::<SegmentGroupKey, (u32, GroupWriter)>::new());
   ```

3. **Update group access** (lines 2543-2562):
   ```rust
   // NEW: DashMap entry API
   let pack_opt = {
       let mut entry = groups.entry(seg_info.key.clone())
           .or_insert_with(|| {
               let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
               (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
           });

       // Only this entry is locked, not entire HashMap!
       entry.1.add_segment(seg_info.segment, &config)?
   };
   ```

4. **Update flush** (line 2642-2686):
   ```rust
   // NEW: DashMap iteration
   for mut entry in groups.iter_mut() {
       let (group_id, group_writer) = entry.pair_mut();
       // ...
   }
   ```

### Expected Results

| Metric | Current | Expected After Fix | Improvement |
|--------|---------|-------------------|-------------|
| **6-thread speedup** | 0.93x | 3.5x | **~4x faster** |
| **12-thread speedup** | 0.93x | 4.5x | **~5x faster** |
| **Wall time (6T)** | 136s | ~36s | **-100s** |
| **CPU utilization** | 111% | 600%+ | **5-6 cores active** |
| **vs C++ AGC** | 22x slower | ~6x slower | Still slower but acceptable |

---

## Why Still Slower Than C++ AGC?

Even with DashMap, RAGC will likely remain slower than C++ AGC due to:

1. **High system time** (81-95s vs 0.4s)
   - 120M page faults vs 5M
   - Memory allocation churn (49M Vec growth calls in LZDiff)
   - Rust allocator overhead

2. **Compression performance**
   - ZSTD binding overhead
   - Buffer management differences

3. **I/O patterns**
   - Different FASTA reading approaches
   - Buffer sizes

### Future Optimizations (After DashMap Fix)

1. **Reduce system time**:
   - Pre-allocate Vec capacities in LZDiff
   - Try jemalloc or mimalloc
   - Reduce allocation churn

2. **Optimize splitter phase**:
   - 153 MB peak during HashMap growth
   - Pre-size HashMap
   - Stream splitters differently

3. **Buffer pooling**:
   - Reuse ZSTD contexts per thread
   - Pool segment buffers

---

## Files to Modify

- `ragc-core/Cargo.toml` - Add dashmap dependency
- `ragc-core/src/compressor_streaming.rs`:
  - Line 2364: Change groups type
  - Lines 2543-2562: Update group access
  - Lines 2642-2686: Update flush
  - Lines 2495-2611: Review worker loop

---

## Testing Plan

1. **Functional test**: Verify output matches current implementation
   ```bash
   cargo test
   ```

2. **Benchmark (yeast235)**:
   ```bash
   /usr/bin/time -v ./target/release/ragc create -o test.agc -k 21 -s 10000 -m 20 -t 6 /home/erik/scrapy/yeast235_split/*.fa
   ```

   Expect:
   - CPU: 500-600% (was 111%)
   - Wall: ~35-40s (was 136s)

3. **Compatibility test**:
   ```bash
   agc listset test.agc
   agc getset test.agc sample1
   ```

4. **Scaling test** (1, 6, 12 threads):
   - Confirm near-linear speedup
   - Target: 3.5x with 6 threads, 4.5x with 12 threads

---

## Conclusion

The parallelization bug is caused by using `Arc<Mutex<HashMap>>` which creates a global bottleneck. **Every segment** acquisition requires an exclusive lock, completely serializing the workload.

C++ AGC uses:
- `shared_mutex` for concurrent reads
- Per-segment internal mutexes for fine-grained locking

The solution is to use `DashMap` which provides:
- Concurrent HashMap with per-shard locking
- Similar performance to C++ AGC's approach
- Minimal code changes

**Expected improvement**: 3.5-4.5x speedup with 6-12 threads, making RAGC competitive with C++ AGC's parallelization.

---

## References

- C++ AGC source: `/home/erik/agc/src/core/agc_compressor.{h,cpp}`
- RAGC source: `/home/erik/ragc/ragc-core/src/compressor_streaming.rs`
- DashMap crate: https://docs.rs/dashmap
- Benchmark logs: `/tmp/yeast235_bench/`
