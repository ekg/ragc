# RAGC Optimization Plan - Starting from Working Code (4a7f895)

**Date**: 2025-11-02
**Current Status**: Working code with 100% correctness, but poor parallelization
**Goal**: Improve from 102% CPU (essentially single-threaded) to proper multi-threaded performance

## Current Performance (4a7f895)

```
Yeast10 (10 samples, 114M bases):
- Time: ~15-20s
- CPU: 102% (with 15 threads configured!)
- Memory: Unknown (need to measure)
- Archive: 7.6MB
- Correctness: ✅ 100% (0% data loss)
```

**Problem**: Sequential processing in main thread, only writer thread runs in parallel (doing almost nothing)

## Current Architecture Analysis

```
┌─────────────────────────────────────────────────────────┐
│ Main Thread (100% of CPU time)                         │
├─────────────────────────────────────────────────────────┤
│ for each contig:                                        │
│   1. Read from file (sequential)                        │
│   2. Segment at splitters (CPU-intensive, sequential)   │
│   3. For each segment:                                  │
│      - Check for inline splitting (sequential)          │
│      - Compress with ZSTD (CPU-intensive, sequential)   │
│      - Send to writer thread                            │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│ Writer Thread (2% of CPU time)                          │
├─────────────────────────────────────────────────────────┤
│ while packs available:                                  │
│   - Receive pack from channel                           │
│   - Write to archive (I/O, mostly waiting)              │
│   - Register metadata                                   │
└─────────────────────────────────────────────────────────┘
```

**Bottleneck**: Main thread does ALL segmentation and compression sequentially!

## Optimization Strategy

Based on C++ AGC study, we need to parallelize the CPU-intensive work:

### Phase 1: Quick Win - Parallelize Compression Only

Keep segmentation sequential but parallelize compression:

```rust
// Sequential: Read + Segment
let all_segments: Vec<UncompressedSegment> = Vec::new();
while let Some(contig) = iterator.next_contig()? {
    let segments = split_at_splitters(contig);
    for seg in segments {
        all_segments.push(seg);  // Buffer instead of compressing immediately
    }
}

// Parallel: Group + Compress
let packs: Vec<CompressedPack> = all_segments
    .par_iter()  // Rayon parallel iterator
    .group_by(|seg| seg.group_key)
    .map(|(key, segments)| {
        let gid = get_or_create_group(key);
        compress_segments(gid, segments)
    })
    .collect();

// Sequential: Write
for pack in packs {
    archive.write(pack)?;
}
```

**Expected improvement**: 3-5x speedup (compression is ~50-70% of runtime)

**Pros**:
- Small change, low risk
- Keeps file reading sequential (simpler)
- Parallelizes the most CPU-intensive part

**Cons**:
- Still sequential segmentation
- Buffers all segments in memory (but they're small ~10KB each)

### Phase 2: Full Parallelization - Match C++ AGC

Parallelize both segmentation AND compression:

```rust
// 1. Load all contigs into memory first
let all_contigs: Vec<(String, String, Vec<u8>)> = ...;

// 2. Parallel segmentation + buffering
let segments: DashMap<SegmentGroupKey, Vec<UncompressedSegment>> = DashMap::new();

all_contigs.par_iter().for_each(|(sample, contig_name, seq)| {
    let segs = split_at_splitters(seq, splitters);
    for seg in segs {
        segments.entry(seg.key).or_default().push(seg);
    }
});

// 3. Registration (thread-safe, sequential is fine)
let mut group_ids = HashMap::new();
for (key, _) in &segments {
    let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
    group_ids.insert(key, gid);
    register_streams(gid);
}

// 4. Parallel compression
segments.par_iter().for_each(|(key, segs)| {
    let gid = group_ids[key];
    let pack = compress_segments(gid, segs);
    send_to_writer(pack);
});
```

**Expected improvement**: 8-12x speedup (both segmentation and compression parallelized)

**Pros**:
- Maximum parallelization
- Matches C++ AGC architecture
- Should get close to C++ AGC performance

**Cons**:
- Loads all contigs in memory (could be ~500MB for yeast10)
- More complex change

### Phase 3: Memory-Efficient Streaming (Optional)

Match C++ AGC exactly with BoundedPriorityQueue:

```rust
// 1. Worker threads pull from queue
let queue = BoundedPriorityQueue::new(capacity);

for _ in 0..num_threads {
    let queue = queue.clone();
    thread::spawn(move || {
        while let Some(contig) = queue.pop() {
            let segments = split_at_splitters(contig);
            for seg in segments {
                buffer.add(seg);  // Thread-safe buffer
            }
        }
    });
}

// 2. Main thread feeds queue
for contig in contigs {
    queue.push(contig);
}

// 3. Registration + compression (same as Phase 2)
```

**Expected improvement**: Same as Phase 2, but lower memory

**Pros**:
- Memory-efficient (streaming)
- Exact C++ AGC architecture
- Handles huge datasets

**Cons**:
- Most complex implementation
- Overkill for small datasets

## Recommendation: Start with Phase 1

### Implementation Plan

1. **Measure current bottlenecks** (5 min)
   ```bash
   cargo build --release
   time ./target/release/ragc create ...
   # Look at perf/htop to confirm single-thread bottleneck
   ```

2. **Implement Phase 1 optimization** (30-60 min)
   - Modify `add_segments_with_inline_splits` to buffer all segments first
   - Use Rayon to parallelize compression
   - Test correctness with yeast10
   - Measure speedup

3. **If Phase 1 works well, implement Phase 2** (2-3 hours)
   - Load all contigs in memory
   - Parallelize segmentation with Rayon
   - Test and measure

4. **If memory is a problem, implement Phase 3** (1 day)
   - Implement BoundedPriorityQueue-based streaming
   - Match C++ AGC architecture exactly

## Expected Results

### Phase 1: Parallel Compression
```
Yeast10 (10 samples, 114M bases):
- Time: ~5-7s (3-4x faster)
- CPU: 400-600% (4-6 threads utilized)
- Memory: ~100-200MB (buffering segments)
- Archive: 7.6MB (unchanged)
- Correctness: Should be ✅ 100%
```

### Phase 2: Full Parallelization
```
Yeast10 (10 samples, 114M bases):
- Time: ~2-3s (5-8x faster, close to C++ AGC's 3.0s)
- CPU: 800-1200% (8-12 threads utilized)
- Memory: ~300-500MB (all contigs + segments)
- Archive: 7.6MB (unchanged)
- Correctness: Should be ✅ 100%
```

### Phase 3: Streaming (if needed)
```
Yeast10 (10 samples, 114M bases):
- Time: ~2-3s (same as Phase 2)
- CPU: 800-1200% (8-12 threads utilized)
- Memory: ~200-300MB (streaming, lower than Phase 2)
- Archive: 7.6MB (unchanged)
- Correctness: Should be ✅ 100%
```

## Key Principles

1. ✅ **Preserve correctness**: Every change must maintain 100% correct output
2. ✅ **Buffer before processing**: Don't compress immediately, buffer first
3. ✅ **Clear phase boundaries**: Segment → Group → Compress → Write
4. ✅ **Test incrementally**: Measure after each phase
5. ✅ **Verify with C++ AGC**: Archives must be readable by C++ AGC

## Testing Strategy

After each optimization:

```bash
# 1. Create archive
time ./target/release/ragc create -o test.agc -k 21 -s 10000 -m 20 -t 15 yeast10/*.fa

# 2. Test extraction (RAGC)
for sample in samples; do
    ./target/release/ragc getset test.agc "$sample" > extracted.fa
    diff original.fa extracted.fa
done

# 3. Test C++ AGC compatibility
/home/erik/agc/bin/agc listset test.agc
for sample in samples; do
    /home/erik/agc/bin/agc getset test.agc "$sample" > extracted_cpp.fa
    diff original.fa extracted_cpp.fa
done

# 4. Compare archive size
ls -lh test.agc  # Should be ~7.6MB

# 5. Measure performance
time ./target/release/ragc create ...  # Wall time
/usr/bin/time -v ./target/release/ragc create ...  # Peak memory
```

## Code Locations

Files to modify:
- `ragc-core/src/compressor_streaming.rs` - Main compression logic
  - Line 1398: Main loop (currently sequential)
  - Line 1418: Segmentation (target for parallelization)
  - Line 1428: Compression (target for parallelization)

## Next Actions

1. [ ] Check that background bash finished (scripts/test_multithreading_performance.sh)
2. [ ] Measure current CPU usage to confirm bottleneck
3. [ ] Implement Phase 1 (parallel compression only)
4. [ ] Test and verify correctness
5. [ ] Measure speedup
6. [ ] If good results, proceed to Phase 2
