# Streaming Architecture Design (Option B)

**Goal**: Reduce memory from 731 MB to ~235 MB by eliminating batch loading

---

## Current Architecture (Batch Mode)

```
Phase 1: Load ALL segments into Vec
  ↓ (216 MB consumed - ALL segment data in memory)
Phase 2: Group by k-mer keys using HashMap
  ↓ (still 216 MB + grouping overhead)
Phase 3: Rayon par_iter() + .collect()
  ↓ (216 MB segments + 200 MB packs = 416 MB peak!)
Phase 4: Write all packs
```

## New Architecture (Streaming Mode)

```
Main Thread: Stream contigs from FASTA
  ↓
Channel → Worker Pool (group-based routing)
  ↓
Workers: Process assigned groups, maintain buffered state
  ↓
Write Channel → Single Writer Thread
  ↓
Archive (immediate write, no buffering)
```

### Key Changes

1. **No Vec<PreparedSegment>**: Stream and process immediately
2. **Group Writers**: Each worker handles specific k-mer key ranges
3. **Immediate Writes**: Packs written as produced, not collected
4. **Memory Limit**: Use bounded channels (memory-based, not count-based)

---

## Implementation Plan

### Step 1: Streaming Segmentation

Replace:
```rust
// OLD: Load all
let mut all_segments: Vec<PreparedSegment> = Vec::new();
for file in files {
    for contig in read_file(file) {
        let segments = segment(contig);
        all_segments.extend(segments); // ← Accumulates everything
    }
}
```

With:
```rust
// NEW: Stream
let (segment_tx, segment_rx) = crossbeam::channel::bounded(100);
thread::spawn(|| {
    for file in files {
        for contig in read_file(file) {
            let segments = segment(contig);
            for seg in segments {
                segment_tx.send(seg); // ← Send immediately
            }
        }
    }
    drop(segment_tx); // Signal completion
});
```

**Memory saved**: -216 MB (no Vec<PreparedSegment>)

### Step 2: Parallel Group Processing

Instead of grouping ALL segments first, route segments to workers by group key:

```rust
// Spawn workers (one per group key range, or hash-based routing)
for worker_id in 0..num_threads {
    let rx = segment_rx.clone();
    let write_tx = write_channel_tx.clone();

    thread::spawn(move || {
        // Each worker maintains its own groups
        let mut my_groups: HashMap<SegmentGroupKey, GroupWriter> = HashMap::new();

        while let Ok(segment) = rx.recv() {
            // Route: only handle segments for this worker's key range
            if segment.key.hash() % num_threads != worker_id {
                continue; // Wrong worker, skip
            }

            // Get or create group writer
            let group = my_groups.entry(segment.key).or_insert_with(|| {
                GroupWriter::new(assign_group_id())
            });

            // Add segment, may produce packs
            if let Some(pack) = group.add_segment(segment) {
                write_tx.send(pack); // ← Write immediately!
            }
        }

        // Flush remaining
        for (_, group) in my_groups {
            if let Some(pack) = group.flush() {
                write_tx.send(pack);
            }
        }
    });
}
```

**Memory saved**: Workers only hold their subset of groups, not all groups

### Step 3: Immediate Writing

Replace Phase 4's iteration over `all_packs` with streaming writes:

```rust
// Writer thread
let (write_tx, write_rx) = crossbeam::channel::bounded(10);

thread::spawn(move || {
    let mut stream_id_map = HashMap::new();

    while let Ok(pack) = write_rx.recv() {
        // Register stream on-demand
        let stream_id = get_or_register_stream(&mut stream_id_map, &pack);

        // Write immediately
        archive.write_pack(stream_id, &pack.compressed_data);

        // Register segments
        for seg_meta in pack.segments {
            collection.register(seg_meta);
        }
    } // ← pack is dropped immediately after writing!
});
```

**Memory saved**: -200 MB (no Vec<CompressedPack>)

---

## Memory Calculation

### Before (Batch)
```
ALL segments:          216 MB
HashMap grouping:       50 MB
CompressedPacks Vec:   200 MB
Rayon overhead:         80 MB
Other:                 185 MB (k-mers, archive, etc.)
                      -------
Total:                 731 MB
```

### After (Streaming)
```
Channel buffers:        20 MB  (100 segments × ~200 KB)
Worker group state:     50 MB  (per-worker groups)
Write channel:          10 MB  (10 packs × ~1 MB)
Rayon overhead:         70 MB  (reduced - simpler threading)
Other:                 185 MB  (same)
                      -------
Total:                ~335 MB  (closer, but still need more optimization)
```

Wait, this doesn't get us to 235 MB yet. Need additional optimizations:

### Additional Optimizations

1. **Reduce channel buffer sizes** (memory-based limits)
   - Current: 100 segments = ~20 MB
   - New: 10 segments = ~2 MB
   - Savings: -18 MB

2. **Use thread pool instead of Rayon**
   - Simpler threading model
   - Less overhead
   - Savings: -20 MB

3. **Reuse buffers** (per-thread Vec pooling)
   - LZ encoding buffer
   - Compression output buffer
   - Savings: -20 MB

4. **Stream FASTA reading** (already doing this)

**Revised Total**: ~335 - 58 = **~277 MB**

To get to 235 MB, we'd also need to optimize the "Other" category (k-mers, collection, archive state).

---

## Implementation Steps

### Step 1: Add streaming infrastructure
- [x] Identify current batch points
- [ ] Add crossbeam dependency
- [ ] Create segment streaming channel
- [ ] Create write channel

### Step 2: Refactor add_fasta_files_parallel_non_adaptive
- [ ] Replace Vec<PreparedSegment> with channel
- [ ] Modify workers to route by group key
- [ ] Remove .collect(), use for_each() with channel send

### Step 3: Test and measure
- [ ] Build and test correctness
- [ ] Measure memory with `/usr/bin/time -v`
- [ ] Verify C++ AGC compatibility

### Step 4: Additional optimizations (if needed)
- [ ] Memory-based channel limits
- [ ] Per-thread buffer pooling
- [ ] Simplified threading

---

## Risk Assessment

### Challenges

1. **Group Assignment**: Need deterministic group IDs
   - Solution: Use atomic counter

2. **Write Ordering**: Packs may arrive out-of-order
   - Solution: Writer thread reorders or AGC format is order-independent

3. **Channel Deadlock**: If workers block on full channels
   - Solution: Bounded channels with appropriate sizing

4. **Performance**: More thread coordination
   - Solution: Profile and optimize critical paths

### Rollback Plan

If streaming doesn't work:
- Keep changes in separate branch
- Can revert to batch mode
- Or implement hybrid approach

---

## Success Criteria

1. **Memory**: < 300 MB (target: ~235-277 MB)
2. **Correctness**: All tests pass
3. **Compatibility**: C++ AGC can read archives
4. **Performance**: ≤ 25s wall time (allow slight slowdown for memory savings)
