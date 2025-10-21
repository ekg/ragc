# C++ AGC Architecture: Per-Group Buffering Analysis

## Problem Summary

Our initial Rust implementation had catastrophic performance issues:
- **Real-time writer approach**: 12+ minutes runtime, 108% CPU
- **Arc<DashMap> buffering**: 4.6GB memory usage
- **Target (C++ AGC)**: 61s runtime, 437% CPU, 103MB memory

## Key Discovery: Per-Group Buffering with Batched Writes

Analysis of C++ AGC source code (`~/agc/src/common/segment.cpp`) reveals the critical architectural pattern:

### 1. Per-Group Mutex Locking

Each `CSegment` object (representing one segment group) has its own mutex:

```cpp
uint32_t CSegment::add(const contig_t& s, ZSTD_CCtx* zstd_cctx, ZSTD_DCtx* zstd_dctx)
{
    lock_guard<mutex> lck(mtx);  // Line 36: PER-GROUP MUTEX
    // ...
}
```

**Impact**: Multiple threads can write to different groups concurrently without blocking each other. Only threads writing to the *same* group contend for the lock.

### 2. Buffering Strategy

Segments buffer in a vector (`v_lzp`) until 50 accumulate:

```cpp
if (v_lzp.size() == contigs_in_pack)  // contigs_in_pack = 50
{
    store_in_archive(v_lzp, zstd_cctx);  // BATCHED WRITE
    v_lzp.clear();
}

// ... LZ differential encoding ...

v_lzp.emplace_back(move(delta));  // BUFFER segment
```

**Impact**:
- For 93,000 segments → ~1,860 archive write operations (93K/50)
- vs. our broken real-time writer: 93,000 write operations
- **50x reduction in write overhead**

### 3. Archive Write Pattern

The archive write happens within the per-group mutex:
- Thread acquires group lock
- Adds segment to buffer
- If buffer full (50 segments), compresses pack and writes to archive
- Releases lock
- Other threads continue working on different groups

## Why Our Approaches Failed

### ❌ Approach 1: Arc<DashMap> with Post-Processing
```rust
// Workers buffer ALL segments in shared map
segment_map.insert(key, segments);

// Main thread writes everything at the end
for (key, segments) in segment_map.iter() {
    for segment in segments {
        write_to_archive(segment);  // 93,000 writes!
    }
}
```

**Problem**: Buffers all 93,000 segments in memory = 4.6GB

### ❌ Approach 2: Real-Time Channel-Based Writing
```rust
// Workers send segments via channel
segment_tx.send(segment)?;

// Main thread receives and writes immediately
loop {
    match segment_rx.try_recv() {
        Ok(seg) => write_to_archive(seg),  // SINGLE-THREADED BOTTLENECK
        Err(_) => std::thread::yield_now(),  // BUSY-WAIT SPIN
    }
}
```

**Problems**:
- All 93,000 segments serialized through single main thread
- Lost all parallelism (108% CPU vs 437% target)
- Busy-wait spin loop wastes CPU cycles
- 93,000 archive writes instead of ~1,860

## Correct Rust Implementation Strategy

### Architecture

```rust
// Per-group metadata with internal buffer
struct GroupMetadata {
    pending_segments: Vec<SegmentInfo>,  // Buffer up to 50
    // ... other metadata ...
}

// Concurrent map with per-group mutexes (implicit in DashMap)
type GroupMap = Arc<DashMap<SegmentGroupKey, Mutex<GroupMetadata>>>;

// Worker thread
fn process_contig(group_map: &GroupMap, segment: SegmentInfo) {
    let key = segment.group_key();

    // Get or create group (per-group locking)
    let mut group = group_map.entry(key).or_insert_with(|| {
        Mutex::new(GroupMetadata::new())
    });

    let mut guard = group.lock();
    guard.pending_segments.push(segment);

    // Flush when buffer reaches 50
    if guard.pending_segments.len() >= 50 {
        let pack = std::mem::take(&mut guard.pending_segments);
        drop(guard);  // Release lock before I/O
        write_pack_to_archive(pack);  // Batched write!
    }
}
```

### Key Properties

1. **Bounded memory**: Each group buffers max 50 segments
   - For yeast235 (~5,000 groups): 5,000 × 50 × ~1KB = ~250MB worst case
   - In practice: groups flush as they fill, much less memory

2. **True parallelism**: Workers write to different groups concurrently
   - No single-threaded bottleneck
   - Lock contention only when multiple threads hit same group
   - Expected CPU utilization: 400%+

3. **Reduced archive writes**: ~1,860 operations instead of 93,000
   - Each write handles 50 compressed segments
   - Lower I/O overhead

4. **Clean final flush**: After all workers finish, flush remaining segments
   - Each group may have 0-49 remaining segments
   - Write final partial packs

## Performance Expectations

| Metric | Broken Real-time | Target C++ AGC | Expected with Fix |
|--------|------------------|----------------|-------------------|
| Time | 12+ min | 61s | ~60-80s |
| CPU | 108% | 437% | 400%+ |
| Memory | ~340MB | 103MB | <200MB |
| Archive writes | 93,000 | ~1,860 | ~1,860 |

## Implementation Checklist

- [ ] Add `pending_segments: Vec<SegmentInfo>` to group metadata
- [ ] Use `Arc<DashMap<SegmentGroupKey, Mutex<GroupMetadata>>>` for concurrent access
- [ ] Workers write to group buffers with per-group locking
- [ ] Flush pack of 50 segments when buffer fills
- [ ] Compress pack before writing to archive
- [ ] Final flush of partial packs after workers finish
- [ ] Track segment metadata in collection
- [ ] Test memory usage (<200MB target)
- [ ] Test runtime (~60s target)
- [ ] Test CPU utilization (400%+ target)
- [ ] Verify round-trip compression correctness
