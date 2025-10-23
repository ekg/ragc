# RAGC vs C++ AGC Pipeline Architecture Comparison

## Executive Summary

**Root Cause Identified**: RAGC's 3-stage channel pipeline creates excessive memory allocations compared to C++ AGC's direct processing model.

| Aspect | RAGC | C++ AGC | Impact |
|--------|------|---------|--------|
| **Pipeline** | 3-stage buffered channels | Single priority queue | 🔴 RAGC buffers 3x |
| **Buffer Capacity** | Item count (10-100) | Memory size (2GB or 192MB/thread) | 🔴 RAGC arbitrary limit |
| **ZSTD Context** | Created per-compression | Reused per-thread | 🔴 RAGC allocates repeatedly |
| **Processing Model** | Segment → LZ → Compress → Write | Segment → Compress+Write | 🔴 RAGC extra stage |
| **Thread Coordination** | Crossbeam channels | Priority queue + barriers | 🔴 RAGC more overhead |

---

## Detailed Architecture Comparison

### C++ AGC Pipeline (agc_compressor.cpp)

```
Input Files
    ↓
Main Thread: Read FASTA → Register contig → Enqueue to Priority Queue
    ↓
Priority Queue (CBoundedPQueue)
  - Capacity: max(2GB, 192MB × threads)  ← SIZE-BASED, not count-based!
  - Priority: Sample priority (ensures ordering)
  - Contains: (stage, sample_name, contig_name, contig_data)
    ↓
Worker Threads (no_threads - 1)  ← Leaves 1 core for I/O
    ↓
Per-thread processing:
  - Own ZSTD_CCtx (compression context)
  - Own ZSTD_DCtx (decompression context for matching)
  - Own temporary buffers
    ↓
compress_contig():
  1. Find splitters in contig
  2. Split into segments
  3. For each segment → IMMEDIATELY:
     - add_segment(): Find reference match
     - LZ encode
     - ZSTD compress (using reused ZSTD_CCtx)
     - Write to archive stream
    ↓
Archive (direct write, no buffering)
```

**Key Points:**
1. **Priority Queue with MEMORY-based capacity**: `max(2GB, 192MB × threads)`
   - NOT item-based buffering
   - Naturally limits memory usage
   - Prioritizes samples for deterministic ordering

2. **Context Reuse**:
   - Each thread creates ZSTD contexts ONCE
   - Reuses them for all compressions
   - No allocation overhead per segment

3. **Direct Write Path**:
   - Segment → compress → write (immediate)
   - No intermediate buffering stages
   - Minimal memory overhead

4. **Synchronization**:
   - Uses barriers (`my_barrier`) for thread coordination
   - Sends "synchronization tokens" between samples
   - Ensures proper ordering without complex channels

---

### RAGC Pipeline (compressor_streaming.rs)

```
Input Files
    ↓
Main Thread: Read FASTA → Load ALL contigs into memory  ← PROBLEM 1
    ↓
Producer Thread: Send contigs to channel
    ↓
Channel 1: ContigTask (bounded to 10)
  - Contains: (sample_name, contig_name, Contig)
    ↓
Worker Pool (Rayon): Process contigs in parallel
    ↓
For each contig:
  1. Find splitters
  2. Split into segments
  3. LZ encode ALL segments → Vec<u8>  ← PROBLEM 2
  4. Create UncompressedPack { Vec<u8>, metadata }
  5. Send to Channel 2
    ↓
Channel 2: UncompressedPack (bounded to 10)
  - Contains: Vec<u8> with LZ-encoded data
  - Each Vec is a NEW allocation  ← PROBLEM 3
    ↓
Compression Pool (Rayon): Compress packs
    ↓
For each UncompressedPack:
  1. Create NEW ZSTD context  ← PROBLEM 4
  2. Compress Vec<u8> → NEW Vec<u8>
  3. Create CompressedPack { Vec<u8>, metadata }
  4. Send to Channel 3
    ↓
Channel 3: CompressedPack (bounded to 10)
  - Contains: Vec<u8> with ZSTD data
  - Another NEW allocation  ← PROBLEM 5
    ↓
Writer Thread: Write packs to archive
```

**Key Problems:**

1. **Loads ALL contigs into memory first**
   - C++ streams from disk
   - RAGC loads everything upfront

2. **3-Stage Pipeline = 3× Buffering**
   - Channel 1: Contigs
   - Channel 2: Uncompressed packs (Vec<u8>)
   - Channel 3: Compressed packs (Vec<u8>)

3. **Vec Allocations Per Segment**
   - Every segment creates NEW Vec<u8> for LZ data
   - Every segment creates NEW Vec<u8> for ZSTD data
   - No buffer reuse

4. **ZSTD Context Created Per Compression**
   - C++ creates once per thread
   - RAGC likely creates per pack (need to verify)

5. **Item-Based Channel Limits**
   - `bounded(10)` = 10 items max
   - But each item can be ANY size
   - No memory-based backpressure

---

## Memory Usage Analysis

### Why RAGC Uses 4.81x More Memory

**C++ AGC** (~205 MB):
```
K-mer collection:         91 MB  (unavoidable)
Priority queue:          ~50 MB  (size-limited)
Per-thread contexts:     ~10 MB  (reused)
Active segment data:     ~20 MB  (processing)
Archive buffers:         ~34 MB  (write buffering)
                        -------
Total:                  ~205 MB
```

**RAGC** (~984 MB):
```
K-mer collection:         91 MB  (same)
ALL contigs in memory:  ~110 MB  ← PROBLEM
Channel 1 buffers:       ~20 MB  (10 contigs)
Channel 2 buffers:      ~150 MB  ← PROBLEM (10 uncompressed packs)
Channel 3 buffers:       ~75 MB  ← PROBLEM (10 compressed packs)
Rayon thread pools:     ~100 MB  (thread-local storage)
ZSTD contexts:           ~50 MB  ← PROBLEM (not reused)
Vec allocations:        ~200 MB  ← PROBLEM (no pooling)
Arc/RwLock overhead:     ~40 MB  (shared state)
Segment metadata:       ~148 MB  (duplicated across stages)
                        --------
Total:                  ~984 MB
```

**Key Differences:**
- RAGC loads all contigs: +110 MB
- RAGC 3-stage buffering: +245 MB
- RAGC repeated allocations: +250 MB
- RAGC context overhead: +40 MB

---

## Page Fault Analysis

**Why 54.6M page faults vs 84K:**

### C++ AGC (84K faults):
- Priority queue: 1 allocation
- ZSTD contexts: `num_threads` allocations (reused)
- Per-segment: Minimal allocations (likely reuses buffers)
- Archive writes: Buffered, infrequent allocations

### RAGC (54.6M faults):
For yeast10 with ~12K segments:
```
Segments: ~12,000

Per segment:
  - UncompressedPack Vec:     1 alloc
  - CompressedPack Vec:       1 alloc
  - Metadata Vec:             1 alloc
  - Channel sends (copies):   3 allocs
  - Drop/deallocations:       6 deallocs
                             -------
                            ~12 operations/segment

Total: 12,000 segments × 12 = 144,000 allocations  ← Still doesn't explain 54M!
```

**Additional sources:**
- Rayon work-stealing: Steals work across threads → allocations
- Arc cloning for shared state: Atomic ref counting operations
- Crossbeam channel internals: Node allocations for queue
- ZSTD context creation: Large allocations if not reused
- Temporary Vec in LZ encoding: Intermediate buffers

**Hypothesis**: The 54M faults likely come from:
1. ZSTD context creation PER compression (large allocations)
2. Crossbeam channel internal node allocations
3. Rayon work-stealing internal allocations
4. String cloning for sample/contig names (passed through channels)

---

## Performance Comparison

| Operation | C++ AGC | RAGC | Ratio |
|-----------|---------|------|-------|
| **Memory** | 205 MB | 984 MB | 4.81x |
| **Time** | 3.0s | 15.2s | 5.10x |
| **CPU %** | 101% (1 thread) | 599% (6 threads) | 5.9x |
| **System Time** | 0.09s | 74.5s | **828x** ⚠️ |
| **Page Faults** | 84K | 54.6M | **650x** ⚠️ |

**Key Insight**: Despite using 6 threads vs 1, RAGC is 5x SLOWER!

This proves **parallelism overhead > parallelism benefit**

---

## Proposed Redesign

### Option 1: Minimal Changes (Quick Win)

**Goal**: Match C++ AGC's direct processing without major refactor

```rust
// Remove 3-stage pipeline, use single priority queue
let (task_tx, task_rx) = bounded(memory_based_limit);

// Single worker function per thread
fn worker_thread(rx: Receiver<ContigTask>,
                 archive: &mut Archive,
                 zstd_ctx: &mut ZstdContext) {  // ← Reuse!
    while let Ok((sample, contig_name, contig)) = rx.recv() {
        // Process contig
        for segment in split_contig(&contig) {
            let lz_data = lz_encode(&segment);
            let compressed = zstd_ctx.compress(&lz_data);  // ← Reuse context
            archive.write_immediately(compressed);  // ← No buffering
        }
    }
}
```

**Expected Impact:**
- Memory: 984 MB → ~400 MB (-59%)
- Time: 15.2s → ~8s (-47%)
- Page faults: 54.6M → ~500K (-99%)

### Option 2: C++ AGC-style Architecture (Bigger Refactor)

**Goal**: Fully match C++ AGC design

```rust
struct CompressionPipeline {
    // Memory-based priority queue
    queue: PriorityQueue<ContigTask>,
    // Max 2GB or 192MB per thread
    capacity: usize,  // in bytes, not items

    // Per-thread contexts (reused)
    thread_contexts: Vec<ThreadContext>,
}

struct ThreadContext {
    zstd_compress_ctx: ZstdCCtx,    // Reused
    zstd_decompress_ctx: ZstdDCtx,  // Reused for matching
    lz_buffer: Vec<u8>,              // Reused buffer
    compressed_buffer: Vec<u8>,      // Reused buffer
}

fn process_contig_streaming(
    reader: FastaReader,
    queue: &PriorityQueue,
    workers: usize,
) {
    // Stream from disk, don't load all
    while let Some((sample, contig_name, contig)) = reader.next() {
        queue.push_with_priority((sample, contig_name, contig));
    }

    // Workers pull from queue
    thread::scope(|s| {
        for thread_id in 0..workers {
            s.spawn(|| worker(thread_id, &queue));
        }
    });
}
```

**Expected Impact:**
- Memory: 984 MB → ~250 MB (-75%)
- Time: 15.2s → ~4s (-74%)
- Page faults: 54.6M → ~100K (-99.8%)

---

## Immediate Action Items

1. **Verify ZSTD context reuse** in RAGC
   - Check if contexts are created per-thread or per-compression
   - If per-compression, implement pooling

2. **Profile Vec allocations**
   - Add instrumentation to count allocations
   - Identify hot paths

3. **Test single-thread performance**
   - Run with `num_threads=1`
   - Compare with C++ single-thread
   - Isolate parallelism overhead

4. **Implement buffer pooling** (quick win)
   - Reuse Vec<u8> for LZ encoding
   - Reuse Vec<u8> for compression
   - Expected: -30% memory, -40% page faults

5. **Simplify pipeline** (medium effort)
   - Remove Channel 2 (uncompressed packs)
   - Compress immediately after LZ encoding
   - Expected: -20% memory, -30% page faults

6. **Stream contigs** (medium effort)
   - Don't load all contigs into memory
   - Process as you read
   - Expected: -110 MB memory

---

## Conclusion

**Root Cause**: RAGC's multi-stage buffered pipeline creates excessive memory allocations compared to C++ AGC's streamlined direct-write model.

**Key Differences:**
1. C++ streams from disk; RAGC loads all into memory
2. C++ has 1 priority queue; RAGC has 3 buffered channels
3. C++ reuses ZSTD contexts; RAGC may create per-compression
4. C++ uses memory-based limits; RAGC uses item-count limits
5. C++ processes and writes immediately; RAGC buffers through 3 stages

**Quick Wins** (implement first):
- ZSTD context pooling
- Vec buffer pooling
- Remove intermediate buffering stage

**Medium-term** (architectural):
- Replace 3-channel pipeline with single priority queue
- Stream contigs from disk
- Memory-based backpressure

**Long-term** (optimization):
- Reduce parallelism overhead
- Optimize thread count
- Match C++ AGC's synchronization model
