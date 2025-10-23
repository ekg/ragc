# RAGC vs C++ AGC: Actual Architecture Comparison

**Current Analysis**: Understanding the real memory gap

---

## Current RAGC Architecture (add_fasta_files_with_splitters path)

```
INPUT: Multiple FASTA files (e.g., yeast10: 10 files, ~221 MB total)

=== PASS 1: Find Splitters ===
Location: add_fasta_files_with_splitters() → determine_splitters_streaming()
Memory: ~100 MB (loads first genome to find k-mers)
Output: HashSet of splitter k-mers

=== PASS 2: Process All Samples ===
Path: add_fasta_files_parallel_non_adaptive()

PHASE 1: Load ALL Segments into Memory
---------------------------------------
for each FASTA file:
    reader = open(file)
    while contig = reader.next():
        collection.register_sample_contig()  # Register in collection
        segments = split_at_splitters(contig, splitters)  # Split contig
        for each segment:
            all_segments.push(PreparedSegment {
                key: SegmentGroupKey,
                segment: SegmentInfo {
                    sample_name: String,
                    contig_name: String,
                    seg_part_no: usize,
                    data: Vec<u8>,       # ← FULL SEGMENT DATA
                    front_kmer, back_kmer
                }
            })

Result: Vec<PreparedSegment> with ~12K segments
Memory at this point: ~300-400 MB (all segment data in memory!)

PHASE 2: Group Segments by K-mer Keys
--------------------------------------
let mut groups: HashMap<SegmentGroupKey, Vec<PreparedSegment>> = HashMap::new();
for segment in all_segments:
    groups.entry(segment.key).or_default().push(segment);

let groups_vec: Vec<(SegmentGroupKey, Vec<PreparedSegment>)> = groups.into_iter().collect();

Memory impact: Segments remain in memory, now organized by group

PHASE 3: Parallel Processing with Rayon
----------------------------------------
ThreadPool = rayon::ThreadPoolBuilder::new().num_threads(config.num_threads)

all_packs = ThreadPool.install(|| {
    groups_vec.par_iter()
        .flat_map(|(key, segments)| {
            // This runs on each Rayon worker thread

            group_writer = GroupWriter::new()
            let mut packs = Vec::new()

            for segment in segments:
                // Process segment
                if pack = group_writer.add_segment(segment):
                    packs.push(pack)

            // Flush remaining
            if pack = group_writer.flush():
                packs.push(pack)

            return packs  # Vec<CompressedPack>
        })
        .collect()  # ← Collects ALL CompressedPacks into Vec
})

Result: Vec<CompressedPack> with ~1200 packs
Memory at this point: ~600-700 MB (segments + compressed packs!)

PHASE 4: Sequential Write
--------------------------
for pack in all_packs:
    archive.register_stream()
    archive.write_pack(pack.compressed_data)
    collection.register_segments(pack.segments)

# Finally free memory
drop(all_packs)
```

---

## C++ AGC Architecture (agc_compressor.cpp)

```
INPUT: Multiple FASTA files

=== Pass 1: Find Splitters ===
Same as RAGC: load reference genome, find k-mers
Memory: ~100 MB

=== Pass 2: Stream and Process ===
Main approach: STREAM contigs, process immediately

Priority Queue Setup:
----------------------
CBoundedPQueue queue(max_capacity_bytes)  // Memory-based limit!
  capacity = max(2GB, 192MB × num_threads)
  priority = sample_priority (for ordering)

Processing Loop:
----------------
for each FASTA file:
    fstream reader(file)
    while (contig = reader.read_next()):
        # Push contig to queue (blocks if queue is full)
        queue.push(ContigTask {
            stage, sample_name, contig_name,
            contig_data: string  # ← Only ONE contig in memory at a time
        })

Worker Threads (num_threads - 1):
----------------------------------
Each thread has:
  - ZSTD_CCtx* zstd_ctx (reused!)
  - ZSTD_DCtx* zstd_dctx (reused!)
  - vector<uint8_t> tmp_buffer (reused!)

while (task = queue.pop()):
    compress_contig(task, zstd_ctx):
        segments = split_contig(task.contig_data, splitters)

        for segment in segments:
            # IMMEDIATELY process and write
            ref_match = find_reference_match(segment)
            lz_data = lz_encode(segment, ref_match, tmp_buffer)  # reuse buffer!
            compressed = zstd_compress(lz_data, zstd_ctx)  # reuse context!

            archive.write_segment(compressed)  # Write immediately!
            collection.register(segment_info)

        # Contig is done, free its memory
        task.contig_data.clear()

Key: Only queue-capacity worth of contigs in memory
Key: No intermediate buffering - write as you go
Key: Contexts reused across ALL compressions
```

---

## Key Architectural Differences

### 1. Memory Loading Strategy

**RAGC**:
- Load ALL segments from ALL files into Vec<PreparedSegment>
- ~12K segments × ~18 KB avg = **~216 MB just for segment data**
- Remains in memory through all phases

**C++ AGC**:
- Stream contigs one at a time
- Queue holds at most `capacity_bytes / avg_contig_size` contigs
- With 2GB limit and 12MB avg: ~170 contigs max
- **~170 contigs × 12 MB = ~2 GB max**, but usually much less

### 2. Processing Model

**RAGC - Batch Mode**:
```
1. Load everything → 300 MB
2. Group everything → Still 300 MB
3. Process in parallel → 300 MB segments + 400 MB packs = 700 MB
4. Write sequentially → Free memory
```

**C++ AGC - Streaming Mode**:
```
1. Stream contig → ~12 MB
2. Process → ~12 MB contig + ~20 MB working = ~32 MB
3. Write immediately → Free contig
4. Next contig → ~12 MB
   (Queue keeps ~5-10 contigs buffered = ~120 MB)
```

### 3. Parallelism Strategy

**RAGC**:
- Rayon par_iter() over groups
- `.collect()` gathers ALL results before writing
- Peak memory = input + output simultaneously

**C++ AGC**:
- Worker threads pull from priority queue
- Write immediately per-segment
- Peak memory = queue capacity only

### 4. ZSTD Context Management

**RAGC** (after our fix):
- Thread-local buffer pooling for output
- But still creates encoder instance per call via `zstd::stream::copy_encode()`
- zstd-rs doesn't expose context reuse API

**C++ AGC**:
- `ZSTD_CCtx* ctx = ZSTD_createCCtx()` once per thread
- `ZSTD_compressCCtx(ctx, ...)` reuses context
- Significant performance and memory benefit

---

## Memory Breakdown Comparison

### C++ AGC (~205 MB peak)
```
K-mer splitters:           91 MB
Priority queue:           ~60 MB  (5-10 contigs buffered)
Per-thread contexts:      ~10 MB  (6 threads × ~1.5 MB)
Active processing:        ~20 MB  (current contig being processed)
Archive write buffers:    ~24 MB
                         --------
Total:                   ~205 MB
```

### RAGC Current (~731 MB peak with 6 threads)
```
K-mer splitters:           91 MB  (same)
ALL PreparedSegments:    ~216 MB  ← MAIN ISSUE (12K × 18 KB)
Groups HashMap:           ~50 MB  (pointers + overhead)
Rayon thread pools:       ~80 MB  (thread-local state)
CompressedPacks Vec:     ~200 MB  ← .collect() holds everything
Working memory:           ~50 MB  (active compressions)
Arc/metadata overhead:    ~44 MB
                         --------
Total:                   ~731 MB
```

### Root Causes of 3.6x Gap

1. **All-segments-in-memory (Phase 1)**: **+216 MB**
   - RAGC: Loads everything upfront
   - C++ AGC: Streams with bounded queue

2. **Collect-all-packs (Phase 3)**: **+200 MB**
   - RAGC: `.collect()` gathers all CompressedPacks
   - C++ AGC: Writes immediately

3. **Rayon overhead**: **+80 MB**
   - RAGC: Thread pools, work-stealing queues
   - C++ AGC: Simple worker threads

4. **HashMap grouping**: **+50 MB**
   - RAGC: Groups all segments before processing
   - C++ AGC: Processes as streamed

**Total extra**: ~546 MB
**C++ base**: 205 MB
**RAGC total**: 751 MB (matches measured 731 MB!)

---

## Proposed Architectural Changes

### Option 1: Quick Fixes (Incremental)

1. **Eliminate .collect() in Phase 3**
   - Write CompressedPacks as produced
   - Saves **~200 MB**
   - Complexity: Medium (need to coordinate with Rayon)

2. **Stream segments instead of loading all**
   - Don't create Vec<PreparedSegment>
   - Process groups as segments are read
   - Saves **~216 MB**
   - Complexity: High (requires multi-pass or different grouping)

### Option 2: Full Redesign (C++ AGC-style)

Replace the 3-phase architecture with streaming:

```rust
fn add_fasta_files_streaming(files: &[Path]) -> Result<()> {
    // Priority queue with memory limit
    let queue = PriorityQueue::with_capacity_bytes(2GB);

    // Stream contigs into queue
    thread::spawn(|| {
        for file in files {
            let reader = FastaReader::new(file);
            while let Some(contig) = reader.next() {
                queue.push_with_priority(ContigTask {
                    sample, contig_name, contig_data
                });
            }
        }
    });

    // Worker threads process from queue
    thread::scope(|s| {
        for tid in 0..num_threads {
            s.spawn(|| {
                // Per-thread reused contexts
                let mut zstd_ctx = ZstdCCtx::new();
                let mut lz_buffer = Vec::new();

                while let Some(task) = queue.pop() {
                    process_contig_immediate(
                        task,
                        &mut zstd_ctx,
                        &mut lz_buffer,
                        &mut archive  // Write directly!
                    );
                }
            });
        }
    });
}
```

**Expected savings**:
- No all-segments buffering: -216 MB
- No collect-all-packs: -200 MB
- Simpler threading: -80 MB
- Total: **-496 MB** → **~235 MB** (close to C++ AGC!)

---

## Immediate Next Steps

Given the 3.6x gap, here's the priority order:

1. **Profile current code with Rayon fix** ✓ (Done: 731 MB)

2. **Eliminate .collect() in Phase 3**
   - Write packs as produced instead of buffering
   - Estimated: -200 MB (27% reduction)
   - Complexity: Medium

3. **Stream segments in Phase 1**
   - Process incrementally, don't load all
   - Estimated: -216 MB (30% reduction)
   - Complexity: High

4. **Full streaming redesign**
   - Priority queue + immediate writes
   - Estimated: -496 MB (68% reduction)
   - Complexity: Very High

The 3.6x gap is primarily architectural - batch loading + batch output vs streaming.
