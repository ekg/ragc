# Option B Implementation Roadmap

**Goal**: Reduce memory from 731 MB to ~235-277 MB through full streaming redesign

---

## Current Achievement
- ✅ Fixed Rayon threading: **36% reduction** (1140 MB → 731 MB)
- ✅ Identified root cause: Batch loading architecture
- ✅ Designed streaming solution (Option B)

## Target
- **Memory**: 731 MB → ~235-277 MB (**~68% reduction**)
- **Gap vs C++ AGC**: 3.6x → 1.1-1.3x

---

## Implementation Steps

### Phase 1: Add Streaming Infrastructure

**File**: `ragc-core/src/compressor_streaming.rs`

1. Add imports (already have crossbeam):
```rust
use crossbeam::channel::{bounded, unbounded};
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::Arc;
```

2. Create new method `add_fasta_files_streaming()` that will eventually replace `add_fasta_files_parallel_non_adaptive()`

### Phase 2: Streaming Segmentation (Replace Phase 1)

**Current code (lines 1655-1720)**:
```rust
let mut all_segments: Vec<PreparedSegment> = Vec::new(); // ← 216 MB!
for (sample_name, fasta_path) in fasta_paths {
    let mut reader = GenomeIO::open(fasta_path)?;
    while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
        let segments = split_at_splitters_with_size(...);
        for segment in segments {
            all_segments.push(PreparedSegment { ... }); // ← Accumulates
        }
    }
}
```

**New code**:
```rust
// Create segment streaming channel
let (segment_tx, segment_rx) = bounded::<PreparedSegment>(100);

// Spawn segmentation thread
let segmentation_handle = {
    let fasta_paths = fasta_paths.to_vec();
    let splitters = splitters.clone();
    let config = self.config.clone();
    let mut collection = self.collection.clone(); // For registration

    std::thread::spawn(move || -> Result<(usize, usize)> {
        let mut total_segments = 0;
        let mut total_bases = 0;

        for (sample_name, fasta_path) in &fasta_paths {
            let mut reader = GenomeIO::open(fasta_path)?;

            while let Some((contig_name, sequence)) = reader.read_contig_converted()? {
                collection.register_sample_contig(&sample_name, &contig_name)?;

                let segments = split_at_splitters_with_size(
                    &sequence,
                    &splitters,
                    config.kmer_length as usize,
                    config.segment_size as usize,
                );

                for (seg_idx, segment) in segments.into_iter().enumerate() {
                    total_bases += segment.data.len();
                    total_segments += 1;

                    let seg_info = Self::prepare_segment_info(
                        &config, &sample_name, &contig_name, seg_idx,
                        segment.data, segment.front_kmer, segment.back_kmer,
                    )?;

                    // Send immediately, don't accumulate!
                    segment_tx.send(PreparedSegment {
                        key: seg_info.key.clone(),
                        segment: seg_info,
                    }).context("Failed to send segment")?;
                }
            }
        }

        drop(segment_tx); // Signal completion
        Ok((total_segments, total_bases))
    })
};
```

**Memory saved**: -216 MB (no Vec accumulation)

### Phase 3: Streaming Group Processing (Replace Phase 2 & 3)

**Current code (lines 1737-1885)**:
```rust
// Phase 2: Group all segments
let mut groups: HashMap<SegmentGroupKey, Vec<PreparedSegment>> = HashMap::new();
for segment in all_segments {
    groups.entry(segment.key).or_default().push(segment);
}

// Phase 3: Process with Rayon, collect all packs
let all_packs: Vec<CompressedPack> = pool.install(|| {
    groups_vec.par_iter()
        .flat_map(|(_key, segments)| { ... })
        .collect()  // ← Collects all 200 MB!
});
```

**New code**:
```rust
// Create write channel
let (write_tx, write_rx) = bounded::<CompressedPack>(10);

// Spawn writer thread
let writer_handle = {
    let mut archive = std::mem::take(&mut self.archive);
    let mut collection = std::mem::take(&mut self.collection);

    std::thread::spawn(move || -> Result<(Archive, CollectionV3)> {
        let mut stream_id_map = HashMap::new();
        let archive_version = AGC_FILE_MAJOR * 1000 + AGC_FILE_MINOR;
        let mut packs_written = 0;

        while let Ok(pack) = write_rx.recv() {
            // Register stream on-demand
            let actual_stream_id = if let Some(&id) = stream_id_map.get(&pack.stream_id) {
                id
            } else {
                let stream_name = if pack.stream_id >= 10000 {
                    stream_ref_name(archive_version, (pack.stream_id - 10000) as u32)
                } else {
                    stream_delta_name(archive_version, pack.group_id)
                };
                let id = archive.register_stream(&stream_name);
                stream_id_map.insert(pack.stream_id, id);
                id
            };

            // Write immediately
            archive.write_pack(actual_stream_id, &pack.compressed_data, pack.uncompressed_size)?;

            // Register segments
            for seg_meta in &pack.segments {
                collection.register_segment(
                    pack.group_id,
                    seg_meta.in_group_id,
                    &seg_meta.sample_name,
                    &seg_meta.contig_name,
                    seg_meta.seg_part_no as u32,
                    seg_meta.is_rev_comp,
                    seg_meta.data_len,
                )?;
            }

            packs_written += 1;
        }

        Ok((archive, collection))
    })
};

// Worker threads: Process segments and send packs to writer
let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));
let config = self.config.clone();

let worker_handles: Vec<_> = (0..self.config.num_threads)
    .map(|worker_id| {
        let segment_rx = segment_rx.clone();
        let write_tx = write_tx.clone();
        let next_group_id = next_group_id.clone();
        let config = config.clone();

        std::thread::spawn(move || -> Result<()> {
            // Each worker maintains its own groups
            let mut my_groups: HashMap<SegmentGroupKey, (u32, GroupWriter)> = HashMap::new();

            while let Ok(prepared_seg) = segment_rx.recv() {
                // Route by hash to this worker
                let key_hash = {
                    use std::collections::hash_map::DefaultHasher;
                    use std::hash::{Hash, Hasher};
                    let mut hasher = DefaultHasher::new();
                    prepared_seg.key.hash(&mut hasher);
                    hasher.finish()
                };

                if (key_hash as usize) % config.num_threads != worker_id {
                    continue; // Not for this worker
                }

                // Get or create group
                let (group_id, group_writer) = my_groups.entry(prepared_seg.key.clone())
                    .or_insert_with(|| {
                        let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                        let stream_id = gid as usize;
                        let ref_stream_id = if gid >= 16 {
                            Some(10000 + gid as usize)
                        } else {
                            None
                        };
                        (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
                    });

                // Add segment, may produce pack
                if let Some(pack) = group_writer.add_segment(prepared_seg.segment, &config)? {
                    write_tx.send(pack).context("Failed to send pack to writer")?;
                }
            }

            // Flush all groups
            for (_key, (_gid, mut group_writer)) in my_groups {
                if let Some(pack) = group_writer.prepare_pack(&config)? {
                    write_tx.send(pack).context("Failed to send pack to writer")?;
                }
            }

            Ok(())
        })
    })
    .collect();

// Drop sender so writer knows when workers are done
drop(segment_rx);
drop(write_tx);

// Wait for segmentation to complete
let (total_segments, total_bases) = segmentation_handle.join().unwrap()?;

// Wait for all workers
for handle in worker_handles {
    handle.join().unwrap()?;
}

// Wait for writer and restore archive/collection
let (archive, collection) = writer_handle.join().unwrap()?;
self.archive = archive;
self.collection = collection;

self.next_group_id = next_group_id.load(Ordering::SeqCst);
self.total_segments = total_segments;
self.total_bases_processed = total_bases;
```

**Memory saved**: -200 MB (no Vec<CompressedPack> collection)

### Phase 4: Remove Old Phase 4 Write Loop

Lines 1896-1998 can be deleted - writing now happens in writer thread.

---

## Testing Plan

### 1. Unit Tests
```bash
cargo test --release
```

### 2. Memory Measurement
```bash
/usr/bin/time -v ./target/release/ragc create \
  -o /tmp/test.agc \
  -k 21 -s 10000 -m 20 -v 1 -t 6 \
  /home/erik/scrapy/yeast10_test/*.fa
```

Expected: **< 300 MB** (target: ~235-277 MB)

### 3. C++ AGC Compatibility
```bash
# Extract with C++ AGC
agc listset /tmp/test.agc

# Verify all samples present
agc getset /tmp/test.agc sample1 > /tmp/extracted.fa

# Compare checksums
```

### 4. Performance Benchmark
```bash
# Measure wall time (allow slight slowdown for memory savings)
time ./target/release/ragc create -o /tmp/test.agc ...
```

Expected: **≤ 25s** (vs current 22s)

---

## Rollback Plan

If streaming implementation has issues:

1. **Keep in separate branch**: `feature/streaming-redesign`
2. **Revert if needed**: Can go back to batch mode (731 MB)
3. **Hybrid approach**: Could implement Phase-by-Phase:
   - First: Just eliminate .collect() (-200 MB → 531 MB)
   - Then: Add streaming segmentation (-216 MB → 315 MB)

---

## Success Criteria

1. ✅ **Memory**: < 300 MB (stretch: ~235 MB)
2. ✅ **Correctness**: All unit tests pass
3. ✅ **Compatibility**: C++ AGC can read/write archives
4. ✅ **Performance**: ≤ 25s wall time (22s currently)
5. ✅ **Maintainability**: Code remains clear and documented

---

## Next Session

Start with:
1. Create `add_fasta_files_streaming()` skeleton
2. Implement segmentation thread (Phase 1)
3. Test segmentation correctness
4. Implement worker threads (Phase 2)
5. Test full pipeline
6. Measure and iterate

Estimated implementation time: **2-3 hours** for full streaming redesign
