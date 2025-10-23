# C++ AGC Exact Implementation Plan

**Goal**: Match C++ AGC architecture precisely using `std::thread`, not Rayon

---

## C++ AGC Architecture (from agc_compressor.cpp)

```cpp
// Priority queue with memory limit
CBoundedPQueue<ContigTask> queue(max_capacity_bytes);

// Main thread: Read FASTAs, push to queue
for (fasta_file : files) {
    while (contig = read_next()) {
        queue.push(ContigTask { sample, contig_name, contig_data });
    }
}

// Worker threads: Pop from queue, process, write immediately
for (worker_id in 0..num_threads-1) {
    std::thread([&]() {
        ZSTD_CCtx* ctx = ZSTD_createCCtx();  // Reused!

        while (task = queue.pop()) {
            segments = split_contig(task.contig_data, splitters);

            for (segment : segments) {
                // Find which group this segment belongs to
                group_id = get_or_create_group(segment.front_kmer, segment.back_kmer);

                // Add to group (groups are shared across workers with mutex)
                group = groups[group_id];
                if (pack = group.add_segment(segment)) {
                    // Group full, write immediately
                    archive.write_pack(pack);  // Thread-safe
                }
            }
        }

        ZSTD_freeCCtx(ctx);
    });
}

// All threads join
// Flush remaining groups
```

---

## Key Insights

1. **Contig-level parallelism**, not segment-level
2. **Shared groups** with mutex protection (Arc<Mutex<HashMap>>)
3. **Immediate writes** when packs are ready
4. **No batching** - process as data arrives
5. **Bounded queue** limits memory (only N contigs in flight)

---

## Rust Implementation

```rust
use std::sync::{Arc, Mutex};
use std::thread;
use crossbeam::channel::bounded;

struct ContigTask {
    sample_name: String,
    contig_name: String,
    sequence: Contig,
}

fn add_fasta_files_cpp_agc_style(
    &mut self,
    fasta_paths: &[(String, &Path)],
    splitters: &HashSet<u64>,
) -> Result<()> {
    // Bounded queue (capacity = number of contigs, not bytes)
    // C++ AGC uses memory-based, but contig count is simpler in Rust
    let queue_capacity = self.config.num_threads * 4;  // 4 contigs per thread
    let (task_tx, task_rx) = bounded::<ContigTask>(queue_capacity);

    // Shared state (groups and archive) - Arc<Mutex> for thread safety
    let groups = Arc::new(Mutex::new(HashMap::<SegmentGroupKey, GroupWriter>::new()));
    let archive = Arc::new(Mutex::new(std::mem::take(&mut self.archive)));
    let collection = Arc::new(Mutex::new(std::mem::take(&mut self.collection)));
    let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));

    // Main thread: Read FASTAs and push to queue
    let reader_handle = thread::spawn(move || {
        for (sample_name, fasta_path) in fasta_paths {
            let reader = GenomeIO::open(fasta_path)?;
            while let Some((contig_name, sequence)) = reader.read_contig()? {
                task_tx.send(ContigTask {
                    sample_name: sample_name.clone(),
                    contig_name,
                    sequence,
                })?;
            }
        }
        drop(task_tx);  // Signal workers
        Ok(())
    });

    // Worker threads: Pop from queue and process
    let worker_handles: Vec<_> = (0..self.config.num_threads - 1)
        .map(|_| {
            let task_rx = task_rx.clone();
            let groups = groups.clone();
            let archive = archive.clone();
            let collection = collection.clone();
            let next_group_id = next_group_id.clone();
            let splitters = splitters.clone();
            let config = self.config.clone();

            thread::spawn(move || {
                // Per-thread reused buffers (like C++ AGC)
                let mut lz_buffer = Vec::new();

                while let Ok(task) = task_rx.recv() {
                    // Split contig into segments
                    let segments = split_at_splitters(
                        &task.sequence,
                        &splitters,
                        config.kmer_length as usize,
                        config.segment_size as usize,
                    );

                    for (seg_idx, segment) in segments.into_iter().enumerate() {
                        let seg_info = prepare_segment_info(
                            &config,
                            &task.sample_name,
                            &task.contig_name,
                            seg_idx,
                            segment.data,
                            segment.front_kmer,
                            segment.back_kmer,
                        )?;

                        // Lock groups, find/create group for this segment
                        let mut groups_lock = groups.lock().unwrap();
                        let group = groups_lock.entry(seg_info.key.clone())
                            .or_insert_with(|| {
                                let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                                GroupWriter::new(gid, ...)
                            });

                        // Add segment, get pack if ready
                        if let Some(pack) = group.add_segment(seg_info, &config)? {
                            // Release groups lock before writing
                            drop(groups_lock);

                            // Write pack immediately
                            let mut archive_lock = archive.lock().unwrap();
                            let mut collection_lock = collection.lock().unwrap();
                            write_pack_to_archive(
                                &mut *archive_lock,
                                &mut *collection_lock,
                                pack,
                            )?;
                        }
                        // groups_lock released here if not already dropped
                    }
                    // Contig processed, memory freed
                }
                Ok(())
            })
        })
        .collect();

    // Wait for reader
    reader_handle.join()??;

    // Wait for all workers
    for handle in worker_handles {
        handle.join()??;
    }

    // Flush remaining groups
    let groups_lock = groups.lock().unwrap();
    for (_, group) in groups_lock.iter() {
        if let Some(pack) = group.flush()? {
            write_pack_to_archive(&mut *archive.lock().unwrap(), pack)?;
        }
    }

    // Restore state
    self.archive = Arc::try_unwrap(archive).unwrap().into_inner().unwrap();
    self.collection = Arc::try_unwrap(collection).unwrap().into_inner().unwrap();
    self.next_group_id = next_group_id.load(Ordering::SeqCst);

    Ok(())
}
```

---

## Memory Profile (Expected)

| Component | Memory | Notes |
|-----------|--------|-------|
| Queue (contigs) | ~48 MB | 24 contigs × 2 MB each |
| Shared groups | ~50 MB | Active GroupWriters |
| Per-thread buffers | ~6 MB | 6 threads × 1 MB |
| Splitters | 91 MB | Same as before |
| Archive | ~50 MB | Write buffers |
| **Total** | **~245 MB** | Close to C++ AGC! |

---

## Key Differences from Rayon Approach

1. ✅ **No Vec<PreparedSegment>** - process on the fly
2. ✅ **No Rayon overhead** - simple `std::thread`
3. ✅ **Bounded queue** - limits memory
4. ✅ **Immediate writes** - no Vec<CompressedPack>
5. ✅ **Shared groups** - one GroupWriter per key across all workers

---

## Implementation Steps

1. Create `ContigTask` struct
2. Implement bounded queue with `crossbeam::channel`
3. Reader thread pushes contigs to queue
4. Worker threads with Arc<Mutex<HashMap<GroupWriter>>>
5. Immediate pack writing with mutex
6. Test and measure

This matches C++ AGC's architecture exactly.
