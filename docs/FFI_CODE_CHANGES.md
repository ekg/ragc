# Exact Code Changes for FFI Integration

## File 1: ragc-core/src/ffi/agc_index.cpp

### Add to extern "C" block (after line 73):

```cpp
// ========== Archive Stream Functions ==========

size_t agc_archive_register_stream(void* archive_ptr, const char* name) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    return (*archive)->RegisterStream(std::string(name));
}

bool agc_archive_add_part(
    void* archive_ptr,
    size_t stream_id,
    const uint8_t* data,
    size_t len,
    uint64_t metadata
) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    std::vector<uint8_t> vec(data, data + len);
    return (*archive)->AddPart(stream_id, vec, metadata);
}
```

## File 2: ragc-core/src/ffi/agc_index.rs

### Add to extern "C" block (after line 40):

```rust
    fn agc_archive_register_stream(
        archive_ptr: *mut std::ffi::c_void,
        name: *const c_char,
    ) -> usize;

    fn agc_archive_add_part(
        archive_ptr: *mut std::ffi::c_void,
        stream_id: usize,
        data: *const u8,
        len: usize,
        metadata: u64,
    ) -> bool;
```

### Add to CppArchive impl block (after line 97):

```rust
    pub fn register_stream(&mut self, name: &str) -> anyhow::Result<usize> {
        let c_name = CString::new(name)?;
        unsafe {
            Ok(agc_archive_register_stream(self.ptr, c_name.as_ptr()))
        }
    }

    pub fn add_part(&mut self, stream_id: usize, data: &[u8], metadata: u64) -> anyhow::Result<()> {
        unsafe {
            if agc_archive_add_part(self.ptr, stream_id, data.as_ptr(), data.len(), metadata) {
                Ok(())
            } else {
                Err(anyhow::anyhow!("Failed to add archive part to stream {}", stream_id))
            }
        }
    }
```

## File 3: ragc-core/src/streaming_compressor_queue.rs

### Change 1: Add imports (after line 9)

```rust
#[cfg(feature = "ffi_cost")]
use crate::ffi::agc_index::{CppArchive, CppCollection, SegmentPlacement};
```

### Change 2: Update struct fields (lines 212-215)

**BEFORE:**
```rust
pub struct StreamingQueueCompressor {
    queue: Arc<MemoryBoundedQueue<ContigTask>>,
    workers: Vec<JoinHandle<Result<()>>>,
    collection: Arc<Mutex<CollectionV3>>,
    splitters: Arc<HashSet<u64>>,
    config: StreamingQueueConfig,
    archive: Arc<Mutex<Archive>>,
```

**AFTER:**
```rust
pub struct StreamingQueueCompressor {
    queue: Arc<MemoryBoundedQueue<ContigTask>>,
    workers: Vec<JoinHandle<Result<()>>>,

    #[cfg(feature = "ffi_cost")]
    collection: Arc<Mutex<CppCollection>>,
    #[cfg(not(feature = "ffi_cost"))]
    collection: Arc<Mutex<CollectionV3>>,

    splitters: Arc<HashSet<u64>>,
    config: StreamingQueueConfig,

    #[cfg(feature = "ffi_cost")]
    archive: Arc<Mutex<CppArchive>>,
    #[cfg(not(feature = "ffi_cost"))]
    archive: Arc<Mutex<Archive>>,

    segment_groups: Arc<Mutex<BTreeMap<SegmentGroupKey, SegmentGroupBuffer>>>,
    group_counter: Arc<AtomicU32>,
    raw_group_counter: Arc<AtomicU32>,
    reference_sample_name: Arc<Mutex<Option<String>>>,

    // Placement buffer for FFI batching
    #[cfg(feature = "ffi_cost")]
    placement_buffer: Arc<Mutex<Vec<SegmentPlacement>>>,
```

### Change 3: Update constructor - Archive creation (lines 270-272)

**BEFORE:**
```rust
        // Create archive
        let mut archive = Archive::new_writer();
        archive.open(output_path)?;
```

**AFTER:**
```rust
        // Create archive
        #[cfg(feature = "ffi_cost")]
        let mut archive = {
            let mut arch = CppArchive::new_writer();
            arch.open(output_path.to_str().unwrap())?;
            arch
        };
        #[cfg(not(feature = "ffi_cost"))]
        let mut archive = {
            let mut arch = Archive::new_writer();
            arch.open(output_path)?;
            arch
        };
```

### Change 4: Update constructor - Collection creation (lines 274-280)

**BEFORE:**
```rust
        // Create collection
        let mut collection = CollectionV3::new();
        collection.set_config(config.segment_size as u32, config.k as u32, None);

        // CRITICAL: Register collection streams FIRST (C++ AGC compatibility)
        // C++ AGC expects collection-samples at stream 0, collection-contigs at 1, collection-details at 2
        collection.prepare_for_compression(&mut archive)?;
```

**AFTER:**
```rust
        // Create collection
        #[cfg(feature = "ffi_cost")]
        let mut collection = {
            let mut coll = CppCollection::new();
            // set_archives() registers streams internally (0-2) and sets config
            coll.set_archives(
                &mut archive,
                config.num_threads as u32,
                1000,  // batch_size (matches C++ AGC default)
                config.segment_size as u32,
                config.k as u32,
            )?;
            coll
        };
        #[cfg(not(feature = "ffi_cost"))]
        let mut collection = {
            let mut coll = CollectionV3::new();
            coll.set_config(config.segment_size as u32, config.k as u32, None);
            coll.prepare_for_compression(&mut archive)?;
            coll
        };
```

### Change 5: Update constructor - Stream registration (lines 282-339)

**AFTER collection creation, BEFORE let collection = Arc::new(Mutex::new(collection)):**

Add conditional stream registration:
```rust
        // Write metadata streams (file_type_info, params, splitters, segment-splitters)
        // These are written AFTER collection streams (0-2) for C++ AGC compatibility

        #[cfg(feature = "ffi_cost")]
        {
            // file_type_info stream
            let mut data = Vec::new();
            let append_str = |data: &mut Vec<u8>, s: &str| {
                data.extend_from_slice(s.as_bytes());
                data.push(0);
            };
            append_str(&mut data, "producer");
            append_str(&mut data, "ragc");
            append_str(&mut data, "producer_version_major");
            append_str(&mut data, &ragc_common::AGC_FILE_MAJOR.to_string());
            append_str(&mut data, "producer_version_minor");
            append_str(&mut data, &ragc_common::AGC_FILE_MINOR.to_string());
            append_str(&mut data, "producer_version_build");
            append_str(&mut data, "0");
            append_str(&mut data, "file_version_major");
            append_str(&mut data, &ragc_common::AGC_FILE_MAJOR.to_string());
            append_str(&mut data, "file_version_minor");
            append_str(&mut data, &ragc_common::AGC_FILE_MINOR.to_string());
            append_str(&mut data, "comment");
            append_str(&mut data, &format!("RAGC v.{}.{}", ragc_common::AGC_FILE_MAJOR, ragc_common::AGC_FILE_MINOR));

            let stream_id = archive.register_stream("file_type_info")?;
            archive.add_part(stream_id, &data, 7)?;

            // params stream
            let mut params_data = Vec::new();
            params_data.extend_from_slice(&(config.k as u32).to_le_bytes());
            params_data.extend_from_slice(&(config.min_match_len as u32).to_le_bytes());
            params_data.extend_from_slice(&50u32.to_le_bytes());
            params_data.extend_from_slice(&(config.segment_size as u32).to_le_bytes());
            let params_stream_id = archive.register_stream("params")?;
            archive.add_part(params_stream_id, &params_data, 0)?;

            // Empty splitters and segment-splitters streams (written in finalize())
            let splitters_stream_id = archive.register_stream("splitters")?;
            archive.add_part(splitters_stream_id, &[], 0)?;

            let seg_splitters_stream_id = archive.register_stream("segment-splitters")?;
            archive.add_part(seg_splitters_stream_id, &[], 0)?;
        }
        #[cfg(not(feature = "ffi_cost"))]
        {
            // Original Rust implementation (lines 282-339)
            // ... (keep existing code)
        }
```

### Change 6: Initialize placement buffer (after line 350)

```rust
        // Placement buffer for batched FFI calls
        #[cfg(feature = "ffi_cost")]
        let placement_buffer = Arc::new(Mutex::new(Vec::new()));
```

And in the return statement:
```rust
        Ok(Self {
            queue,
            workers,
            collection,
            splitters,
            config,
            archive,
            segment_groups,
            group_counter,
            raw_group_counter,
            reference_sample_name,
            map_segments,
            map_segments_terminators,
            #[cfg(feature = "ffi_cost")]
            grouping_engine,
            reference_segments,
            split_offsets,
            sample_priorities,
            next_priority,
            #[cfg(feature = "ffi_cost")]
            placement_buffer,
        })
```

### Change 7: Add helper method for flushing placements

Add new method to impl block (before finalize()):

```rust
    #[cfg(feature = "ffi_cost")]
    fn flush_placement_buffer(&self) -> Result<()> {
        let mut buffer = self.placement_buffer.lock().unwrap();
        if !buffer.is_empty() {
            let mut coll = self.collection.lock().unwrap();
            coll.add_segments_placed(&buffer)
                .context("Failed to flush placement buffer")?;
            buffer.clear();
        }
        Ok(())
    }
```

### Change 8: Update flush_pack() - Reference segment registration (line 859)

**BEFORE:**
```rust
        // Register reference in collection with in_group_id = 0
        {
            let mut coll = collection.lock().unwrap();
            coll.add_segment_placed(
                &ref_seg.sample_name,
                &ref_seg.contig_name,
                ref_seg.seg_part_no,
                buffer.group_id,
                0, // Reference is always at position 0
                ref_seg.is_rev_comp,
                ref_seg.data.len() as u32,
            )
            .context("Failed to register reference")?;
        }
```

**AFTER:**
```rust
        // Register reference in collection with in_group_id = 0
        #[cfg(feature = "ffi_cost")]
        {
            // For FFI: Cannot access placement_buffer from flush_pack
            // Solution: Add placement directly via collection (individual call acceptable for references)
            let placement = SegmentPlacement {
                sample_name: ref_seg.sample_name.clone(),
                contig_name: ref_seg.contig_name.clone(),
                seg_part_no: ref_seg.seg_part_no as u32,
                group_id: buffer.group_id,
                in_group_id: 0,
                is_rev_comp: ref_seg.is_rev_comp,
                data_size: ref_seg.data.len() as u32,
            };
            let mut coll = collection.lock().unwrap();
            coll.add_segments_placed(&[placement])
                .context("Failed to register reference")?;
        }
        #[cfg(not(feature = "ffi_cost"))]
        {
            let mut coll = collection.lock().unwrap();
            coll.add_segment_placed(
                &ref_seg.sample_name,
                &ref_seg.contig_name,
                ref_seg.seg_part_no,
                buffer.group_id,
                0,
                ref_seg.is_rev_comp,
                ref_seg.data.len() as u32,
            )
            .context("Failed to register reference")?;
        }
```

### Change 9: Update flush_pack() - Delta segment registration (line 927-944)

**BEFORE:**
```rust
    // Register segments in collection with their delta indices
    for (seg, &delta_idx) in buffer.segments.iter().zip(segment_delta_indices.iter()) {
        // in_group_id represents which delta this segment uses
        // 0 = reference, 1+ = delta index (offset by buffer.segments_written and +1 for reference)
        // NOTE: Raw groups (0-15) don't have references, but still use 1+ indexing to match C++ AGC
        let in_group_id = buffer.segments_written + delta_idx + 1;

        let mut coll = collection.lock().unwrap();
        coll.add_segment_placed(
            &seg.sample_name,
            &seg.contig_name,
            seg.seg_part_no,
            buffer.group_id,
            in_group_id,
            seg.is_rev_comp,
            seg.data.len() as u32,
        )
        .context("Failed to register segment")?;
    }
```

**AFTER:**
```rust
    // Register segments in collection with their delta indices
    #[cfg(feature = "ffi_cost")]
    {
        let placements: Vec<SegmentPlacement> = buffer.segments.iter()
            .zip(segment_delta_indices.iter())
            .map(|(seg, &delta_idx)| {
                let in_group_id = buffer.segments_written + delta_idx + 1;
                SegmentPlacement {
                    sample_name: seg.sample_name.clone(),
                    contig_name: seg.contig_name.clone(),
                    seg_part_no: seg.seg_part_no as u32,
                    group_id: buffer.group_id,
                    in_group_id,
                    is_rev_comp: seg.is_rev_comp,
                    data_size: seg.data.len() as u32,
                }
            })
            .collect();

        let mut coll = collection.lock().unwrap();
        coll.add_segments_placed(&placements)
            .context("Failed to register segments")?;
    }
    #[cfg(not(feature = "ffi_cost"))]
    {
        for (seg, &delta_idx) in buffer.segments.iter().zip(segment_delta_indices.iter()) {
            let in_group_id = buffer.segments_written + delta_idx + 1;
            let mut coll = collection.lock().unwrap();
            coll.add_segment_placed(
                &seg.sample_name,
                &seg.contig_name,
                seg.seg_part_no,
                buffer.group_id,
                in_group_id,
                seg.is_rev_comp,
                seg.data.len() as u32,
            )
            .context("Failed to register segment")?;
        }
    }
```

### Change 10: Update finalize() - Metadata writing (lines 690-708)

**BEFORE:**
```rust
        // Write collection metadata to archive
        {
            let mut archive = self.archive.lock().unwrap();
            let mut collection = self.collection.lock().unwrap();

            // Write sample names
            collection
                .store_batch_sample_names(&mut archive)
                .context("Failed to write sample names")?;

            // Write contig names and segment details
            collection
                .store_contig_batch(&mut archive, 0, num_samples)
                .context("Failed to write contig batch")?;

            if self.config.verbosity > 0 {
                eprintln!("Collection metadata written successfully");
            }
        }
```

**AFTER:**
```rust
        // Write collection metadata to archive
        #[cfg(feature = "ffi_cost")]
        {
            // Flush any remaining placements before finalizing
            // (Not needed - we call add_segments_placed immediately in flush_pack)

            let mut collection = self.collection.lock().unwrap();

            // C++ Collection handles sample names internally during add_segments_placed
            // Just need to finalize and write contig batch
            collection.complete_serialization();
            collection.store_contig_batch(0, num_samples);

            if self.config.verbosity > 0 {
                eprintln!("Collection metadata written successfully");
            }
        }
        #[cfg(not(feature = "ffi_cost"))]
        {
            let mut archive = self.archive.lock().unwrap();
            let mut collection = self.collection.lock().unwrap();

            collection
                .store_batch_sample_names(&mut archive)
                .context("Failed to write sample names")?;

            collection
                .store_contig_batch(&mut archive, 0, num_samples)
                .context("Failed to write contig batch")?;

            if self.config.verbosity > 0 {
                eprintln!("Collection metadata written successfully");
            }
        }
```

### Change 11: Update finalize() - Splitters/segment-splitters (lines 710-762)

**AFTER metadata writing:**

```rust
        // Write splitters and segment-splitters metadata (AFTER collection data written)
        #[cfg(feature = "ffi_cost")]
        {
            let mut archive = self.archive.lock().unwrap();

            // Write splitters stream (update existing empty stream)
            if self.config.verbosity > 0 {
                eprintln!("Writing splitters metadata ({} k-mers)...", self.splitters.len());
            }

            let mut splitters_data = Vec::new();
            let mut splitters_vec: Vec<u64> = self.splitters.iter().copied().collect();
            splitters_vec.sort_unstable();
            for kmer in splitters_vec.iter() {
                splitters_data.extend_from_slice(&kmer.to_le_bytes());
            }

            let stream_id = archive.register_stream("splitters")?;
            archive.add_part(stream_id, &splitters_data, splitters_vec.len() as u64)?;

            // Write segment-splitters stream
            let map_segs = self.map_segments.lock().unwrap();
            if self.config.verbosity > 0 {
                eprintln!("Writing segment-splitters metadata ({} groups)...", map_segs.len());
            }

            let mut seg_map_vec: Vec<((u64, u64), u32)> = map_segs
                .iter()
                .map(|(key, &group_id)| ((key.kmer_front, key.kmer_back), group_id))
                .collect();
            seg_map_vec.sort_unstable();

            let mut seg_splitters_data = Vec::new();
            for ((kmer_front, kmer_back), group_id) in seg_map_vec.iter() {
                seg_splitters_data.extend_from_slice(&kmer_front.to_le_bytes());
                seg_splitters_data.extend_from_slice(&kmer_back.to_le_bytes());
                seg_splitters_data.extend_from_slice(&group_id.to_le_bytes());
            }

            let stream_id = archive.register_stream("segment-splitters")?;
            archive.add_part(stream_id, &seg_splitters_data, map_segs.len() as u64)?;
        }
        #[cfg(not(feature = "ffi_cost"))]
        {
            // Original implementation (lines 710-762)
            // ... (keep existing code)
        }
```

## Summary of Changes

**Files modified:**
1. `ragc-core/src/ffi/agc_index.cpp` - Add 2 functions (register_stream, add_part)
2. `ragc-core/src/ffi/agc_index.rs` - Add FFI bindings + wrapper methods
3. `ragc-core/src/streaming_compressor_queue.rs` - Major refactoring with conditional compilation

**Key patterns:**
- All changes guarded by `#[cfg(feature = "ffi_cost")]`
- Rust implementation preserved as fallback
- FFI uses batched `add_segments_placed(&[SegmentPlacement])`
- Archive streams still written via FFI (not through Collection)

**Testing:**
```bash
# Build with FFI
cargo build --release --features ffi_cost

# Test
./target/release/ragc create -o ffi_test.agc -k 21 -s 10000 -m 20 samples/chr5*.fa

# Compare with C++ AGC
/home/erik/agc/bin/agc create -o cpp_test.agc -k 21 -s 10000 -l 20 samples/chr5*.fa
sha256sum ffi_test.agc cpp_test.agc

# Verify segment layout
./target/release/ragc inspect ffi_test.agc --segment-layout > ffi_layout.csv
./target/release/ragc inspect cpp_test.agc --segment-layout > cpp_layout.csv
diff ffi_layout.csv cpp_layout.csv
```
