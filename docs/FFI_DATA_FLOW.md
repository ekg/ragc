# FFI Integration Data Flow

## Current Rust-only Flow

```
┌─────────────────────────────────────────────────────────────┐
│ StreamingQueueCompressor::with_splitters()                  │
│                                                              │
│  1. Archive::new_writer()                                   │
│     └─> archive.open(path)                                  │
│                                                              │
│  2. CollectionV3::new()                                     │
│     ├─> collection.set_config(segment_size, k)             │
│     └─> collection.prepare_for_compression(&mut archive)   │
│         └─> Registers streams 0-2 in archive               │
│                                                              │
│  3. archive.register_stream("file_type_info") → stream 3   │
│     archive.register_stream("params")         → stream 4   │
│     archive.register_stream("splitters")      → stream 5   │
│     archive.register_stream("segment-splitters") → stream 6│
│                                                              │
│  4. Spawn worker threads                                    │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│ worker_thread() [runs in parallel]                          │
│                                                              │
│  1. queue.pull() → ContigTask                               │
│  2. split_at_splitters_with_size() → Vec<Segment>          │
│  3. For each segment:                                       │
│     ├─> Normalize k-mers (front <= back)                   │
│     ├─> Try to split (via GroupingEngine.find_middle())    │
│     └─> Buffer in SegmentGroupBuffer                       │
│                                                              │
│  4. When buffer reaches PACK_CARDINALITY (50):             │
│     └─> flush_pack()                                        │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│ flush_pack(buffer, collection, archive, ...)                │
│                                                              │
│  1. Sort buffer.segments (alphabetically)                   │
│                                                              │
│  2. Write reference segment:                                │
│     ├─> compress_reference_segment(data)                   │
│     ├─> archive.add_part(ref_stream_id, compressed)        │
│     └─> collection.add_segment_placed(                     │ ◄── REPLACE
│             sample, contig, part_no, group_id, 0, ...)     │     with FFI
│                                                              │
│  3. For each delta segment:                                 │
│     ├─> lz_diff.encode(data)                               │
│     └─> Deduplicate deltas                                 │
│                                                              │
│  4. Pack unique deltas:                                     │
│     ├─> compress_segment_configured(packed_data)           │
│     └─> archive.add_part(stream_id, compressed)            │
│                                                              │
│  5. Register delta segments:                                │
│     └─> For each segment:                                  │
│         └─> collection.add_segment_placed(...)             │ ◄── REPLACE
│                                                              │     with FFI
│  6. Register group in map_segments                          │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│ finalize()                                                   │
│                                                              │
│  1. Close queue, join workers                               │
│  2. Flush all remaining segment groups                      │
│  3. Write metadata:                                         │
│     ├─> collection.store_batch_sample_names(&mut archive)  │ ◄── REPLACE
│     └─> collection.store_contig_batch(&mut archive, ...)   │     with FFI
│                                                              │
│  4. Update splitters streams:                               │
│     ├─> archive.register_stream("splitters")               │
│     ├─> archive.add_part(stream_id, splitters_data)        │
│     ├─> archive.register_stream("segment-splitters")       │
│     └─> archive.add_part(stream_id, seg_map_data)          │
│                                                              │
│  5. archive.close()                                         │
└─────────────────────────────────────────────────────────────┘
```

## New FFI-based Flow

```
┌─────────────────────────────────────────────────────────────┐
│ StreamingQueueCompressor::with_splitters()                  │
│                                                              │
│  1. CppArchive::new_writer()                                │ ◄── C++ FFI
│     └─> archive.open(path)                                  │
│                                                              │
│  2. CppCollection::new()                                    │ ◄── C++ FFI
│     └─> collection.set_archives(                            │
│             &mut archive, num_threads, batch_size,          │
│             segment_size, kmer_length)                      │
│         └─> Registers streams 0-2 internally                │
│         └─> Sets up batch processing (1000 segments/batch)  │
│                                                              │
│  3. archive.register_stream("file_type_info") → stream 3   │ ◄── C++ FFI
│     archive.add_part(stream_id, data, metadata)            │ ◄── C++ FFI
│     archive.register_stream("params")         → stream 4   │
│     archive.add_part(stream_id, params_data, 0)            │
│     archive.register_stream("splitters")      → stream 5   │
│     archive.add_part(stream_id, [], 0)  // Empty initially │
│     archive.register_stream("segment-splitters") → stream 6│
│     archive.add_part(stream_id, [], 0)  // Empty initially │
│                                                              │
│  4. Spawn worker threads                                    │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│ worker_thread() [UNCHANGED - same logic]                    │
│                                                              │
│  1. queue.pull() → ContigTask                               │
│  2. split_at_splitters_with_size() → Vec<Segment>          │
│  3. For each segment:                                       │
│     ├─> Normalize k-mers (front <= back)                   │
│     ├─> GroupingEngine.find_middle() ◄── Already C++ FFI!  │
│     └─> Buffer in SegmentGroupBuffer                       │
│                                                              │
│  4. When buffer reaches PACK_CARDINALITY (50):             │
│     └─> flush_pack()                                        │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│ flush_pack(buffer, collection, archive, ...)                │
│                                                              │
│  1. Sort buffer.segments (alphabetically) [UNCHANGED]       │
│                                                              │
│  2. Write reference segment:                                │
│     ├─> compress_reference_segment(data) [Rust - unchanged]│
│     ├─> archive.add_part(ref_stream_id, compressed)        │ ◄── C++ FFI
│     │                                                        │
│     └─> NEW: Create SegmentPlacement struct                │ ◄── NEW
│         ┌────────────────────────────────────────┐         │
│         │ SegmentPlacement {                     │         │
│         │   sample_name: String,                 │         │
│         │   contig_name: String,                 │         │
│         │   seg_part_no: u32,                    │         │
│         │   group_id: u32,                       │         │
│         │   in_group_id: u32,  // 0 for ref      │         │
│         │   is_rev_comp: bool,                   │         │
│         │   data_size: u32,                      │         │
│         │ }                                       │         │
│         └────────────────────────────────────────┘         │
│         └─> collection.add_segments_placed(&[placement])   │ ◄── C++ FFI
│                                                              │
│  3. For each delta segment: [UNCHANGED]                     │
│     ├─> lz_diff.encode(data)                               │
│     └─> Deduplicate deltas                                 │
│                                                              │
│  4. Pack unique deltas:                                     │
│     ├─> compress_segment_configured(packed_data)           │
│     └─> archive.add_part(stream_id, compressed)            │ ◄── C++ FFI
│                                                              │
│  5. Register delta segments:                                │
│     └─> NEW: Batch convert to Vec<SegmentPlacement>       │ ◄── NEW
│         ┌────────────────────────────────────────┐         │
│         │ placements = segments.iter()           │         │
│         │   .map(|seg| SegmentPlacement {        │         │
│         │     sample_name: seg.sample_name,      │         │
│         │     contig_name: seg.contig_name,      │         │
│         │     seg_part_no: seg.seg_part_no,      │         │
│         │     group_id: buffer.group_id,         │         │
│         │     in_group_id: 1 + delta_idx + ...,  │         │
│         │     is_rev_comp: seg.is_rev_comp,      │         │
│         │     data_size: seg.data.len(),         │         │
│         │   })                                    │         │
│         │   .collect()                            │         │
│         └────────────────────────────────────────┘         │
│         └─> collection.add_segments_placed(&placements)    │ ◄── C++ FFI
│                                                              │
│  6. Register group in map_segments [UNCHANGED]              │
└─────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│ finalize()                                                   │
│                                                              │
│  1. Close queue, join workers [UNCHANGED]                   │
│  2. Flush all remaining segment groups [UNCHANGED]          │
│                                                              │
│  3. Write metadata:                                         │
│     ├─> collection.complete_serialization()                │ ◄── C++ FFI
│     │   └─> Writes sample names internally                  │
│     │   └─> Writes collection-samples stream (0)            │
│     │                                                        │
│     └─> collection.store_contig_batch(0, num_samples)      │ ◄── C++ FFI
│         └─> Writes collection-contigs stream (1)            │
│         └─> Writes collection-details stream (2)            │
│                                                              │
│  4. Update splitters streams:                               │
│     ├─> archive.register_stream("splitters")               │ ◄── C++ FFI
│     ├─> archive.add_part(stream_id, splitters_data, len)   │ ◄── C++ FFI
│     ├─> archive.register_stream("segment-splitters")       │ ◄── C++ FFI
│     └─> archive.add_part(stream_id, seg_map_data, len)     │ ◄── C++ FFI
│                                                              │
│  5. archive.close()                                         │ ◄── C++ FFI
└─────────────────────────────────────────────────────────────┘
```

## Key Differences

| Aspect | Rust-only | FFI |
|--------|-----------|-----|
| **Archive type** | `ragc_common::Archive` | `CppArchive` (C++ AGC) |
| **Collection type** | `ragc_common::CollectionV3` | `CppCollection` (C++ AGC) |
| **Stream registration** | Rust: `collection.prepare_for_compression()` | C++: `collection.set_archives()` |
| **Segment registration** | One-at-a-time: `add_segment_placed(...)` | Batched: `add_segments_placed(&[...])` |
| **Metadata write** | `store_batch_sample_names()` + `store_contig_batch()` | `complete_serialization()` + `store_contig_batch()` |
| **Archive writes** | Rust `add_part()` | C++ FFI `add_part()` |

## Data Structure Conversions

### BufferedSegment → SegmentPlacement

**Purpose**: Convert from Rust internal representation (with data) to C++ FFI metadata-only struct.

```rust
// BufferedSegment (Rust internal - has segment data)
struct BufferedSegment {
    sample_name: String,
    contig_name: String,
    seg_part_no: usize,
    data: Vec<u8>,          // ◄── Segment sequence data
    is_rev_comp: bool,
}

// SegmentPlacement (FFI - metadata only)
struct SegmentPlacement {
    sample_name: String,
    contig_name: String,
    seg_part_no: u32,
    group_id: u32,          // ◄── Added during flush
    in_group_id: u32,       // ◄── Position within group
    is_rev_comp: bool,
    data_size: u32,         // ◄── Only size, not data!
}

// Conversion:
let placement = SegmentPlacement {
    sample_name: buffered.sample_name.clone(),
    contig_name: buffered.contig_name.clone(),
    seg_part_no: buffered.seg_part_no as u32,
    group_id: group_id,                      // From SegmentGroupBuffer
    in_group_id: computed_in_group_id,       // 0 for ref, 1+ for deltas
    is_rev_comp: buffered.is_rev_comp,
    data_size: buffered.data.len() as u32,   // Size only
};
```

**Key insight**: Segment **data** is written to archive via `add_part()`, **metadata** is registered via `add_segments_placed()`.

## FFI Call Frequencies

### High-frequency calls (per segment pack):
- `archive.add_part()` - 2x per pack (reference + deltas)
- `collection.add_segments_placed()` - 2x per pack (reference + batch of deltas)

### Low-frequency calls (once per archive):
- `archive.register_stream()` - 4x (file_type_info, params, splitters, segment-splitters)
- `archive.add_part()` for metadata streams - 4x
- `collection.complete_serialization()` - 1x
- `collection.store_contig_batch()` - 1x
- `archive.close()` - 1x

**Performance impact**: Minimal - FFI overhead only for 2-4 calls per pack (50 segments).

## Thread Safety

**All FFI types wrapped in `Arc<Mutex<T>>`:**
- `Arc<Mutex<CppArchive>>` - Ensures sequential archive writes
- `Arc<Mutex<CppCollection>>` - Ensures sequential collection updates

**Worker threads:**
- Each worker buffers segments independently
- Lock collection/archive ONLY during `flush_pack()`
- Matches C++ AGC's mutex-protected `CArchive` and `CCollection_V3`

## Error Handling

**C++ functions return `bool`:**
```cpp
bool agc_archive_open(void* archive_ptr, const char* path);
bool agc_archive_close(void* archive_ptr);
bool agc_archive_add_part(...);
```

**Rust wrappers convert to `Result<T>`:**
```rust
pub fn open(&mut self, path: &str) -> anyhow::Result<()> {
    unsafe {
        if agc_archive_open(self.ptr, c_path.as_ptr()) {
            Ok(())
        } else {
            Err(anyhow::anyhow!("Failed to open archive: {}", path))
        }
    }
}
```

**Every FFI call includes error context:**
```rust
archive.add_part(stream_id, &data, metadata)
    .context("Failed to write reference segment to archive")?;
```

## Memory Management

**C++ objects wrapped in smart pointers:**
```cpp
void* agc_archive_create_writer() {
    auto* archive = new std::shared_ptr<CArchive>(std::make_shared<CArchive>(false));
    return archive;
}

void agc_archive_destroy(void* archive_ptr) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    delete archive;  // Decrements shared_ptr refcount
}
```

**Rust Drop ensures cleanup:**
```rust
impl Drop for CppArchive {
    fn drop(&mut self) {
        unsafe {
            agc_archive_destroy(self.ptr);
        }
    }
}
```

**No manual memory management needed** - RAII handles everything.

## Testing Strategy

### 1. Unit test FFI bindings
```rust
#[test]
fn test_cpp_archive_roundtrip() {
    let mut archive = CppArchive::new_writer();
    archive.open("/tmp/test.agc").unwrap();
    let stream_id = archive.register_stream("test").unwrap();
    archive.add_part(stream_id, b"data", 4).unwrap();
    archive.close().unwrap();
}
```

### 2. Integration test with small dataset
```bash
cargo build --release --features ffi_cost
./target/release/ragc create -o test_ffi.agc samples/chr5_001.fa
sha256sum test_ffi.agc
```

### 3. Compare with C++ AGC
```bash
/home/erik/agc/bin/agc create -o test_cpp.agc samples/chr5_001.fa
diff <(xxd test_ffi.agc) <(xxd test_cpp.agc)
```

### 4. Segment layout verification
```bash
./target/release/ragc inspect test_ffi.agc --segment-layout > ffi.csv
./target/release/ragc inspect test_cpp.agc --segment-layout > cpp.csv
diff ffi.csv cpp.csv  # MUST be identical
```

### 5. Multi-sample test
```bash
./target/release/ragc create -o yeast_ffi.agc samples/chr5*.fa
/home/erik/agc/bin/agc create -o yeast_cpp.agc samples/chr5*.fa
sha256sum yeast_ffi.agc yeast_cpp.agc  # Should match
```

## Success Metrics

1. **Compiles** without errors with `--features ffi_cost`
2. **No segfaults** during execution
3. **Archives produced** and readable
4. **SHA256 matches** C++ AGC byte-for-byte
5. **Segment layouts identical** line-by-line
6. **Performance** within 20% of Rust-only implementation
