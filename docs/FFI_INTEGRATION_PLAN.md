# FFI Integration Plan for streaming_compressor_queue.rs

## Overview

Replace Rust `Archive` and `CollectionV3` with C++ FFI versions (`CppArchive` and `CppCollection`) in the streaming compressor to achieve byte-identical archives with C++ AGC.

## Current Architecture

### Key Data Structures

1. **BufferedSegment** (line 115)
   - Stores segment metadata and data
   - Fields: `sample_name`, `contig_name`, `seg_part_no`, `data`, `is_rev_comp`
   - Sorted before writing (matches C++ AGC)

2. **SegmentGroupBuffer** (line 151)
   - Manages packs of 50 segments
   - Tracks reference segment and delta segments
   - Flushes when full

3. **StreamingQueueCompressor** (line 209)
   - Main compressor structure
   - Currently uses: `Arc<Mutex<Archive>>` and `Arc<Mutex<CollectionV3>>`
   - **TARGET**: Replace with `Arc<Mutex<CppArchive>>` and `Arc<Mutex<CppCollection>>`

### Current Rust API Calls

| Location | Rust Call | Purpose |
|----------|-----------|---------|
| Line 271 | `Archive::new_writer()` | Create archive |
| Line 272 | `archive.open(output_path)` | Open archive file |
| Line 275 | `CollectionV3::new()` | Create collection |
| Line 276 | `collection.set_config()` | Configure collection |
| Line 280 | `collection.prepare_for_compression(&mut archive)` | Register collection streams |
| Line 697 | `collection.store_batch_sample_names(&mut archive)` | Write sample metadata |
| Line 702 | `collection.store_contig_batch(&mut archive, 0, num_samples)` | Write contig metadata |
| Line 859, 934, 1011 | `coll.add_segment_placed(...)` | Register segment placement |
| Line 765 | `archive.close()` | Finalize archive |

### Critical Data Flow Points

**Worker threads** (line 1256):
1. Pull ContigTask from queue
2. Split into segments
3. Group by k-mer keys
4. Buffer segments (BufferedSegment)
5. Flush packs when full → calls `flush_pack()` → calls `add_segment_placed()`

**flush_pack()** (line 796):
1. Writes reference segment (line 852: `arch.add_part()`)
2. Writes delta segments (line 961: `arch.add_part()`)
3. Registers segments (lines 859, 934: `coll.add_segment_placed()`)

## FFI API Available

From `/home/erik/ragc/ragc-core/src/ffi/agc_index.rs`:

### CppArchive
```rust
pub struct CppArchive {
    pub fn new_writer() -> Self
    pub fn open(&mut self, path: &str) -> Result<()>
    pub fn close(&mut self) -> Result<()>
    pub fn as_ptr(&mut self) -> *mut c_void  // For passing to Collection
}
```

### CppCollection
```rust
pub struct CppCollection {
    pub fn new() -> Self
    pub fn set_archives(&mut self, archive: &mut CppArchive, num_threads: u32,
                        batch_size: usize, segment_size: u32, kmer_length: u32) -> Result<()>
    pub fn add_segments_placed(&mut self, placements: &[SegmentPlacement]) -> Result<()>
    pub fn complete_serialization(&mut self)
    pub fn store_contig_batch(&mut self, id_from: u32, id_to: u32)
}
```

### SegmentPlacement
```rust
pub struct SegmentPlacement {
    pub sample_name: String,
    pub contig_name: String,
    pub seg_part_no: u32,
    pub group_id: u32,
    pub in_group_id: u32,
    pub is_rev_comp: bool,
    pub data_size: u32,
}
```

## Integration Strategy

### Phase 1: Type Replacements

**In StreamingQueueCompressor struct (line 209):**

```rust
// BEFORE:
archive: Arc<Mutex<Archive>>,
collection: Arc<Mutex<CollectionV3>>,

// AFTER:
#[cfg(feature = "ffi_cost")]
archive: Arc<Mutex<crate::ffi::agc_index::CppArchive>>,
#[cfg(not(feature = "ffi_cost"))]
archive: Arc<Mutex<Archive>>,

#[cfg(feature = "ffi_cost")]
collection: Arc<Mutex<crate::ffi::agc_index::CppCollection>>,
#[cfg(not(feature = "ffi_cost"))]
collection: Arc<Mutex<CollectionV3>>,
```

### Phase 2: Initialization Changes

**In `with_splitters()` constructor (lines 270-342):**

#### Archive Creation (line 271):
```rust
// BEFORE:
let mut archive = Archive::new_writer();
archive.open(output_path)?;

// AFTER:
#[cfg(feature = "ffi_cost")]
let mut archive = {
    let mut arch = crate::ffi::agc_index::CppArchive::new_writer();
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

#### Collection Creation (lines 275-280):
```rust
// BEFORE:
let mut collection = CollectionV3::new();
collection.set_config(config.segment_size as u32, config.k as u32, None);
collection.prepare_for_compression(&mut archive)?;

// AFTER:
#[cfg(feature = "ffi_cost")]
let mut collection = {
    let mut coll = crate::ffi::agc_index::CppCollection::new();
    // C++ AGC set_archives() handles stream registration internally
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

#### Stream Registration (lines 282-339):
**CRITICAL**: C++ AGC Collection handles stream 0-2 internally. We still need:
- file_type_info (stream 3)
- params (stream 4)
- splitters (stream 5)
- segment-splitters (stream 6)

**Problem**: Rust `Archive::register_stream()` and `Archive::add_part()` are NOT in FFI!

**Solution**: Need to add FFI functions:
```rust
// In ffi/agc_index.rs extern "C" block:
fn agc_archive_register_stream(archive_ptr: *mut c_void, name: *const c_char) -> usize;
fn agc_archive_add_part(archive_ptr: *mut c_void, stream_id: usize,
                        data: *const u8, len: usize, metadata: u64) -> bool;
```

```cpp
// In ffi/agc_index.cpp:
size_t agc_archive_register_stream(void* archive_ptr, const char* name) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    return (*archive)->RegisterStream(std::string(name));
}

bool agc_archive_add_part(void* archive_ptr, size_t stream_id,
                          const uint8_t* data, size_t len, uint64_t metadata) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    std::vector<uint8_t> vec(data, data + len);
    return (*archive)->AddPart(stream_id, vec, metadata);
}
```

### Phase 3: Segment Registration Changes

**In flush_pack() (lines 859, 934, 1011):**

Current code calls `add_segment_placed()` for EACH segment individually.

**CRITICAL CHANGE**: Batch segments into `Vec<SegmentPlacement>` and flush periodically!

#### Option A: Immediate Conversion (Simple but less efficient)
```rust
// BEFORE (line 859):
coll.add_segment_placed(
    &ref_seg.sample_name,
    &ref_seg.contig_name,
    ref_seg.seg_part_no,
    buffer.group_id,
    0,
    ref_seg.is_rev_comp,
    ref_seg.data.len() as u32,
)?;

// AFTER:
#[cfg(feature = "ffi_cost")]
{
    let placement = crate::ffi::agc_index::SegmentPlacement {
        sample_name: ref_seg.sample_name.clone(),
        contig_name: ref_seg.contig_name.clone(),
        seg_part_no: ref_seg.seg_part_no as u32,
        group_id: buffer.group_id,
        in_group_id: 0,
        is_rev_comp: ref_seg.is_rev_comp,
        data_size: ref_seg.data.len() as u32,
    };
    coll.add_segments_placed(&[placement])?;
}
#[cfg(not(feature = "ffi_cost"))]
{
    coll.add_segment_placed(
        &ref_seg.sample_name,
        &ref_seg.contig_name,
        ref_seg.seg_part_no,
        buffer.group_id,
        0,
        ref_seg.is_rev_comp,
        ref_seg.data.len() as u32,
    )?;
}
```

#### Option B: Batched (Better performance, matches worker.rs pattern)

**Problem**: `flush_pack()` is called from multiple places. Need thread-local buffers.

**Solution**: Add placement buffer to `StreamingQueueCompressor`:
```rust
// In struct (line 209):
placement_buffers: Arc<Mutex<Vec<crate::ffi::agc_index::SegmentPlacement>>>,
```

**In flush_pack():**
```rust
// Accumulate placements
let placement = SegmentPlacement { ... };
let mut buffers = placement_buffers.lock().unwrap();
buffers.push(placement);

// Flush when buffer reaches 1000 (matches C++ AGC batch_size)
if buffers.len() >= 1000 {
    coll.add_segments_placed(&buffers)?;
    buffers.clear();
}
```

**In finalize() before metadata write:**
```rust
// Flush remaining placements
{
    let mut buffers = self.placement_buffers.lock().unwrap();
    if !buffers.is_empty() {
        let mut coll = self.collection.lock().unwrap();
        coll.add_segments_placed(&buffers)?;
        buffers.clear();
    }
}
```

### Phase 4: Finalization Changes

**In finalize() (lines 697-702):**

```rust
// BEFORE:
collection.store_batch_sample_names(&mut archive)?;
collection.store_contig_batch(&mut archive, 0, num_samples)?;

// AFTER:
#[cfg(feature = "ffi_cost")]
{
    collection.complete_serialization();
    collection.store_contig_batch(0, num_samples);
}
#[cfg(not(feature = "ffi_cost"))]
{
    collection.store_batch_sample_names(&mut archive)?;
    collection.store_contig_batch(&mut archive, 0, num_samples)?;
}
```

**Note**: C++ `store_contig_batch()` calls `complete_serialization()` internally, but calling it explicitly is safer.

## Data Structure Conversions

### BufferedSegment → SegmentPlacement

**Direct mapping:**
```rust
let placement = SegmentPlacement {
    sample_name: buffered.sample_name.clone(),
    contig_name: buffered.contig_name.clone(),
    seg_part_no: buffered.seg_part_no as u32,
    group_id: group_id,           // from SegmentGroupBuffer
    in_group_id: in_group_id,     // calculated during flush
    is_rev_comp: buffered.is_rev_comp,
    data_size: buffered.data.len() as u32,
};
```

**Key differences:**
- `BufferedSegment` has segment **data** (Vec<u8>)
- `SegmentPlacement` has **data_size** only (metadata)
- Data is written to archive separately via `add_part()`

## Implementation Checklist

### Step 1: Add FFI Functions
- [ ] Add `agc_archive_register_stream()` to ffi/agc_index.cpp
- [ ] Add `agc_archive_add_part()` to ffi/agc_index.cpp
- [ ] Add Rust bindings in ffi/agc_index.rs
- [ ] Add wrapper methods to `CppArchive` struct
- [ ] Test compilation

### Step 2: Update StreamingQueueCompressor
- [ ] Add `#[cfg(feature = "ffi_cost")]` conditional types for archive/collection
- [ ] Add placement buffer field
- [ ] Update constructor `with_splitters()`
- [ ] Add helper method `flush_placements()`

### Step 3: Update flush_pack()
- [ ] Replace individual `add_segment_placed()` with batched placements
- [ ] Update reference segment registration (line 859)
- [ ] Update delta segment registration (line 934)
- [ ] Add placement buffer flushing

### Step 4: Update finalize()
- [ ] Flush remaining placements before metadata write
- [ ] Replace Rust collection methods with FFI versions
- [ ] Keep archive stream writes for params/splitters/etc

### Step 5: Testing
- [ ] Compile with `--features ffi_cost`
- [ ] Run minimal test case (chrV yeast)
- [ ] Compare archive SHA256 with C++ AGC
- [ ] Verify segment layout CSV matches exactly

## Potential Issues

### Issue 1: Archive Stream Registration Order

**Problem**: C++ Collection registers streams 0-2 internally. Rust code registers file_type_info, params, splitters, segment-splitters.

**Risk**: Stream IDs may conflict or be out of order.

**Solution**:
- Check C++ AGC stream registration order in collection_v3.cpp
- Ensure Rust follows same order
- Consider registering ALL streams via FFI for consistency

### Issue 2: Placement Batching Performance

**Problem**: Batching adds locking overhead for placement buffer.

**Risk**: May slow down worker threads.

**Solution**:
- Use thread-local buffers (one per worker)
- Flush to global buffer periodically
- Benchmark before/after

### Issue 3: Error Handling

**Problem**: C++ functions return `bool` for errors, Rust uses `Result<T>`.

**Risk**: Silent failures if error checking is incomplete.

**Solution**:
- Always check bool returns and convert to `anyhow::Error`
- Add context to all FFI calls
- Log C++ errors if possible (add error message FFI)

## Alternative: Minimal FFI Surface

Instead of replacing Archive/Collection entirely, keep Rust versions but use FFI for:
1. GroupingEngine (already done)
2. Segment compression (LZ/ZSTD)
3. Archive writing (add_part)

**Pros**: Less invasive, easier to debug
**Cons**: May not achieve byte-identical archives due to subtle differences

## Recommendation

**Start with Phase 1-2** (Archive/Collection FFI) as a feature flag experiment:
- Implement conditional compilation
- Keep Rust version as fallback
- Test on minimal dataset
- Compare archives byte-by-byte
- If successful, proceed to Phase 3-4

**If FFI integration fails:**
- Document exact divergence point
- Consider hybrid approach (Rust archive + C++ collection)
- Focus on algorithmic parity rather than FFI

## Files to Modify

1. `/home/erik/ragc/ragc-core/src/ffi/agc_index.rs` - Add archive stream functions
2. `/home/erik/ragc/ragc-core/src/ffi/agc_index.cpp` - Implement FFI bridge
3. `/home/erik/ragc/ragc-core/src/streaming_compressor_queue.rs` - Main integration
4. `/home/erik/ragc/ragc-core/Cargo.toml` - Ensure ffi_cost feature enabled
5. `/home/erik/ragc/build.rs` - Verify C++ compilation includes agc_index.cpp

## Success Criteria

1. **Compiles** with `cargo build --release --features ffi_cost`
2. **Runs** without crashes on test dataset
3. **Archives match** C++ AGC byte-for-byte (SHA256 identical)
4. **Segment layout CSV** matches exactly
5. **Performance** within 10% of current RAGC (FFI overhead acceptable)
