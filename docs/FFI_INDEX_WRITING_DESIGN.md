# FFI Integration Design: C++ AGC Index Writing from Rust

## Executive Summary

This document describes the design for replacing RAGC's index writing implementation with C++ AGC's via FFI, while keeping RAGC's segmentation and compression pipeline intact. This hybrid approach allows us to achieve byte-identical archives by delegating the complex metadata serialization to the proven C++ implementation.

**Goal**: Make RAGC produce **byte-for-byte identical archives** to C++ AGC by using C++ AGC's `CCollection_V3` for all index writing operations.

---

## Current Architecture Analysis

### RAGC Index Writing Pipeline

RAGC's index writing happens in `ragc-core/src/worker.rs::create_agc_archive()`:

```
1. Create Archive → Archive::new_writer()
2. Write file_type_info stream (version metadata)
3. Write params stream (k, segment_size, etc.)
4. Create CollectionV3 for metadata
5. Run compression pipeline → compress_samples_streaming_with_archive()
   - Segments are added to collection via collection.add_segment()
6. Serialize collection metadata:
   - collection.store_batch_sample_names()
   - collection.store_contig_batch()
7. Close archive → archive.close()
```

**Key Data Structures (Rust)**:
- `Archive` (`ragc-common/src/archive.rs`): Stream-based binary format
- `CollectionV3` (`ragc-common/src/collection.rs`): Sample/contig/segment metadata
- `SegmentDesc`: `{group_id, in_group_id, is_rev_comp, raw_length}`

### C++ AGC Index Writing

C++ AGC's index writing in `src/core/agc_compressor.cpp` and `src/common/collection_v3.cpp`:

```
1. Create CArchive (output mode)
2. Create CCollection_V3
3. Set archives: collection_desc->set_archives(in_archive, out_archive, ...)
4. Run compression → AddSampleFiles()
   - Segments added via collection_desc->add_segment_placed()
   - Or batched: collection_desc->add_segments_placed()
5. Finalize:
   - collection_desc->complete_serialization()  // Writes sample names
   - store_metadata() → collection_desc->store_contig_batch()  // Writes contig metadata
   - store_file_type_info()
6. Archive closes automatically
```

**Key C++ Types**:
- `CArchive`: Binary archive with stream/part model
- `CCollection_V3`: Metadata manager
- `segment_desc_t`: Same fields as Rust `SegmentDesc`

---

## Divergence Analysis

### Archive Format Compatibility

**✓ Compatible**: Both use identical stream-based binary format:
- Stream registration: `archive.register_stream(name)` (Rust) ↔ `RegisterStream(name)` (C++)
- Part writing: `archive.add_part(stream_id, data, metadata)` (Rust) ↔ `AddPartBuffered(stream_id, data, metadata)` (C++)
- Footer serialization: Identical varint encoding for offsets/sizes

**⚠ Potential Differences**:
1. **Stream naming**: Rust uses `ragc_common::stream_naming`, C++ has its own helpers
2. **Buffering**: C++ uses `AddPartBuffered()` with deferred writes, Rust writes immediately
3. **Metadata order**: Stream registration and part writing order must match exactly

### Collection Serialization

**Known Differences** (source of archive divergence):
1. **Encoding algorithms**: Subtle differences in varint encoding, zigzag encoding
2. **Compression order**: ZSTD compression may differ in context reuse patterns
3. **Batch boundaries**: C++ uses 1MB batches, Rust uses all-at-once serialization
4. **String encoding**: Contig name delta encoding may differ

**Solution**: Replace entire Rust `CollectionV3` serialization with C++ `CCollection_V3`.

---

## FFI Design

### Strategy: Thin C API Layer

**Approach**: Extend existing `ragc-core/src/ffi/cost.cpp` with C++ AGC collection management.

**Why not `cxx` crate?**
- Current FFI uses raw `extern "C"` (proven pattern in this codebase)
- `cxx` requires complex bridging layer for C++ objects
- Manual FFI gives precise control over data ownership and lifetime

**Why not full C++ AGC linking?**
- Don't want to replace RAGC's segmentation/compression (already correct)
- Only need index writing (Collection) functionality
- Minimize FFI surface area for maintainability

### FFI Layer Components

#### 1. C++ Bridge: `CCollection_V3` Wrapper

**File**: `ragc-core/src/ffi/collection.cpp` (new file)

This will expose C++ AGC's `CCollection_V3` via C ABI:

```cpp
// Opaque handles for C API
extern "C" {
    // Lifecycle
    void* agc_collection_create(uint32_t no_threads, size_t batch_size,
                                 uint32_t segment_size, uint32_t kmer_length);
    void agc_collection_destroy(void* collection);

    // Archive binding
    void agc_collection_set_archive(void* collection, void* out_archive);
    void agc_collection_prepare_for_compression(void* collection);

    // Sample/contig registration (matches C++ API)
    int agc_collection_register_sample_contig(void* collection,
                                               const char* sample_name,
                                               const char* contig_name);

    // Segment placement (batched for performance)
    void agc_collection_add_segments_placed(void* collection,
                                             const char** sample_names,
                                             const char** contig_names,
                                             const uint32_t* places,
                                             const uint32_t* group_ids,
                                             const uint32_t* in_group_ids,
                                             const uint8_t* is_rev_comps,
                                             const uint32_t* raw_lengths,
                                             size_t count);

    // Serialization
    void agc_collection_complete_serialization(void* collection);
    void agc_collection_store_contig_batch(void* collection,
                                            uint32_t id_from, uint32_t id_to);
}
```

**Implementation Strategy**:
- Wrap `CCollection_V3` pointer as opaque `void*`
- C++ `shared_ptr<CArchive>` managed internally
- String marshaling: Copy C strings into C++ `std::string`
- Batch API for segment placement (avoid FFI call overhead)

#### 2. C++ Archive Wrapper

**Problem**: C++ `CCollection_V3` expects `shared_ptr<CArchive>`, but Rust owns `Archive`.

**Solutions Considered**:

**Option A: Expose Rust Archive to C++** ❌
- Would need to wrap Rust `Archive` methods as C API
- Complex bidirectional FFI (Rust → C++ → Rust)
- Lifetime management nightmare

**Option B: Parallel C++ Archive** ⚠️
- Open same file in C++ and Rust
- Rust writes compressed data, C++ writes metadata
- Requires careful stream ID synchronization
- Risk: Interleaved writes to same file

**Option C: C++ Archive Only** ✅ **RECOMMENDED**
- Replace Rust `Archive` with thin wrapper around C++ `CArchive`
- All archive operations go through FFI
- Simpler ownership model
- Matches C++ AGC behavior exactly

**C++ Archive Bridge** (`ragc-core/src/ffi/archive.cpp`):

```cpp
extern "C" {
    // Archive lifecycle
    void* agc_archive_create_writer();
    void agc_archive_destroy(void* archive);

    // Operations
    int agc_archive_open(void* archive, const char* path);
    int agc_archive_close(void* archive);
    void agc_archive_flush_buffers(void* archive);

    // Stream management
    int agc_archive_register_stream(void* archive, const char* name);
    int agc_archive_get_stream_id(void* archive, const char* name);

    // Part writing
    int agc_archive_add_part_buffered(void* archive, int stream_id,
                                       const uint8_t* data, size_t data_len,
                                       uint64_t metadata);

    void agc_archive_set_raw_size(void* archive, int stream_id, uint64_t raw_size);
}
```

#### 3. Rust Wrapper Types

**File**: `ragc-core/src/ffi/mod.rs` (extend existing `ragc_ffi` module)

```rust
pub struct CppArchive {
    ptr: *mut std::ffi::c_void,
}

impl CppArchive {
    pub fn new_writer() -> Self { ... }
    pub fn open(&mut self, path: &str) -> Result<()> { ... }
    pub fn register_stream(&mut self, name: &str) -> usize { ... }
    pub fn add_part(&mut self, stream_id: usize, data: &[u8], metadata: u64) -> Result<()> { ... }
    pub fn close(&mut self) -> Result<()> { ... }
}

pub struct CppCollection {
    ptr: *mut std::ffi::c_void,
}

impl CppCollection {
    pub fn new(no_threads: u32, batch_size: usize, segment_size: u32, kmer_length: u32) -> Self { ... }
    pub fn set_archive(&mut self, archive: &mut CppArchive) { ... }
    pub fn register_sample_contig(&mut self, sample_name: &str, contig_name: &str) -> bool { ... }

    // Batched segment placement
    pub fn add_segments_placed(&mut self, segments: &[SegmentPlacement]) { ... }

    pub fn complete_serialization(&mut self) { ... }
    pub fn store_contig_batch(&mut self, id_from: u32, id_to: u32) { ... }
}

pub struct SegmentPlacement {
    pub sample_name: String,
    pub contig_name: String,
    pub place: u32,
    pub group_id: u32,
    pub in_group_id: u32,
    pub is_rev_comp: bool,
    pub raw_length: u32,
}
```

---

## Data Marshaling Strategy

### Segment Metadata Batching

**Problem**: FFI calls are expensive. RAGC processes millions of segments.

**Solution**: Batch segment placement calls.

**Rust Side** (`ragc-core/src/streaming_compressor_queue.rs`):
```rust
// Buffer segment placements
let mut segment_buffer = Vec::new();
const BATCH_SIZE: usize = 1024;

// During compression
segment_buffer.push(SegmentPlacement { ... });

if segment_buffer.len() >= BATCH_SIZE {
    collection.add_segments_placed(&segment_buffer);
    segment_buffer.clear();
}

// Final flush
if !segment_buffer.is_empty() {
    collection.add_segments_placed(&segment_buffer);
}
```

**C++ Bridge** (pre-allocate C arrays):
```cpp
void agc_collection_add_segments_placed(void* coll_ptr,
    const char** sample_names,
    const char** contig_names,
    const uint32_t* places,
    const uint32_t* group_ids,
    const uint32_t* in_group_ids,
    const uint8_t* is_rev_comps,
    const uint32_t* raw_lengths,
    size_t count)
{
    auto* coll = static_cast<CCollection_V3*>(coll_ptr);
    vector<segments_to_place_t> segments;
    segments.reserve(count);

    for (size_t i = 0; i < count; ++i) {
        segments.emplace_back(
            sample_names[i], contig_names[i],
            places[i], group_ids[i], in_group_ids[i],
            is_rev_comps[i] != 0, raw_lengths[i]
        );
    }

    coll->add_segments_placed(segments);
}
```

### String Handling

**Ownership Model**:
- **Rust → C++**: Copy strings via `CString::new()`, pass raw pointers
- **C++ stores**: `std::string` owns copies
- **Lifetime**: C strings freed immediately after C++ copies

**Pattern**:
```rust
let c_sample = CString::new(sample_name).unwrap();
let c_contig = CString::new(contig_name).unwrap();

unsafe {
    agc_collection_register_sample_contig(
        self.ptr,
        c_sample.as_ptr(),
        c_contig.as_ptr()
    );
}
// c_sample/c_contig dropped here
```

### Archive Pointer Sharing

**Problem**: Rust needs to keep `CppArchive` alive while `CppCollection` uses it.

**Solution**: Reference counting with manual lifetime management:

```rust
pub struct CppArchive {
    ptr: *mut c_void,
    // Track if collection holds reference
    collection_bound: bool,
}

impl CppArchive {
    pub fn bind_to_collection(&mut self) {
        self.collection_bound = true;
    }
}

impl Drop for CppArchive {
    fn drop(&mut self) {
        if self.collection_bound {
            panic!("Cannot drop CppArchive while bound to collection!");
        }
        unsafe { agc_archive_destroy(self.ptr); }
    }
}
```

**Usage**:
```rust
let mut archive = CppArchive::new_writer();
archive.open(output_path)?;

let mut collection = CppCollection::new(...);
collection.set_archive(&mut archive);
archive.bind_to_collection();

// ... compression ...

collection.complete_serialization();
collection.drop(); // Releases archive
archive.close()?;
```

---

## Build System Integration

### Cargo Configuration

**File**: `ragc-core/Cargo.toml`

```toml
[features]
default = ["ffi_cost", "ffi_collection"]
ffi_cost = []
ffi_collection = []

[build-dependencies]
cc = "1.0"
```

### Build Script

**File**: `ragc-core/build.rs`

```rust
fn main() {
    // Existing cost FFI
    #[cfg(feature = "ffi_cost")]
    {
        println!("cargo:rerun-if-changed=src/ffi/cost.cpp");
        cc::Build::new()
            .cpp(true)
            .flag_if_supported("-std=c++17")
            .file("src/ffi/cost.cpp")
            .compile("agc_cost");
    }

    // NEW: Collection FFI
    #[cfg(feature = "ffi_collection")]
    {
        println!("cargo:rerun-if-changed=src/ffi/collection.cpp");
        println!("cargo:rerun-if-changed=src/ffi/archive.cpp");

        // Find C++ AGC source directory
        let agc_src = std::env::var("AGC_SRC_DIR")
            .unwrap_or_else(|_| "/home/erik/agc/src".to_string());

        cc::Build::new()
            .cpp(true)
            .flag_if_supported("-std=c++17")
            .include(format!("{}/common", agc_src))
            .include(format!("{}/core", agc_src))
            // FFI bridge files
            .file("src/ffi/collection.cpp")
            .file("src/ffi/archive.cpp")
            // C++ AGC implementation files
            .file(format!("{}/common/collection_v3.cpp", agc_src))
            .file(format!("{}/common/collection.cpp", agc_src))
            .file(format!("{}/common/archive.cpp", agc_src))
            .compile("agc_collection");

        // Link ZSTD (required by C++ AGC)
        println!("cargo:rustc-link-lib=zstd");
    }

    // Link C++ stdlib
    #[cfg(target_os = "linux")]
    println!("cargo:rustc-link-lib=stdc++");
}
```

**Environment Variable**:
```bash
export AGC_SRC_DIR=/home/erik/agc/src
cargo build --release
```

---

## Integration Points in RAGC

### Modified Files

#### 1. `ragc-core/src/worker.rs`

**Function**: `create_agc_archive()`

**Changes**:
```rust
pub fn create_agc_archive(...) -> anyhow::Result<()> {
    // OLD: Create Rust Archive
    // let mut archive = Archive::new_writer();
    // archive.open(output_path)?;

    // NEW: Create C++ Archive
    let mut archive = crate::ffi::CppArchive::new_writer();
    archive.open(output_path)?;

    // OLD: Create Rust CollectionV3
    // let mut collection = CollectionV3::new();
    // collection.set_config(segment_size, kmer_length as u32, None);
    // collection.prepare_for_compression(&mut archive)?;

    // NEW: Create C++ Collection
    let mut collection = crate::ffi::CppCollection::new(
        num_threads as u32,
        1 << 20,  // 1MB batch size
        segment_size,
        kmer_length as u32
    );
    collection.set_archive(&mut archive);
    collection.prepare_for_compression();

    // Write file_type_info (keep Rust version for now)
    // ... unchanged ...

    // Write params
    // ... unchanged ...

    // Compression pipeline (unchanged)
    compress_samples_streaming_with_archive(
        sample_files,
        splitters,
        ...,
        Some(archive),      // Now CppArchive wrapped in Arc<Mutex<>>
        Some(collection),   // Now CppCollection wrapped in Arc<Mutex<>>
    )?;

    // OLD: Serialize with Rust
    // collection.store_batch_sample_names(&mut archive)?;
    // collection.store_contig_batch(&mut archive, 0, num_samples)?;

    // NEW: Serialize with C++
    collection.complete_serialization();
    collection.store_contig_batch(0, num_samples as u32);

    // Close archive
    archive.close()?;

    Ok(())
}
```

#### 2. `ragc-core/src/streaming_compressor_queue.rs`

**Function**: `write_reference_immediately()` and delta writing

**Changes**: Replace `collection.add_segment()` calls with batched `CppCollection::add_segments_placed()`.

**Before**:
```rust
collection.add_segment(
    sample_name,
    contig_name,
    segment_index,
    group_id,
    in_group_id,
    is_rev_comp,
    raw_length
);
```

**After** (with batching):
```rust
// In struct:
segment_buffer: Vec<SegmentPlacement>,

// During compression:
self.segment_buffer.push(SegmentPlacement {
    sample_name: sample_name.clone(),
    contig_name: contig_name.clone(),
    place: segment_index,
    group_id,
    in_group_id,
    is_rev_comp,
    raw_length,
});

if self.segment_buffer.len() >= 1024 {
    collection.add_segments_placed(&self.segment_buffer);
    self.segment_buffer.clear();
}
```

---

## Testing Strategy

### Phase 1: Minimal FFI Test

**Goal**: Verify C++ collection can be created and used from Rust.

**Test**: `ragc-core/tests/test_cpp_collection.rs`
```rust
#[test]
fn test_cpp_collection_creation() {
    let mut archive = CppArchive::new_writer();
    archive.open("/tmp/test_cpp_ffi.agc").unwrap();

    let mut coll = CppCollection::new(1, 1<<20, 10000, 21);
    coll.set_archive(&mut archive);
    coll.prepare_for_compression();

    assert!(coll.register_sample_contig("sample1", "chr1"));

    let segments = vec![
        SegmentPlacement {
            sample_name: "sample1".to_string(),
            contig_name: "chr1".to_string(),
            place: 0,
            group_id: 16,
            in_group_id: 0,
            is_rev_comp: false,
            raw_length: 10021,
        }
    ];

    coll.add_segments_placed(&segments);
    coll.complete_serialization();
    coll.store_contig_batch(0, 1);

    archive.close().unwrap();
}
```

### Phase 2: Byte-Identical Archive Test

**Goal**: Verify RAGC with C++ collection produces identical archives to C++ AGC.

**Test Setup**:
```bash
# Create with RAGC (using C++ collection FFI)
./target/release/ragc create -o ragc_cpp_ffi.agc -k 21 -s 10000 -m 20 minimal_test.fa

# Create with C++ AGC
/home/erik/agc/bin/agc create -o cpp.agc -k 21 -s 10000 -l 20 minimal_test.fa

# Compare
diff ragc_cpp_ffi.agc cpp.agc  # Should be identical
sha256sum ragc_cpp_ffi.agc cpp.agc
```

### Phase 3: Multi-Sample Test

**Dataset**: Yeast chr5 (3 samples)

**Verification**:
1. Size comparison (should match exactly)
2. Extraction verification (all samples byte-identical)
3. C++ AGC can read RAGC archive
4. RAGC can read C++ AGC archive

---

## Implementation Phases

### Phase 1: C++ Bridge Implementation (Week 1)

**Deliverables**:
- [ ] `ragc-core/src/ffi/archive.cpp` - C++ Archive wrapper
- [ ] `ragc-core/src/ffi/collection.cpp` - C++ Collection wrapper
- [ ] Updated `ragc-core/build.rs` - Build integration
- [ ] Rust wrapper types in `ragc-core/src/ffi/mod.rs`

**Acceptance**: FFI test passes (Phase 1 testing)

### Phase 2: Integration into RAGC (Week 2)

**Deliverables**:
- [ ] Modified `ragc-core/src/worker.rs::create_agc_archive()`
- [ ] Batched segment placement in `streaming_compressor_queue.rs`
- [ ] Remove Rust `CollectionV3` serialization code (dead code)

**Acceptance**: RAGC creates valid archives (readable by C++ AGC)

### Phase 3: Byte-Identical Verification (Week 2-3)

**Deliverables**:
- [ ] Minimal test case (single sample, single chromosome)
- [ ] Multi-sample test case (yeast chr5)
- [ ] CI integration (automated byte-identical checks)

**Acceptance**: Archives match C++ AGC byte-for-byte

### Phase 4: Cleanup & Documentation (Week 3)

**Deliverables**:
- [ ] Remove dead Rust collection code
- [ ] Performance benchmarking
- [ ] Documentation updates
- [ ] FFI safety audit

---

## Risk Analysis

### Memory Safety

**Risk**: Dangling pointers, use-after-free in FFI boundary

**Mitigation**:
- RAII wrappers (Drop trait) for C++ objects
- Clear ownership rules documented
- Valgrind testing
- AddressSanitizer builds

### Thread Safety

**Risk**: C++ collection accessed from multiple Rust threads

**Mitigation**:
- Wrap in `Arc<Mutex<CppCollection>>`
- C++ `CCollection_V3` already uses internal mutex
- Single-threaded FFI calls (mutex serializes)

### Build Complexity

**Risk**: C++ AGC source dependency breaks portability

**Mitigation**:
- Feature flag (`ffi_collection`) to make optional
- Fallback to Rust collection (existing code)
- Document required C++ AGC version
- Consider vendoring C++ AGC sources

### String Encoding

**Risk**: UTF-8 (Rust) vs implementation-defined (C++) encoding mismatches

**Mitigation**:
- Both codebases assume UTF-8 for file paths
- Test with non-ASCII sample names
- Validate with international character sets

---

## Performance Considerations

### FFI Call Overhead

**Concern**: Millions of `add_segment()` calls → slow FFI transitions

**Solution**: Batching (1024 segments per FFI call)

**Expected Impact**:
- 1M segments / 1024 = ~1000 FFI calls (negligible overhead)
- String copies: ~100MB/s allocation rate (acceptable)

### Archive Buffering

**C++ AGC uses `AddPartBuffered()`** with deferred writes to minimize I/O.

**RAGC must match** this pattern via FFI to avoid performance regression.

### Memory Usage

**C++ Collection**: Uses 1MB batches (configurable)

**RAGC segment buffer**: 1024 × ~100 bytes = ~100KB (minimal)

**Total overhead**: < 5MB (acceptable)

---

## Alternative Approaches Considered

### Alternative 1: Port C++ Collection to Rust

**Pros**:
- Pure Rust (no FFI)
- Easier to maintain long-term

**Cons**:
- High effort (2000+ lines of complex C++)
- Subtle bugs likely (serialization is tricky)
- Still won't guarantee byte-identical output

**Verdict**: ❌ Rejected (too risky, defeats purpose)

### Alternative 2: Use `cxx` Crate

**Pros**:
- Type-safe FFI bridge
- Automatic C++/Rust type conversions

**Cons**:
- Requires significant boilerplate
- Doesn't support `shared_ptr` transparently
- Existing FFI pattern works well

**Verdict**: ⚠️ Consider for future refactor (not now)

### Alternative 3: Replace Entire Archive System

**Pros**:
- Ultimate compatibility

**Cons**:
- Defeats purpose (RAGC's compression is correct!)
- Massive scope creep
- Lose all RAGC improvements

**Verdict**: ❌ Rejected (wrong direction)

---

## Success Criteria

1. ✅ **Byte-Identical Archives**: RAGC with C++ FFI produces archives matching C++ AGC exactly
2. ✅ **Cross-Compatibility**: C++ AGC can read RAGC archives, RAGC can read C++ AGC archives
3. ✅ **Performance**: No more than 10% slowdown vs pure Rust implementation
4. ✅ **Reliability**: No crashes, memory leaks, or undefined behavior
5. ✅ **Maintainability**: Clear FFI boundary, well-documented ownership model

---

## Future Enhancements

### Post-FFI Optimization

Once byte-identical archives are achieved:

1. **Profile FFI overhead**: Identify hotspots
2. **Optimize batching**: Tune batch size for best performance
3. **Reduce string copies**: Use string interning if needed
4. **Lazy archive writes**: Defer writes until batch boundary

### Long-Term: Pure Rust Port

If FFI proves successful, consider:

1. **Reverse-engineer C++ logic**: Document exact algorithms
2. **Incremental port**: One component at a time
3. **Keep FFI as reference**: Validate Rust implementation against C++
4. **Eventually remove FFI**: Once Rust version proven equivalent

---

## Appendix: FFI Function Reference

### Archive Functions

| C Function | Rust Wrapper | Purpose |
|------------|-------------|---------|
| `agc_archive_create_writer()` | `CppArchive::new_writer()` | Create archive handle |
| `agc_archive_open(path)` | `CppArchive::open(path)` | Open file for writing |
| `agc_archive_register_stream(name)` | `CppArchive::register_stream(name)` | Register stream, get ID |
| `agc_archive_add_part_buffered(...)` | `CppArchive::add_part(...)` | Write compressed data |
| `agc_archive_flush_buffers()` | `CppArchive::flush()` | Force write pending data |
| `agc_archive_close()` | `CppArchive::close()` | Finalize and close |
| `agc_archive_destroy()` | `Drop::drop()` | Free C++ object |

### Collection Functions

| C Function | Rust Wrapper | Purpose |
|------------|-------------|---------|
| `agc_collection_create(...)` | `CppCollection::new(...)` | Create collection |
| `agc_collection_set_archive(...)` | `CppCollection::set_archive(...)` | Bind to archive |
| `agc_collection_register_sample_contig(...)` | `CppCollection::register_sample_contig(...)` | Register sample/contig |
| `agc_collection_add_segments_placed(...)` | `CppCollection::add_segments_placed(...)` | Add segment batch |
| `agc_collection_complete_serialization()` | `CppCollection::complete_serialization()` | Write sample names |
| `agc_collection_store_contig_batch(...)` | `CppCollection::store_contig_batch(...)` | Write contig metadata |
| `agc_collection_destroy()` | `Drop::drop()` | Free C++ object |

---

## References

- RAGC archive implementation: `/home/erik/ragc/ragc-common/src/archive.rs`
- RAGC collection implementation: `/home/erik/ragc/ragc-common/src/collection.rs`
- C++ AGC archive: `/home/erik/agc/src/common/archive.cpp`
- C++ AGC collection: `/home/erik/agc/src/common/collection_v3.cpp`
- Existing FFI: `/home/erik/ragc/ragc-core/src/ffi/cost.cpp`
- Build script: `/home/erik/ragc/ragc-core/build.rs`

---

**Document Version**: 1.0
**Date**: 2025-11-18
**Author**: Claude (Anthropic)
**Status**: Design Review
