# FFI Integration Implementation Checklist

## Phase 1: FFI Infrastructure

### File: ragc-core/src/ffi/agc_index.cpp

- [ ] Add `agc_archive_register_stream()` function (after line 73)
  ```cpp
  size_t agc_archive_register_stream(void* archive_ptr, const char* name)
  ```

- [ ] Add `agc_archive_add_part()` function (after line 73)
  ```cpp
  bool agc_archive_add_part(void* archive_ptr, size_t stream_id,
                           const uint8_t* data, size_t len, uint64_t metadata)
  ```

- [ ] Test compilation: `cargo build --release --features ffi_cost`

### File: ragc-core/src/ffi/agc_index.rs

- [ ] Add extern "C" declarations (after line 40)
  - [ ] `fn agc_archive_register_stream(...) -> usize`
  - [ ] `fn agc_archive_add_part(...) -> bool`

- [ ] Add wrapper methods to `CppArchive` impl (after line 97)
  - [ ] `pub fn register_stream(&mut self, name: &str) -> Result<usize>`
  - [ ] `pub fn add_part(&mut self, stream_id: usize, data: &[u8], metadata: u64) -> Result<()>`

- [ ] Test compilation: `cargo build --release --features ffi_cost`

## Phase 2: Streaming Compressor Changes

### File: ragc-core/src/streaming_compressor_queue.rs

#### 2.1 Imports (line 9)

- [ ] Add conditional import:
  ```rust
  #[cfg(feature = "ffi_cost")]
  use crate::ffi::agc_index::{CppArchive, CppCollection, SegmentPlacement};
  ```

#### 2.2 Struct Fields (lines 212-240)

- [ ] Change `archive` field to conditional type
  ```rust
  #[cfg(feature = "ffi_cost")]
  archive: Arc<Mutex<CppArchive>>,
  #[cfg(not(feature = "ffi_cost"))]
  archive: Arc<Mutex<Archive>>,
  ```

- [ ] Change `collection` field to conditional type
  ```rust
  #[cfg(feature = "ffi_cost")]
  collection: Arc<Mutex<CppCollection>>,
  #[cfg(not(feature = "ffi_cost"))]
  collection: Arc<Mutex<CollectionV3>>,
  ```

- [ ] Add placement buffer field (optional optimization)
  ```rust
  #[cfg(feature = "ffi_cost")]
  placement_buffer: Arc<Mutex<Vec<SegmentPlacement>>>,
  ```

#### 2.3 Constructor - Archive Creation (lines 270-272)

- [ ] Wrap archive creation in conditional blocks
  ```rust
  #[cfg(feature = "ffi_cost")]
  let mut archive = {
      let mut arch = CppArchive::new_writer();
      arch.open(output_path.to_str().unwrap())?;
      arch
  };
  #[cfg(not(feature = "ffi_cost"))]
  let mut archive = { /* existing Rust code */ };
  ```

#### 2.4 Constructor - Collection Creation (lines 274-280)

- [ ] Wrap collection creation in conditional blocks
  ```rust
  #[cfg(feature = "ffi_cost")]
  let mut collection = {
      let mut coll = CppCollection::new();
      coll.set_archives(&mut archive, num_threads, 1000, segment_size, k)?;
      coll
  };
  #[cfg(not(feature = "ffi_cost"))]
  let mut collection = { /* existing Rust code */ };
  ```

#### 2.5 Constructor - Stream Registration (lines 282-339)

- [ ] Wrap file_type_info stream in conditional
  ```rust
  #[cfg(feature = "ffi_cost")]
  {
      // Use archive.register_stream() and archive.add_part()
  }
  #[cfg(not(feature = "ffi_cost"))]
  {
      // Keep existing code
  }
  ```

- [ ] Wrap params stream in conditional
- [ ] Wrap splitters stream in conditional
- [ ] Wrap segment-splitters stream in conditional

#### 2.6 Constructor - Return Statement (line 350+)

- [ ] Add placement_buffer to return struct (if using)
  ```rust
  Ok(Self {
      // ... existing fields ...
      #[cfg(feature = "ffi_cost")]
      placement_buffer,
  })
  ```

#### 2.7 Helper Method (optional)

- [ ] Add placement buffer flush method
  ```rust
  #[cfg(feature = "ffi_cost")]
  fn flush_placement_buffer(&self) -> Result<()> {
      // Implementation
  }
  ```

## Phase 3: Segment Registration Changes

### File: ragc-core/src/streaming_compressor_queue.rs (continued)

#### 3.1 flush_pack() - Reference Segment (line 859)

- [ ] Wrap reference registration in conditional
  ```rust
  #[cfg(feature = "ffi_cost")]
  {
      let placement = SegmentPlacement { /* ... */ };
      coll.add_segments_placed(&[placement])?;
  }
  #[cfg(not(feature = "ffi_cost"))]
  {
      coll.add_segment_placed(/* ... */)?;
  }
  ```

#### 3.2 flush_pack() - Delta Segments (lines 927-944)

- [ ] Wrap delta registration in conditional
  ```rust
  #[cfg(feature = "ffi_cost")]
  {
      let placements: Vec<SegmentPlacement> = segments.iter()
          .map(|seg| SegmentPlacement { /* ... */ })
          .collect();
      coll.add_segments_placed(&placements)?;
  }
  #[cfg(not(feature = "ffi_cost"))]
  {
      // Keep existing loop
  }
  ```

## Phase 4: Finalization Changes

### File: ragc-core/src/streaming_compressor_queue.rs (continued)

#### 4.1 finalize() - Metadata Writing (lines 690-708)

- [ ] Wrap metadata writing in conditional
  ```rust
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

#### 4.2 finalize() - Splitters Streams (lines 710-762)

- [ ] Wrap splitters write in conditional
  ```rust
  #[cfg(feature = "ffi_cost")]
  {
      let stream_id = archive.register_stream("splitters")?;
      archive.add_part(stream_id, &splitters_data, len)?;
  }
  #[cfg(not(feature = "ffi_cost"))]
  {
      // Keep existing code
  }
  ```

- [ ] Wrap segment-splitters write in conditional

## Phase 5: Testing

### 5.1 Compilation Tests

- [ ] Clean build: `cargo clean`
- [ ] Build without FFI: `cargo build --release`
  - [ ] Should compile successfully (uses Rust implementations)

- [ ] Build with FFI: `cargo build --release --features ffi_cost`
  - [ ] Should compile successfully (uses C++ FFI)
  - [ ] No linker errors
  - [ ] No undefined symbols

### 5.2 Minimal Test Case

- [ ] Create test archive (single sample):
  ```bash
  ./target/release/ragc create -o test_ffi.agc -k 21 -s 10000 -m 20 \
      /home/erik/ragc/samples/chr5_001.fa
  ```

- [ ] Verify archive opens:
  ```bash
  ./target/release/ragc list test_ffi.agc
  ```

- [ ] Extract and verify:
  ```bash
  ./target/release/ragc extract test_ffi.agc sample_name > extracted.fa
  diff /home/erik/ragc/samples/chr5_001.fa extracted.fa
  ```

### 5.3 Multi-Sample Test

- [ ] Create multi-sample archive:
  ```bash
  ./target/release/ragc create -o multi_ffi.agc -k 21 -s 10000 -m 20 \
      /home/erik/ragc/samples/chr5*.fa
  ```

- [ ] Verify all samples listed:
  ```bash
  ./target/release/ragc list multi_ffi.agc | wc -l
  # Should match number of input files
  ```

### 5.4 Comparison with C++ AGC

- [ ] Create C++ AGC archive:
  ```bash
  /home/erik/agc/bin/agc create -o cpp_test.agc -k 21 -s 10000 -l 20 \
      /home/erik/ragc/samples/chr5*.fa
  ```

- [ ] Compare sizes:
  ```bash
  ls -lh test_ffi.agc cpp_test.agc
  # Should be within 5%
  ```

- [ ] Compare segment layouts:
  ```bash
  ./target/release/ragc inspect test_ffi.agc --segment-layout > ffi.csv
  ./target/release/ragc inspect cpp_test.agc --segment-layout > cpp.csv
  diff ffi.csv cpp.csv
  # GOAL: Should be identical (no output)
  ```

- [ ] Compare SHA256 (stretch goal):
  ```bash
  sha256sum test_ffi.agc cpp_test.agc
  # GOAL: Should match byte-for-byte
  ```

### 5.5 Performance Test

- [ ] Benchmark FFI version:
  ```bash
  /usr/bin/time -v ./target/release/ragc create -o bench_ffi.agc \
      -k 21 -s 10000 -m 20 /home/erik/ragc/samples/chr5*.fa
  ```

- [ ] Benchmark Rust version:
  ```bash
  cargo build --release  # Without ffi_cost feature
  /usr/bin/time -v ./target/release/ragc create -o bench_rust.agc \
      -k 21 -s 10000 -m 20 /home/erik/ragc/samples/chr5*.fa
  ```

- [ ] Compare:
  - [ ] Wall time (should be within 20%)
  - [ ] Memory usage (should be within 10%)
  - [ ] Archive size (should be within 5%)

### 5.6 Stress Test

- [ ] Test with large dataset (if available):
  ```bash
  ./target/release/ragc create -o yeast_ffi.agc -k 21 -s 10000 -m 20 \
      /path/to/yeast235/*.fa
  ```

- [ ] Monitor for:
  - [ ] No segfaults
  - [ ] No memory leaks (use valgrind if needed)
  - [ ] Stable memory usage
  - [ ] Successful completion

## Phase 6: Debugging (If Tests Fail)

### If Archive Won't Open

- [ ] Check archive.close() was called
- [ ] Check all streams registered correctly
- [ ] Verify file exists and has non-zero size
- [ ] Try opening with C++ AGC: `/home/erik/agc/bin/agc getset test_ffi.agc`

### If Segment Layouts Differ

- [ ] Enable verbose logging (verbosity = 2)
- [ ] Compare RAGC_SEGMENT debug output with C++ AGC
- [ ] Check group_id assignments match
- [ ] Check in_group_id assignments match
- [ ] Check is_rev_comp flags match

### If Archives Different Sizes

- [ ] Export segment layouts to CSV
- [ ] Count total segments: `wc -l ffi.csv cpp.csv`
- [ ] Check for missing/extra segments
- [ ] Verify group counts match

### If Segfaults Occur

- [ ] Run with gdb: `gdb --args ./target/release/ragc create ...`
- [ ] Check FFI pointer validity
- [ ] Verify C++ objects not freed prematurely
- [ ] Check thread safety of FFI calls

## Phase 7: Documentation

- [ ] Update CLAUDE.md with results
  - [ ] Success metrics achieved
  - [ ] Performance numbers
  - [ ] Any issues encountered

- [ ] Update README.md
  - [ ] Note ffi_cost feature requirement
  - [ ] Build instructions for FFI variant

- [ ] Git commit with detailed message:
  ```
  feat(ffi): Integrate C++ AGC Archive/Collection for byte-identical archives

  - Replace Rust Archive/CollectionV3 with C++ FFI versions
  - Add register_stream() and add_part() FFI bridges
  - Batch segment registrations using SegmentPlacement
  - Conditional compilation via ffi_cost feature flag

  Results:
  - Archive size: X MB (vs C++ AGC Y MB, Z% difference)
  - Segment layout: [identical/differs by N segments]
  - Performance: X.Xs (vs Rust-only Y.Ys, Z% overhead)
  - SHA256: [matches/differs]
  ```

## Sign-Off Criteria

### Minimum Viable

- [x] Compiles with `--features ffi_cost`
- [x] Creates readable archives
- [x] No crashes or segfaults
- [x] Segments extract correctly

### Target Success

- [ ] Segment layout CSV matches C++ AGC exactly
- [ ] Archive size within 5% of C++ AGC
- [ ] Performance within 20% of Rust-only
- [ ] All samples extract byte-identically

### Stretch Goals

- [ ] SHA256 matches C++ AGC byte-for-byte
- [ ] Performance within 10% of Rust-only
- [ ] Works with all test datasets (chr5, yeast10, yeast235)

---

**Started**: [Date]
**Completed**: [Date]
**Result**: [Success/Partial/Failed]
**Notes**: [Key findings, issues, workarounds]
