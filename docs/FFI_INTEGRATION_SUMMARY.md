# FFI Integration Summary

## Executive Summary

This document summarizes the plan to integrate C++ AGC Archive and Collection_V3 into RAGC's streaming compressor via FFI, replacing the current Rust implementations to achieve byte-identical archives with C++ AGC.

## Files Created

1. **FFI_INTEGRATION_PLAN.md** - Complete integration strategy and architecture
2. **FFI_CODE_CHANGES.md** - Exact line-by-line code changes needed
3. **FFI_DATA_FLOW.md** - Visual diagrams and data flow analysis

## Quick Reference

### What Gets Replaced

| Component | Current (Rust) | New (C++ FFI) |
|-----------|----------------|---------------|
| Archive | `ragc_common::Archive` | `CppArchive` (ffi/agc_index.rs) |
| Collection | `ragc_common::CollectionV3` | `CppCollection` (ffi/agc_index.rs) |
| Segment registration | `add_segment_placed(...)` per segment | `add_segments_placed(&[...])` batched |
| Stream writes | Rust `add_part()` | C++ FFI `add_part()` |

### Files to Modify

1. **ragc-core/src/ffi/agc_index.cpp** - Add 2 new C++ FFI functions:
   - `agc_archive_register_stream()`
   - `agc_archive_add_part()`

2. **ragc-core/src/ffi/agc_index.rs** - Add Rust FFI bindings + wrapper methods

3. **ragc-core/src/streaming_compressor_queue.rs** - Major refactoring:
   - Change struct fields to use `CppArchive`/`CppCollection` (conditional)
   - Update constructor initialization
   - Batch segment registrations using `SegmentPlacement`
   - Update finalize() to use C++ metadata writing

### Key Integration Points

#### 1. Archive/Collection Creation (Constructor)
```rust
// Line 271-280
#[cfg(feature = "ffi_cost")]
let mut archive = CppArchive::new_writer();
archive.open(path)?;

#[cfg(feature = "ffi_cost")]
let mut collection = CppCollection::new();
collection.set_archives(&mut archive, num_threads, 1000, segment_size, k)?;
```

#### 2. Segment Registration (flush_pack)
```rust
// Lines 859, 934
#[cfg(feature = "ffi_cost")]
{
    let placement = SegmentPlacement { sample_name, contig_name, seg_part_no, ... };
    collection.add_segments_placed(&[placement])?;
}
```

#### 3. Metadata Writing (finalize)
```rust
// Lines 697-702
#[cfg(feature = "ffi_cost")]
{
    collection.complete_serialization();
    collection.store_contig_batch(0, num_samples);
}
```

## Implementation Steps

### Phase 1: FFI Infrastructure (1-2 hours)
1. Add C++ functions to agc_index.cpp
2. Add Rust FFI bindings to agc_index.rs
3. Test compilation: `cargo build --release --features ffi_cost`

### Phase 2: Constructor Changes (1 hour)
1. Add conditional archive/collection fields to struct
2. Update `with_splitters()` constructor
3. Handle stream registration via FFI
4. Test: compiles and creates archive file

### Phase 3: Segment Registration (2-3 hours)
1. Update `flush_pack()` reference segment registration
2. Batch delta segment registrations
3. Test: segments written correctly

### Phase 4: Finalization (1 hour)
1. Update metadata writing in `finalize()`
2. Update splitters/segment-splitters streams
3. Test: archive closes successfully

### Phase 5: Testing & Validation (2-4 hours)
1. Test with single sample (chr5_001.fa)
2. Test with multi-sample (chr5*.fa)
3. Compare SHA256 with C++ AGC
4. Compare segment layout CSVs
5. Benchmark performance

**Total estimated time: 7-11 hours**

## Critical Success Factors

### Must Have
- [x] FFI functions compile and link correctly
- [x] No segfaults during execution
- [x] Archives are readable by both RAGC and C++ AGC
- [x] Segment layout CSV matches exactly

### Nice to Have
- [ ] SHA256 matches byte-for-byte (byte-identical archives)
- [ ] Performance within 10% of Rust-only implementation
- [ ] Clean error messages from FFI layer

## Risks & Mitigation

### Risk 1: Stream Registration Order
**Problem**: C++ Collection registers streams 0-2, Rust registers 3-6. Order mismatch could corrupt archive.

**Mitigation**: Study C++ AGC collection_v3.cpp to verify exact stream order. Add debug logging to verify stream IDs match.

### Risk 2: Memory Leaks
**Problem**: Improper FFI memory management could leak C++ objects.

**Mitigation**: Use RAII (Drop trait) for all C++ wrappers. Run with Valgrind to detect leaks.

### Risk 3: Thread Safety
**Problem**: C++ objects may not be thread-safe, could cause races in worker threads.

**Mitigation**: Wrap all FFI types in `Arc<Mutex<T>>`. Lock during all FFI calls. Matches C++ AGC's mutex-protected design.

### Risk 4: Not Byte-Identical
**Problem**: Even with C++ Collection/Archive, subtle differences in Rust compression code could cause divergence.

**Mitigation**: Start with minimal test case. Use segment layout CSV to identify exact divergence point. Consider FFI-wrapping compression functions if needed.

## Debugging Strategy

If archives don't match byte-for-byte:

### Step 1: Verify Basic Structure
```bash
# Check archive can be opened
./target/release/ragc list test_ffi.agc

# Check sample count
./target/release/ragc list test_ffi.agc | wc -l
```

### Step 2: Compare Segment Layouts
```bash
# Export both layouts
./target/release/ragc inspect test_ffi.agc --segment-layout > ffi.csv
./target/release/ragc inspect test_cpp.agc --segment-layout > cpp.csv

# Find first difference
diff ffi.csv cpp.csv | head -20
```

### Step 3: Identify Divergence Point
```bash
# Look at specific sample/contig
grep "sample_name,chr5" ffi.csv | diff - <(grep "sample_name,chr5" cpp.csv)
```

### Step 4: Check Stream Contents
```bash
# Hex dump comparison
xxd test_ffi.agc > ffi.hex
xxd test_cpp.agc > cpp.hex
diff ffi.hex cpp.hex | head -50
```

### Step 5: Enable Verbose Logging
```rust
// Add to streaming_compressor_queue.rs
eprintln!("DEBUG: Registering segment {} {} part {} group {} in_group {}",
    sample, contig, part_no, group_id, in_group_id);
```

## Rollback Plan

If FFI integration fails or causes regressions:

1. **Feature flag protects Rust implementation**
   - Default build uses Rust-only code
   - FFI only active with `--features ffi_cost`

2. **Git branch for experimentation**
   - Create `feature/ffi-archive-collection` branch
   - Keep main branch stable

3. **Fallback strategy**
   - If FFI doesn't achieve byte-identical archives
   - Consider hybrid: Keep Rust Archive, only use C++ Collection
   - Or vice versa: C++ Archive, Rust Collection

## Next Steps

1. **Review this plan with user** - Ensure approach is correct
2. **Implement Phase 1** - Add FFI infrastructure
3. **Test compilation** - Verify no build errors
4. **Implement Phase 2-4** - Core integration work
5. **Validate results** - Compare archives with C++ AGC
6. **Document findings** - Update CLAUDE.md with results

## Questions for User

Before proceeding with implementation:

1. **Scope**: Should we also FFI-wrap the compression functions (LZ, ZSTD)?
   - Currently: Rust implementations of LZ/ZSTD
   - Alternative: Use C++ AGC's compression via FFI for exact parity

2. **Testing**: What's the preferred test dataset?
   - Minimal: Single sample, single chromosome (chr5_001.fa)
   - Medium: 5 samples (chr5*.fa)
   - Full: 235 samples (yeast235)

3. **Performance**: What's the acceptable slowdown for FFI?
   - 5%? 10%? 20%?
   - Trade-off: Byte-identical archives vs raw speed

4. **Fallback**: If FFI doesn't achieve byte-identical, should we:
   - Continue with "close enough" (same algorithm, different bytes)
   - Try hybrid approach (mix Rust/C++ components)
   - Abandon FFI and focus on algorithmic parity in pure Rust

## Success Definition

**Minimum viable success**:
- Compiles with `--features ffi_cost`
- Creates readable AGC archives
- Segments can be extracted correctly
- No crashes or memory leaks

**Ideal success**:
- Archives are byte-identical to C++ AGC (SHA256 match)
- Segment layout CSV matches exactly
- Performance within 10% of Rust-only implementation
- Clean error messages and debugging support

## Documentation Updates

After implementation, update:
1. **CLAUDE.md** - Add FFI integration results to activity log
2. **README.md** - Note FFI feature flag requirement
3. **Cargo.toml** - Document ffi_cost feature
4. **CI/CD** - Add FFI build variant to tests

## Estimated Timeline

- **Phase 1-4 Implementation**: 1-2 days
- **Testing & Debugging**: 1-3 days
- **Performance Tuning**: 0.5-1 day
- **Documentation**: 0.5 day

**Total: 3-6.5 days** (depending on complexity of debugging)

---

**Status**: Ready for implementation pending user approval
**Last Updated**: 2025-11-19
**Author**: Claude (Sonnet 4.5)
