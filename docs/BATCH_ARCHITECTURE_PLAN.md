# Batch Architecture Plan: "Write Packs As They Fill"

## Goal
Implement C++ AGC's exact pack-writing architecture in pure Rust to achieve:
- Byte-identical archives to C++ AGC
- Bounded memory usage (like C++ AGC)
- Better performance (5x+ improvement expected)

## Current Problem

RAGC's streaming mode accumulates ALL segments in groups, then writes everything at finalize():
```
Process sample 1 → segments added to groups
Process sample 2 → segments added to groups
Process sample 3 → segments added to groups
finalize() → write ALL groups at once
```

This causes:
- Unbounded group memory growth
- Different ZSTD frame ordering
- ~4% larger archives

## C++ AGC Architecture

C++ AGC writes packs **as they fill** (when group reaches `contigs_in_pack` segments):
```
Process segment → add to group
  if group.segments.len() >= pack_size:
    write_pack(group)
    group.segments.clear()
```

Key parameters:
- `contigs_in_pack` = 50 (default) - segments per pack
- `pack_cardinality` = 50 (default) - samples per metadata batch

## Implementation Plan

### Phase 1: Add Pack Writing Config
**File**: `ragc-core/src/streaming_compressor_queue.rs`

Add to `StreamingQueueConfig`:
```rust
/// Number of segments per pack (matches C++ AGC's contigs_in_pack)
/// When a group reaches this many segments, write a pack immediately
pub pack_size: usize,  // default: 50
```

### Phase 2: Modify SegmentGroupBuffer
**File**: `ragc-core/src/streaming_compressor_queue.rs`

Add method to check and write when full:
```rust
impl SegmentGroupBuffer {
    /// Check if buffer is full and should write a pack
    fn should_write_pack(&self, pack_size: usize) -> bool {
        self.segments.len() >= pack_size
    }

    /// Write current segments as a pack, clear buffer
    /// Returns the segments that were written
    fn take_segments_for_pack(&mut self) -> Vec<SegmentInfo> {
        std::mem::take(&mut self.segments)
    }
}
```

### Phase 3: Modify Segment Addition Logic
**Location**: Worker thread segment processing

Currently:
```rust
// Add segment to group
group.segments.push(segment_info);
// ... continue processing
```

New:
```rust
// Add segment to group
group.segments.push(segment_info);

// Check if pack is full
if group.should_write_pack(config.pack_size) {
    // Write pack immediately (like C++ AGC)
    flush_pack(group, ...)?;
}
```

### Phase 4: Handle Reference Segment Correctly

The first segment in each group becomes the reference. Current logic:
```rust
if !buffer.ref_written && !buffer.segments.is_empty() {
    let ref_seg = buffer.segments.remove(0);  // First becomes reference
    // ... write reference
}
```

This needs to happen on the FIRST pack write, not just at finalize.

### Phase 5: Thread Safety

Current architecture already uses `Arc<Mutex<BTreeMap<Key, SegmentGroupBuffer>>>`.
Pack writes will happen while holding the group lock, which is correct.

For better parallelism, consider:
- Release lock after extracting segments
- Compress outside lock
- Re-acquire to update group state

### Phase 6: Verification Protocol

After implementation:
```bash
# 1. Create with C++ AGC
/home/erik/agc/bin/agc create -o cpp.agc -k 21 -s 10000 -l 20 -t 1 samples/*.fa

# 2. Create with new RAGC
./target/release/ragc create -o ragc.agc -k 21 -s 10000 -m 20 -t 1 samples/*.fa

# 3. Compare (should be identical or very close)
sha256sum cpp.agc ragc.agc
ls -la cpp.agc ragc.agc

# 4. Verify correctness
./target/release/ragc inspect ragc.agc --segment-layout > ragc_layout.csv
./target/release/ragc inspect cpp.agc --segment-layout > cpp_layout.csv
diff ragc_layout.csv cpp_layout.csv
```

### Phase 7: Optional Testing Flag

During development, add a flag to compare both approaches:
```rust
pub struct StreamingQueueConfig {
    // ...
    /// Write packs as they fill (C++ AGC style) vs all at finalize
    pub immediate_pack_writes: bool,  // default: true (new behavior)
}
```

Once verified correct, remove the flag and make immediate writes the only behavior.

## Implementation Order

1. **Add config parameter** - `pack_size: usize` (low risk)
2. **Add pack-full check** - `should_write_pack()` method (low risk)
3. **Implement immediate write** - Call flush_pack when full (medium risk)
4. **Test single-threaded** - Verify with `-t 1` first
5. **Test multi-threaded** - Verify thread safety
6. **Benchmark** - Compare performance
7. **Verify byte-identical** - Compare to C++ AGC output
8. **Remove old code path** - If new is strictly better

## Expected Results

| Metric | Current | After |
|--------|---------|-------|
| Archive size | +4% vs C++ AGC | ~0% (byte-identical) |
| Peak memory | Unbounded groups | Bounded by pack_size |
| Performance | Slower | 5x+ faster (less memory churn) |

## Key Files to Modify

1. `ragc-core/src/streaming_compressor_queue.rs`
   - `StreamingQueueConfig` - add pack_size
   - `SegmentGroupBuffer` - add pack-full logic
   - Worker thread code - call flush_pack when full

2. `ragc-cli/src/main.rs`
   - Add CLI flag for pack_size (optional)

## Risk Assessment

**Low Risk**:
- Adding config parameter
- Adding check methods

**Medium Risk**:
- Changing when flush_pack is called
- Thread safety with immediate writes

**Mitigation**:
- Keep old behavior available via flag during testing
- Test extensively with single-thread first
- Compare archives byte-by-byte
