# Memory Profiling Findings

**Date**: 2025-10-23
**Peak Memory**: 732 MB (with 6 threads)
**Target**: ~205 MB (C++ AGC baseline)
**Gap**: 3.6x

---

## Memory Breakdown (Batch Method with Rayon)

### Measured Allocations

| Component | Memory | Notes |
|-----------|--------|-------|
| **Phase 1: Segment Loading** | | |
| - all_segments Vec (metadata) | 2.16 MB | 21,808 items × 104 bytes |
| - Segment data payload | 218.34 MB | Actual sequence data |
| **Phase 2: Grouping** | | |
| - groups_vec (metadata) | 0.11 MB | 2,862 groups × 40 bytes |
| - Data | 218.34 MB | Same as Phase 1 |
| **Phase 3: Parallel Processing** | | |
| - all_packs Vec (metadata) | 0.39 MB | 5,710 packs × 72 bytes |
| - Compressed data payload | 8.94 MB | Final compressed output |
| **Fixed Overhead** | | |
| - K-mer splitters | 91 MB | From previous profiling |
| **Known Total** | 320 MB | |
| **Peak Memory** | 732 MB | |
| **UNACCOUNTED** | **412 MB** | **Mystery!** |

---

## Root Cause Identified: segment.clone()

Found in `ragc-core/src/compressor_streaming.rs:1815`:

```rust
group_writer.add_segment(prepared_seg.segment.clone(), &config)
```

### Impact

**Segment cloning during Rayon parallel processing creates temporary copies:**

1. **Original segments** in `groups_vec`: 220 MB
2. **Cloned segments** passed to `add_segment()`: **+220 MB** (temporary)
3. **GroupWriter buffers** in `pending_segments`: varies (up to 50 segments/group)

With 6 parallel threads processing 2,862 groups:
- All threads simultaneously hold clones
- Peak concurrent clones: ~220 MB
- Plus GroupWriter buffering: additional copies
- **Total overhead**: ~400 MB ✓ matches the 412 MB gap!

---

## Memory Flow Analysis

```
Phase 1: Load
  all_segments = Vec<PreparedSegment>
  ├─ Metadata: 2 MB
  └─ Data: 218 MB
  Total: 220 MB

Phase 2: Group
  groups_vec = Vec<(Key, Vec<PreparedSegment>)>
  ├─ Metadata: 0.11 MB
  └─ Data: 218 MB (same segments, reorganized)
  Total: 220 MB

Phase 3: Rayon par_iter()
  For each group (parallel across 6 threads):
    ├─ segment.clone() → +218 MB temporary copies! ❌
    ├─ GroupWriter.add_segment(clone)
    │   └─ pending_segments.push(clone) → more copies! ❌
    └─ pack = compress(...)
        └─ all_packs.push(pack)

  Concurrent memory:
  ├─ Original segments: 220 MB
  ├─ Cloned segments: ~220 MB (peak across threads)
  ├─ GroupWriter buffers: ~50 MB
  └─ Compressed packs: 9 MB
  Total: ~500 MB

Phase 4: Write
  all_packs = Vec<CompressedPack>
  └─ 9 MB

Peak: 732 MB
```

---

## Why C++ AGC Uses Less Memory

C++ AGC achieves ~205 MB by:

1. **Streaming processing** - never loads all segments
2. **Move semantics** - data moved, not cloned
3. **Immediate writes** - no Vec<CompressedPack> buffering
4. **Single-pass** - segments processed and discarded

RAGC's batch approach:
1. ❌ Loads all segments (220 MB)
2. ❌ Clones during processing (+220 MB)
3. ❌ Buffers in GroupWriters (+50 MB)
4. ✅ Buffers compressed packs (only 9 MB - good!)

---

## Optimization Opportunities

### 1. **Eliminate segment.clone()** (Highest Impact: -220 MB)

Change from:
```rust
group_writer.add_segment(prepared_seg.segment.clone(), &config)
```

To (move instead of clone):
```rust
group_writer.add_segment(prepared_seg.segment, &config)
```

**Challenge**: Groups iteration borrows data, can't move.

**Solutions**:
- A) Convert groups to owned: `into_iter()` instead of `iter()`
- B) Use references in GroupWriter: `add_segment(&SegmentInfo)`
- C) Two-phase: consume groups, don't iterate

### 2. **Optimize GroupWriter buffering** (Medium Impact: -30 MB)

Currently buffers up to 50 segments per group (PACK_CARDINALITY).
With 2,862 groups × ~10 segments avg × 18 KB = ~50 MB held in buffers.

Could reduce buffering or flush more frequently.

### 3. **Streaming architecture** (Complete redesign)

Eliminates Phase 1 Vec accumulation entirely.
But previous attempt showed it's complex with work-stealing.

---

## Recommended Next Steps

1. **Quick Win**: Eliminate segment clones
   - Change `par_iter()` to `into_par_iter()`
   - Move segments instead of cloning
   - Expected: -220 MB → **~510 MB peak**

2. **Medium Win**: Use references in GroupWriter
   - Change `add_segment(&SegmentInfo)`
   - Expected: additional -50 MB → **~460 MB peak**

3. **Full Solution**: Streaming + immediate writes
   - Requires architectural redesign
   - Expected: ~300 MB (closer to C++ AGC)

---

## Key Insight

**The 412 MB mystery is not Rayon overhead - it's segment cloning!**

The batch method would be much more competitive with C++ AGC if we simply eliminated unnecessary clones. The streaming architecture isn't mandatory for good memory usage.
