# RAGC vs C++ AGC Segment Ordering Comparison

**Date**: 2025-12-04
**Purpose**: Compare RAGC and C++ AGC segment ordering to identify why `in_group_id` assignments differ

---

## Summary of Divergence

**Test case**: 5-sample yeast dataset
**First divergence**: Position 1202 (BBT#0#chrI seg 4)
- **Group ID**: 17 (SAME in both)
- **in_group_id**: RAGC=5, C++=1 (DIFFERENT)

**Interpretation**: Same segment goes to same group, but is written at different positions within the group.

---

## C++ AGC Flow (from CPP_AGC_SEGMENT_ORDERING.md)

### Phase 1: Segment Collection
- NEW segments → `s_seg_part` (set, automatically sorted)
- KNOWN segments → `vl_seg_part[group_id]` (vector, unsorted)

### Phase 2: Batch Flush (`register_segments()`)
1. **`sort_known(nt)`** (line 974): Sorts ALL segments within EACH group
   - Multi-threaded sorting
   - Applies `std::sort()` to each `vl_seg_part[group_id].l_seg_part`
   - Result: Segments ordered by (sample_name, contig_name, seg_part_no)

2. **`process_new()`** (line 976): Assigns group IDs to new segments
   - Iterates `s_seg_part` (already sorted by set nature)
   - Creates groups for unique k-mer pairs
   - Moves segments to `vl_seg_part[group_id]` **in sorted order**

3. **`distribute_segments()`** (line 986): Distributes segments for raw groups

### Phase 3: Archive Writing (`store_segments()`)
- Retrieves segments from `vl_seg_part[group_id]` via `get_part()`
- Segments are in SORTED order (from `sort_known()`)
- `in_group_id` assigned sequentially as written

**Key invariant**: Segments within each group are written in lexicographic order by (sample_name, contig_name, seg_part_no).

---

## RAGC Flow

### Phase 1: Segment Collection (`flush_batch()`)
```rust
// Line 2513-2533: Get pending segments and sort GLOBALLY
let mut pending = pending_batch_segments.lock().unwrap();
pending.sort();  // Sorts ALL pending segments across ALL groups
```

**CRITICAL**: Pending segments are sorted **globally** across all groups, not per-group.

### Phase 2: Group Assignment and Buffering
```rust
// Lines 2542-2670: Process sorted pending segments
for pend in pending.iter() {
    // Assign group_id
    let group_id = /* lookup or create group */;

    // Get or create buffer for this group
    let buffer = groups_map.entry(pend.key.clone()).or_insert_with(|| ...);

    // Add to buffer
    buffer.segments.push(buffered);

    // Flush if buffer full
    if buffer.should_flush_pack(config.pack_size) {
        flush_pack(buffer, ...)?;
    }
}
```

**Key observations**:
1. Segments are added to `buffer.segments` in **global sorted order**
2. `flush_pack()` is called when buffer reaches `pack_size`
3. Each group has its own buffer that fills independently

### Phase 3: Pack Writing (`flush_pack()`)
```rust
// Line 1302: Sort segments within group BEFORE writing
buffer.segments.sort();

// Lines 1307-1381: Write reference segment (if LZ group)
if use_lz_encoding && !buffer.ref_written {
    let ref_seg = buffer.segments.remove(0);  // First after sort
    // Write reference...
    buffer.ref_written = true;
}

// Lines 1384-onwards: Write remaining segments
for (idx, seg) in buffer.segments.iter().enumerate() {
    let in_group_id = buffer.next_in_group_id;
    buffer.next_in_group_id += 1;
    // Compress and write segment...
    // Register with collection...
}
```

**Key observations**:
1. Segments are sorted **per-pack** before writing
2. `in_group_id` is assigned sequentially within each pack
3. Reference segment is the first segment (after sorting) in the first pack

---

## Potential Divergence Causes

### 1. Pack Boundaries vs Batch Boundaries

**C++ AGC**:
- Sorts ALL segments in a group ONCE (in `sort_known()`)
- Writes all segments from that batch in one go
- Multiple batches → multiple calls to `sort_known()` and `store_segments()`

**RAGC**:
- Sorts segments PER-PACK (in `flush_pack()`)
- A group's segments may span multiple packs
- Each pack is sorted independently

**Example causing divergence**:

Group 17 receives segments in this batch order:
1. AAA#0/chrI/seg0
2. ZZZ#0/chrII/seg0
3. BBT#0/chrI/seg4
4. AAA#0/chrI/seg1

**If pack_size = 2**:

RAGC:
- Pack 1: [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0] → sorted → [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0]
  - AAA#0/chrI/seg0 → in_group_id = 0 (reference)
  - ZZZ#0/chrII/seg0 → in_group_id = 1
- Pack 2: [BBT#0/chrI/seg4, AAA#0/chrI/seg1] → sorted → [AAA#0/chrI/seg1, BBT#0/chrI/seg4]
  - AAA#0/chrI/seg1 → in_group_id = 2
  - BBT#0/chrI/seg4 → in_group_id = 3  ❌ WRONG!

C++ AGC (sorts all before writing):
- [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0, BBT#0/chrI/seg4, AAA#0/chrI/seg1] → sorted
- → [AAA#0/chrI/seg0, AAA#0/chrI/seg1, BBT#0/chrI/seg4, ZZZ#0/chrII/seg0]
  - AAA#0/chrI/seg0 → in_group_id = 0 (reference)
  - AAA#0/chrI/seg1 → in_group_id = 1
  - BBT#0/chrI/seg4 → in_group_id = 2  ✅ CORRECT
  - ZZZ#0/chrII/seg0 → in_group_id = 3

**This is the root cause!**

### 2. Sorting Granularity

**C++ AGC**: Sorts at BATCH granularity (all segments in batch processed together)

**RAGC**: Sorts at PACK granularity (segments in each pack sorted independently)

**Result**: Different `in_group_id` assignments when segments arrive across multiple packs.

---

## Root Cause: Pack-Level Sorting

**Problem**: RAGC sorts segments within each pack, but C++ AGC sorts all segments in a batch before writing.

**Effect**: When a group's segments span multiple packs, RAGC assigns `in_group_id` in pack-by-pack sorted order, which differs from the global sorted order.

**Example from test**:
- BBT#0/chrI seg 4 gets in_group_id=5 in RAGC (position within its pack)
- BBT#0/chrI seg 4 gets in_group_id=1 in C++ AGC (position within entire batch)

---

## Solution: Defer `in_group_id` Assignment Until Final Flush

### Option A: Buffer All Segments Per-Group
1. Collect segments in buffers without assigning `in_group_id`
2. At final flush (end of archive), sort ALL segments in each group
3. Write in sorted order, assigning `in_group_id` sequentially

**Pros**: Matches C++ AGC exactly
**Cons**: Requires buffering all segments in memory (defeats streaming purpose)

### Option B: Global Per-Group Counters with Post-Processing
1. Write segments to packs as currently done
2. Track actual write order in a separate structure
3. After all writing complete, re-number `in_group_id` based on global sort
4. Update collection descriptor with corrected `in_group_id` values

**Pros**: Maintains streaming, corrects afterward
**Cons**: Complex, requires post-processing step

### Option C: Accept Pack-Level Sorting (NO FIX - DOCUMENT DIFFERENCE)
1. Document that RAGC uses pack-level sorting
2. Accept that `in_group_id` assignments differ from C++ AGC
3. Archives are compatible but not byte-identical

**Pros**: No code changes, works correctly
**Cons**: Doesn't achieve byte-identical goal

### Option D: Match C++ AGC's Batch Processing (RECOMMENDED)
1. Remove pack-based flushing from `flush_batch()`
2. Buffer ALL segments from a batch in their respective groups
3. At end of batch, sort each group's segments
4. Write all groups' segments to archive
5. THEN emit pack boundaries

**Pros**: Matches C++ AGC architecture, achieves byte-identical archives
**Cons**: Requires refactoring flush logic

---

## Next Step: Investigate Pack vs Batch Boundaries

**Questions to answer**:
1. When does RAGC call `flush_batch()`? (At pack boundaries? At file boundaries?)
2. When does C++ AGC call `register_segments()`? (Same timing?)
3. Are pack_cardinality and batch boundaries aligned?

**Test to run**:
```bash
# Add logging to see when flush_batch and flush_pack are called
RAGC_TRACE_FLUSH=1 ./target/release/ragc create -o test.agc ...

# Compare with C++ AGC batch logging
AGC_DEBUG_BATCH=1 /home/erik/agc/bin/agc create -o test.agc ...
```

**Expected finding**: RAGC flushes packs mid-batch, C++ AGC flushes entire batches atomically.
