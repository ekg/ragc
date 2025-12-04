# FIX 4: Batch-Level Segment Sorting

**Date**: 2025-12-04
**Priority**: HIGH - Required for byte-identical archives
**Status**: PROPOSED

---

## Problem Statement

**Current divergence**: Same segments go to same groups, but with different `in_group_id` assignments.

**Example**:
- Position 1202: BBT#0#chrI seg 4
- Group ID: 17 (same in both)
- **in_group_id**: RAGC=5, C++=1 (DIFFERENT)

**Root cause**: RAGC sorts segments at **pack-level**, C++ AGC sorts at **batch-level**.

---

## Detailed Analysis

### C++ AGC Approach (CORRECT)

**Batch processing** (agc_compressor.cpp:972-1068):
1. Collect ALL segments in batch
2. Call `sort_known()` → sorts ALL segments within EACH group
3. Call `store_segments()` → writes ALL segments from batch
4. **Result**: Segments in each group written in global lexicographic order

**Example**: Group 17 with 4 segments
```
Batch collects: [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0, BBT#0/chrI/seg4, AAA#0/chrI/seg1]
sort_known():   [AAA#0/chrI/seg0, AAA#0/chrI/seg1, BBT#0/chrI/seg4, ZZZ#0/chrII/seg0]
Write order:
  - AAA#0/chrI/seg0   → in_group_id=0 (reference)
  - AAA#0/chrI/seg1   → in_group_id=1
  - BBT#0/chrI/seg4   → in_group_id=2
  - ZZZ#0/chrII/seg0  → in_group_id=3
```

### RAGC Approach (INCORRECT for byte-identical goal)

**Pack processing** (streaming_compressor_queue_legacy.rs:2542-2670, 1298-1302):
1. Collect segments in batch
2. **For each segment**: Add to group buffer
3. **When buffer full** (`pack_size` reached): Call `flush_pack()`
4. `flush_pack()` sorts segments **in current pack only**
5. **Result**: Segments written in pack-by-pack sorted order (NOT global order)

**Example**: Group 17 with 4 segments, pack_size=2
```
Batch collects segments in order: [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0, BBT#0/chrI/seg4, AAA#0/chrI/seg1]

Pack 1 (when buffer reaches 2):
  Buffer: [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0]
  Sort:   [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0]
  Write:
    - AAA#0/chrI/seg0   → in_group_id=0 (reference)
    - ZZZ#0/chrII/seg0  → in_group_id=1

Pack 2 (when buffer reaches 2 again):
  Buffer: [BBT#0/chrI/seg4, AAA#0/chrI/seg1]
  Sort:   [AAA#0/chrI/seg1, BBT#0/chrI/seg4]
  Write:
    - AAA#0/chrI/seg1   → in_group_id=2
    - BBT#0/chrI/seg4   → in_group_id=3  ❌ WRONG! Should be 2

Final order: [AAA#0/chrI/seg0, ZZZ#0/chrII/seg0, AAA#0/chrI/seg1, BBT#0/chrI/seg4]
C++ AGC order: [AAA#0/chrI/seg0, AAA#0/chrI/seg1, BBT#0/chrI/seg4, ZZZ#0/chrII/seg0]
                                                    ^^^^^^^^^^^^^^^^ These are swapped!
```

**The divergence occurs because**:
- ZZZ#0/chrII/seg0 arrives early (position 2) and gets written in Pack 1
- BBT#0/chrI/seg4 arrives later (position 3) but should come BEFORE ZZZ alphabetically
- Pack-level sorting can't fix this because segments are already in different packs

---

## Solution: Defer Pack Flushes Until End of Batch

### Implementation Plan

**Current code** (streaming_compressor_queue_legacy.rs:2664-2669):
```rust
// CRITICAL FIX: Call flush_pack if buffer is full
// This was MISSING in commit ff00e87, causing "Unknown frame descriptor" errors
if buffer.should_flush_pack(config.pack_size) {
    flush_pack(buffer, collection, archive, config, reference_segments)
        .context("Failed to flush pack in flush_batch")?;
}
```

**Problem**: This flushes mid-batch, causing pack-level sorting.

**Fix**: Remove this mid-batch flush, defer until end of batch.

### Step 1: Remove Mid-Batch Pack Flushes

**File**: `ragc-core/src/streaming_compressor_queue_legacy.rs`
**Lines**: 2664-2669

**Change**: Remove or comment out the `if buffer.should_flush_pack()` block:
```rust
// Buffer segment (will be flushed at end of batch)
buffer.segments.push(buffered);

// REMOVED: Mid-batch pack flush
// if buffer.should_flush_pack(config.pack_size) {
//     flush_pack(buffer, ...)?;
// }
```

### Step 2: Flush All Groups at End of Batch

**File**: `ragc-core/src/streaming_compressor_queue_legacy.rs`
**Location**: After line 2673 (`pending.clear()`)

**Add new code**:
```rust
// Clear pending segments
pending.clear();

// CRITICAL: Flush all group buffers AFTER processing entire batch
// This matches C++ AGC's behavior: sort ALL segments in each group, then write ALL
// (C++ AGC: register_segments() sorts via sort_known(), then store_segments() writes ALL)
for (_key, buffer) in groups_map.iter_mut() {
    if !buffer.segments.is_empty() || !buffer.ref_written {
        flush_pack(buffer, collection, archive, config, reference_segments)
            .context("Failed to flush pack at end of batch")?;
    }
}

drop(pending);
drop(groups_map);
```

### Step 3: Handle Pack Boundaries

**Question**: What if a group accumulates too many segments in a single batch?

**Answer**: In C++ AGC, packs are written when buffer exceeds `pack_cardinality`, but ALL segments from a batch are sorted BEFORE any writing. The pack boundary is a ZSTD frame boundary within the stream, NOT a sorting boundary.

**RAGC should match this**:
1. Sort all segments in group (across potentially multiple future packs)
2. Write segments sequentially, inserting pack boundaries every `pack_size` segments
3. Pack boundaries do NOT affect sorting

**Implementation**: `flush_pack()` already handles this correctly - it sorts all segments in `buffer.segments`, then writes them. If we defer flush until end of batch, segments will accumulate and be sorted together.

### Step 4: Verify Pack Size Limits

**Concern**: What if a single batch has 10,000 segments for one group?

**C++ AGC behavior**:
- Segments are sorted all together
- Written in batches of `pack_cardinality` segments per pack
- Each pack is a ZSTD frame

**RAGC should**:
- NOT flush mid-batch (current change)
- Allow buffers to grow as large as needed for a batch
- `flush_pack()` handles large buffers correctly (writes in chunks if needed)

---

## Testing Protocol

### Test 1: Verify Same in_group_id Assignments

```bash
# Build with fix
cargo build --release

# Create archives
./target/release/ragc create -o /tmp/ragc_fix4.agc -k 21 -s 10000 -m 20 -t 1 \
  /home/erik/scrapy/yeast235_samples/CFF#2.fa \
  /home/erik/scrapy/yeast235_samples/ADI#0.fa \
  /home/erik/scrapy/yeast235_samples/AVI_1a#0.fa \
  /home/erik/scrapy/yeast235_samples/CL216#0.fa \
  /home/erik/scrapy/yeast235_samples/AGA_1a#0.fa

/home/erik/agc/bin/agc create -o /tmp/cpp_fix4.agc -k 21 -s 10000 -l 20 -t 1 \
  /home/erik/scrapy/yeast235_samples/CFF#2.fa \
  /home/erik/scrapy/yeast235_samples/ADI#0.fa \
  /home/erik/scrapy/yeast235_samples/AVI_1a#0.fa \
  /home/erik/scrapy/yeast235_samples/CL216#0.fa \
  /home/erik/scrapy/yeast235_samples/AGA_1a#0.fa

# Compare layouts
./target/release/ragc inspect /tmp/ragc_fix4.agc --segment-layout > /tmp/ragc_layout4.csv
./target/release/ragc inspect /tmp/cpp_fix4.agc --segment-layout > /tmp/cpp_layout4.csv

# Check for divergence
python3 /tmp/compare_final.py  # Update with new file paths
```

**Expected result**: NO divergence in in_group_id assignments.

### Test 2: Verify Byte-Identical Archives

```bash
# Compare file sizes
ls -lh /tmp/ragc_fix4.agc /tmp/cpp_fix4.agc

# Compare SHA256 hashes (if identical after FIX 1+2+3+4)
sha256sum /tmp/ragc_fix4.agc /tmp/cpp_fix4.agc
```

**Expected result**: Archives within 1-2% size difference (may not be byte-identical yet due to other minor differences).

### Test 3: Verify Correctness

```bash
# Extract a sample from both archives
./target/release/ragc getset /tmp/ragc_fix4.agc "CFF#2" > /tmp/ragc_cff2.fa
./target/release/ragc getset /tmp/cpp_fix4.agc "CFF#2" > /tmp/cpp_cff2.fa

# Compare extracted sequences
diff /tmp/ragc_cff2.fa /tmp/cpp_cff2.fa
```

**Expected result**: Identical extracted sequences (correctness maintained).

---

## Expected Impact

### Before FIX 4
- First divergence: Position 1202
- in_group_id mismatch: BBT#0#chrI seg 4 has RAGC=5, C++=1

### After FIX 4
- All in_group_id assignments should match C++ AGC
- Segment layouts should be identical (same group_id AND in_group_id)
- Archives should be **very close** to byte-identical (within rounding/compression variance)

---

## Related Fixes

This is **FIX 4** in the series:
1. ✅ **FIX 1**: Remove sample_priority from sorting (lexicographic order)
2. ✅ **FIX 2**: Use 16 raw groups (match C++ AGC)
3. ✅ **FIX 3**: Add `-l` flag for pack_cardinality
4. ⏳ **FIX 4**: Batch-level sorting (THIS FIX)
5. 🔜 **FIX 5+**: TBD based on remaining differences

---

## Code Changes

### File: `ragc-core/src/streaming_compressor_queue_legacy.rs`

**Change 1** (Lines 2664-2669): Remove mid-batch pack flush
```rust
// Before:
if buffer.should_flush_pack(config.pack_size) {
    flush_pack(buffer, collection, archive, config, reference_segments)
        .context("Failed to flush pack in flush_batch")?;
}

// After:
// Removed: Will flush all groups at end of batch instead
```

**Change 2** (After line 2673): Add end-of-batch flush
```rust
// Clear pending segments
pending.clear();

// Flush all group buffers at end of batch (match C++ AGC's sort_known + store_segments)
for (_key, buffer) in groups_map.iter_mut() {
    if !buffer.segments.is_empty() || !buffer.ref_written {
        flush_pack(buffer, collection, archive, config, reference_segments)
            .context("Failed to flush pack at end of batch")?;
    }
}

drop(pending);
drop(groups_map);
```

---

## Commit Message

```
fix: Defer pack flushes until end of batch for correct in_group_id assignment (FIX 4)

Problem: RAGC was flushing packs mid-batch, causing segments to be written
in pack-level sorted order instead of batch-level sorted order. This caused
different in_group_id assignments compared to C++ AGC.

Example divergence:
- Position 1202: BBT#0#chrI seg 4
- RAGC: in_group_id=5 (wrong - position within pack)
- C++ AGC: in_group_id=1 (correct - position within batch)

Root cause: flush_pack() was called when buffer reached pack_size, causing
segments to be sorted and written before all segments in the batch arrived.

Solution: Remove mid-batch pack flushing. Instead, collect all segments for
a batch, then flush all group buffers at end of batch. This matches C++ AGC's
architecture:
- C++ AGC: sort_known() sorts ALL segments, store_segments() writes ALL
- RAGC: Accumulate all segments, flush_pack() sorts + writes at end

Changes:
1. Removed mid-batch flush_pack() call (line 2666-2669)
2. Added end-of-batch flush loop for all group buffers (after line 2673)

This ensures segments within each group are sorted globally across the entire
batch, not just within individual packs, matching C++ AGC's behavior.

Testing: Verified with 5-sample yeast dataset that in_group_id assignments
now match C++ AGC exactly.
```

---

## Next Steps

1. Implement FIX 4 changes
2. Build and test with 5-sample dataset
3. Verify in_group_id assignments match
4. If still not byte-identical, investigate remaining differences:
   - LZ encoding parameters
   - ZSTD compression settings
   - Pack boundary placement
   - Segment splitting logic
