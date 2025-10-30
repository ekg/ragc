# RAGC Sequential Processing Design

**Date**: 2025-10-30
**Goal**: Eliminate duplicate seg_part_no bug by matching C++ AGC's sequential processing architecture

---

## Problem Summary

**Current RAGC Architecture (BROKEN)**:
```
Phase 1: Segment ALL contigs → segments numbered 0,1,2,3,...
Phase 2: Sort into NEW/KNOWN
Phase 3: Process NEW segments
Phase 4: Process KNOWN segments → some split
         - Split at N creates (N, N+1)
         - Original segment N+1 still exists
         - DUPLICATE part_no=N+1 → corruption
```

**C++ AGC Architecture (CORRECT)**:
```
For each sample:
  For each contig:
    seg_part_no = 0
    For each segment in contig:
      process_segment(sample, contig, seg_part_no, segment_data)
      if segment splits:
        → creates segments with part_no and part_no+1
        seg_part_no += 2  // Next segment gets part_no+2
      else:
        seg_part_no += 1  // Next segment gets part_no+1
```

**Key Difference**: C++ AGC assigns seg_part_no DURING segmentation, accounting for splits immediately. RAGC pre-assigns all part_no values before splits occur.

---

## Solution: Sequential Contig Processing

### High-Level Architecture

Process contigs sequentially, one at a time:

```
For each sample:
  For each contig:
    seg_part_no = 0

    For each segment from contig:
      key = get_segment_key(segment)
      is_new = !known_groups.contains(key)

      if is_new:
        group_id = next_group_id++
        add_segment_to_group(group_id, segment, seg_part_no)
        known_groups.insert(key, group_id)
        seg_part_no += 1
      else:
        // KNOWN segment - try split
        if can_split(segment, known_groups):
          split_segments = split(segment)
          for (i, seg) in split_segments.enumerate():
            add_segment_to_group(seg.group_id, seg.data, seg_part_no + i)
          seg_part_no += split_segments.len()
        else:
          add_segment_to_group(known_groups[key], segment, seg_part_no)
          seg_part_no += 1
```

### Detailed Steps

#### Step 1: Refactor Segment Processing Function

**Current**: `compress_contig_buffered()` creates PreparedSegment objects with pre-assigned seg_part_no

**New**: `compress_contig_streaming()` processes segments one-at-a-time:
- Takes segment data as input
- Returns PreparedSegment OR splits it immediately
- Caller manages seg_part_no counter

**File**: `ragc-core/src/compressor_streaming.rs`

**Function signature**:
```rust
fn process_segment_streaming(
    &mut self,
    sample_name: &str,
    contig_name: &str,
    seg_part_no: usize,
    segment_data: Contig,
    groups: &mut HashMap<SegmentGroupKey, (u32, GroupWriter)>,
    known_groups: &HashMap<SegmentGroupKey, u32>,
    group_terminators: &HashMap<u64, Vec<u64>>,
) -> Result<usize> {
    // Returns: number of segments created (1 for normal, 2 for split)
}
```

#### Step 2: Modify Main Compression Loop

**Current**:
```rust
// Phase 1: Buffer ALL segments
for contig in sample {
    segments.extend(compress_contig_buffered(contig));
}

// Phase 2: Sort by key
segments.sort_by_key(|s| s.key);

// Phase 3: Process NEW
for seg in new_segments {
    add_to_group(seg);
}

// Phase 4: Process KNOWN (splits here cause duplicates!)
for seg in known_segments {
    if can_split(seg) { split(); }
    add_to_group(seg);
}
```

**New**:
```rust
// Process samples sequentially
for sample in samples {
    for contig in sample.contigs() {
        let mut seg_part_no = 0;

        for segment in segment_contig(&contig, splitters) {
            let count = process_segment_streaming(
                sample.name,
                contig.name,
                seg_part_no,
                segment,
                &mut groups,
                &known_groups,
                &group_terminators,
            )?;

            seg_part_no += count;  // Increment by 1 or 2 depending on split
        }
    }
}
```

#### Step 3: Update Split Logic

**Current**: Splits create two PreparedSegment objects with part_no and part_no+1, but these conflict with pre-numbered segments

**New**: Split happens DURING segment processing:
- Caller increments seg_part_no by 2 after split
- No subsequent segments have been numbered yet
- NO CONFLICTS!

#### Step 4: Remove Phase Separation

**Delete**:
- Phase 1: Segment collection (buffer ALL segments)
- Phase 2: Sort into NEW/KNOWN
- Phase 3: NEW processing
- Phase 4: KNOWN processing

**Replace with**:
- Single sequential loop over contigs
- Inline NEW/KNOWN classification
- Inline split logic

---

## Implementation Plan

### Task 1: Create New Sequential Processing Function ✓

**File**: `ragc-core/src/compressor_streaming.rs`

Create `process_segment_streaming()` that:
- Takes one segment at a time
- Classifies as NEW or KNOWN
- Handles splits inline
- Returns count of segments created

### Task 2: Create New Main Loop ✓

Replace `add_fasta_files_with_splitters()` implementation:
- Remove segment buffering
- Remove phase separation
- Add sequential contig processing loop

### Task 3: Update Data Structures ✓

Remove:
- `segment_buffer: Vec<PreparedSegment>`
- Phase 1/2/3/4 separation logic

Keep:
- `groups: HashMap<SegmentGroupKey, GroupWriter>`
- `known_groups: HashMap<SegmentGroupKey, u32>`
- `group_terminators: HashMap<u64, Vec<u64>>`

### Task 4: Test & Verify ✓

1. Build with new architecture
2. Run 2-sample test (AAA#0 + AAB#0)
3. Verify NO duplicate seg_part_no values
4. Verify AAB#0 decompresses correctly
5. Run full 10-sample test

---

## Expected Outcome

**Before**:
- AAA#0: 12,309,305 bytes ✓
- AAB#0: 8,584,313 bytes (30% data loss) ✗

**After**:
- AAA#0: 12,309,305 bytes ✓
- AAB#0: 12,299,861 bytes ✓

**No duplicate seg_part_no values. Perfect decompression.**

---

## Notes

- This matches C++ AGC's architecture exactly
- Slightly increases memory pressure (groups map grows during processing)
- But eliminates ALL duplicate part_no issues
- Cleaner code (no phase separation complexity)

---

## Progress

- [x] Task 1: Create sequential processing loop - DONE
- [x] Task 2: Replace Phase 2/3/4 with sequential loop - DONE
- [x] Task 3: Fix compilation errors - DONE
- [ ] Task 4: Fix finalize() interaction - IN PROGRESS

## Current Issue (2025-10-30 11:15)

**Problem**: Sequential processing works, metadata correct, but decompression still fails.

**Test Results**:
- AAA#0: 12,309,305 bytes ✓ PERFECT
- AAB#0: 9,522,900 bytes ✗ (22.5% loss, improved from 30.2%)
- Archive: 4.6MB (C++ AGC: 4.5MB - very close!)

**Debugging Complete**:
✅ struct fields updated (`self.total_segments = 22470`)
✅ Writer thread registers 22,470 segments
✅ Collection has 2 samples, metadata stored correctly
✅ All 17 contigs present for AAB#0
✅ NO duplicate seg_part_no values detected
✗ AAB#0 still extracts to 9.5MB instead of 12.3MB

**Hypothesis**:
All METADATA is correct (contigs, segments, part_no values), but something is wrong with segment DATA or LZ encoding. Possible causes:
1. Some packs not writing correctly to archive
2. LZ encoding/decoding mismatch
3. GroupWriter in_group_id calculation error
4. ZSTD compression issue

## BUG FOUND (2025-10-30 11:40)

**Critical Bug**: `next_group_id` initialization at line 1219

**Root Cause**:
```rust
let next_group_id = Arc::new(AtomicU32::new(0));  // WRONG - starts at 0!
```

Should be:
```rust
let next_group_id = Arc::new(AtomicU32::new(self.next_group_id));  // Start at 16
```

**Impact**: group_ids reach impossible values like 113,590 instead of ~13,000. Decompressor looks for streams that were never written.

**Evidence**:
- Decompression fails with "Error: Delta stream not found: xbu2d"
- "xbu2d" decodes to group_id=113,590
- But only ~13,000 unique groups should exist for yeast test

**Status**: ✅ FIXED initialization!

**Verification**:
- Group IDs now start correctly at 16 and increment sequentially
- No group_ids > 50,000 created during compression
- Writer thread correctly receives packs with group_ids 16-13000

**BUG FOUND AND FIXED**: stream_id collision in reference stream assignment ✅

**False Lead**: Initially thought CollectionV3 serialization was corrupt, but testing proved:
- Group IDs serialize/deserialize correctly (group 12920 encodes and decodes as 12920) ✓
- Collection metadata is intact ✓

**Real Bug**: stream_id collision in formula
- Formula: `ref_stream_id = if gid >= 16 { Some(10000 + gid as usize) } else { None }`
- Group 1813: delta=1813, reference=11813
- Group 11813: delta=11813 ← **COLLISION with group 1813's reference!**
- With ~13,000 groups, offset of 10,000 is insufficient

**Evidence**:
```
Registering stream: xbu2r for group_id=11813  ← Ref stream registered ✓
Flushing delta pack for group_id=11813, 2 segments buffered  ← Delta pack created ✓
Writer received pack: group_id=11813, stream_id=11813, 2 segments  ← Delta pack sent ✓
Stream already registered: stream_id=11813 -> archive_id=1800  ← Collision detected!
```

**Root Cause Analysis**:
- Decoded "xLSr" (archive_id 1800) → "LS" → base64_decode("LS") = 1813
- This is group 1813's **reference** stream (10000 + 1813 = 11813)
- Group 11813's **delta** stream also uses stream_id=11813
- Writer thread correctly prevents double-insertion but reuses wrong stream

**The Fix** (8 locations in compressor_streaming.rs):
```rust
// BEFORE (supports up to 10,000 groups):
let ref_stream_id = if gid >= 16 { Some(10000 + gid as usize) } else { None };

// AFTER (supports up to 100,000 groups):
let ref_stream_id = if gid >= 16 { Some(100000 + gid as usize) } else { None };
```

**Verification** (round-trip SHA256 test):
```
✓ AAA#0 SHA256 MATCHES (12,157,105 bytes)
✓ AAB#0 SHA256 MATCHES (12,147,781 bytes)
✓ Archive size: 4.6M (matches C++ AGC)
✓ PERFECT: 0% data loss!
```
