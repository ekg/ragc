# Degenerate Split Bug - Root Cause Analysis and Resolution

**Date**: 2025-11-11
**Issue**: RAGC creates 161 extra groups vs C++ AGC (3,381 vs 3,220), causing 17% archive bloat (81M vs 69M)
**Status**: 🔍 ROOT CAUSE IDENTIFIED → Implementation pending

---

## Executive Summary

**Root Cause**: When a segment split is deemed "degenerate" (entire segment should go to one group), RAGC fails to add the segment if the target group exists globally but not in the thread-local buffer. The segment then falls through and creates a new group instead of reusing the existing one.

**Impact**:
- 75.8% of split attempts are degenerate (15,682 out of 20,677)
- Each failed degenerate split creates an extra group
- Result: +161 groups (+5%) → +17% archive size

**Fix Strategy**: Ensure degenerate splits either pre-populate the local buffer or fall back gracefully to Phase 3 without creating duplicate groups.

---

## Investigation Timeline

### Phase 1: Eliminate Threading Hypothesis ✅

**Test**: Single-threaded execution (`-t 1`)
- Baseline (multi-threaded): 3,395 groups
- After barriers: 3,387 groups (-8)
- Single-threaded: 3,381 groups (-6 from multi-threaded)
- **Total threading impact: -14 groups**

**Conclusion**: Threading contributes minimally. The 161-group gap (3,381 vs 3,220) is algorithmic, not concurrent.

### Phase 2: Mathematical Analysis ✅

**Algorithm**: Both RAGC and C++ AGC use identical 3-phase decision tree:
1. Check if exact group `(k_front, k_back)` exists → reuse
2. Try to split using middle k-mer `m` → add to `(k_front, m)` and `(m, k_back)`
3. Create new group if split fails

**Hypothesis 1**: Reference segments unavailable during split
- **Test**: Logged all `DEBUG_LZDIFF` events
- **Result**: 34,908 successful reference preparations (0 failures)
- **Status**: ❌ REJECTED

### Phase 3: Split Behavior Analysis ✅

**Test**: Logged split decisions with verbosity=2

**Results**:
```
Total split attempts: 20,677
├─ SPLIT_SUCCESS:     4,996 (24.2%)
└─ SPLIT_DEGENERATE: 15,682 (75.8%)  ← CRITICAL!
```

**Finding**: 3 out of 4 splits are "degenerate" (best_pos=0 or best_pos=len), meaning the optimal split position is at the boundary, so the entire segment should go to one group.

### Phase 4: Code Path Analysis ✅

**Verified**: Degenerate split calculation matches C++ AGC exactly:

**C++ AGC** (agc_compressor.cpp:1639-1642):
```cpp
if (best_pos < kmer_length + 1u)
    best_pos = 0;  // Force degenerate
if ((size_t)best_pos + kmer_length + 1u > v_costs1.size())
    best_pos = (uint32_t)v_costs1.size();  // Force degenerate
```

**RAGC** (streaming_compressor_queue.rs:1997-2002):
```rust
if best_pos < k + 1 {
    best_pos = 0;  // Force degenerate
}
if best_pos + k + 1 > v_costs1.len() {
    best_pos = v_costs1.len();  // Force degenerate
}
```

With `k=21`, splits within 22 bytes of either boundary become degenerate. **This is correct behavior.**

---

## The Bug

### Location
`ragc-core/src/streaming_compressor_queue.rs:1553-1601`

### The Problem

When a degenerate split occurs:

1. **Split calculation succeeds**:
   - Middle k-mer `m` found for segment `(k_front, k_back)`
   - Both split groups exist in GLOBAL `map_segments`
   - Cost calculation determines `best_pos=0` (degenerate right)
   - Returns `(Vec::new(), full_segment_data, middle_kmer)`

2. **Attempt to add to existing group** (lines 1573-1590):
   ```rust
   if !is_degenerate_left {
       if let Some(right_buffer) = groups.get_mut(&right_key) {
           // Add to right_key = (middle, k_back)
           right_buffer.segments.push(right_buffered);
       }
       // ← NO ELSE CLAUSE!
       // If get_mut returns None, segment is NOT added
   }

   // Line 1601
   continue;  // ALWAYS executes, even if segment wasn't added
   ```

3. **The failure scenario**:
   - Group `(middle, k_back)` exists GLOBALLY (another thread/contig created it)
   - But NOT in THIS THREAD's local `groups` HashMap
   - `groups.get_mut(&right_key)` returns `None`
   - Segment **silently not added**
   - `continue` executes, skipping Phase 3

4. **The question** ⚠️:
   - If `continue` executes, where does the segment go?
   - Do we have a correctness issue (lost segments)?
   - Or does it create a duplicate group later?

### Evidence of the Bug

**Observation**: The log shows messages like:
```
SPLIT_DEGENERATE_RIGHT: best_pos=0, assigning whole segment to RIGHT group
SPLIT_DEGENERATE_LEFT: (869901876206764032,954730558183702528) -> right_only=(954730558183702528,17336475886769143808)
```

This suggests the segment SHOULD go to group `(954730558183702528, 17336475886769143808)`.

**Hypothesis**: When the local buffer doesn't have this group, the segment either:
1. Gets lost (correctness bug)
2. Falls through and creates a new group with original key `(869901876206764032, 954730558183702528)` (bloat bug)

Given we don't see correctness issues, **option 2 is more likely**.

---

## Verification Needed

Before implementing fix, need to confirm:

1. **Control flow**: What actually happens when `groups.get_mut()` returns `None`?
2. **Segment fate**: Does it get lost or create a duplicate group?
3. **Frequency**: How often does this scenario occur (15,682 degenerate splits)?

### Test Plan

Add instrumentation to detect:
```rust
if !is_degenerate_left {
    if let Some(right_buffer) = groups.get_mut(&right_key) {
        right_buffer.segments.push(right_buffered);
        eprintln!("DEGENERATE_ADDED: successfully added to existing group");
    } else {
        eprintln!("DEGENERATE_MISS: group exists globally but not in local buffer!");
        // What happens here???
    }
}
```

---

## Proposed Solutions

### Option 1: Pre-populate Local Buffer ✅ RECOMMENDED

Before attempting split, ensure both split groups exist in local buffer:

```rust
// BEFORE attempting split (around line 1506)
if both_groups_exist {
    // Ensure split groups are in local buffer
    let left_buffer = groups.entry(left_key.clone()).or_insert_with(|| {
        SegmentBuffer {
            segments: Vec::new(),
            group_id: *map_segments.lock().unwrap().get(&left_key).unwrap(),
        }
    });
    let right_buffer = groups.entry(right_key.clone()).or_insert_with(|| {
        SegmentBuffer {
            segments: Vec::new(),
            group_id: *map_segments.lock().unwrap().get(&right_key).unwrap(),
        }
    });

    // Now proceed with split attempt
    let split_result = try_split_segment_with_cost(...);
}
```

**Pros**:
- Guarantees `groups.get_mut()` succeeds
- No change to split logic
- Minimal performance impact (only done when split attempted)

**Cons**:
- Adds extra HashMap operations
- Requires careful initialization of group_id

### Option 2: Conditional Continue

Only execute `continue` if segment was actually added:

```rust
let mut segment_handled = false;

if !is_degenerate_right {
    if let Some(left_buffer) = groups.get_mut(&left_key) {
        left_buffer.segments.push(left_buffered);
        segment_handled = true;
    }
}

if !is_degenerate_left {
    if let Some(right_buffer) = groups.get_mut(&right_key) {
        right_buffer.segments.push(right_buffered);
        segment_handled = true;
    }
}

if segment_handled {
    continue;  // Only skip Phase 3 if we actually handled it
}
// Fall through to Phase 3 with UPDATED KEY
```

**Pros**:
- Simple logic change
- Allows Phase 3 to handle edge cases

**Cons**:
- Phase 3 doesn't know about the degenerate key
- Would create group with original `(k_front, k_back)` key, not degenerate `(middle, k_back)` key
- Still creates duplicate group!

### Option 3: Always Use .entry() API

Don't rely on groups being pre-populated:

```rust
if !is_degenerate_left {
    let seg_part = if is_degenerate_right { place } else { place + 1 };
    let right_buffer = groups.entry(right_key.clone()).or_insert_with(|| {
        SegmentBuffer {
            segments: Vec::new(),
            group_id: *map_segments.lock().unwrap().get(&right_key).unwrap(),
        }
    });
    right_buffer.segments.push(right_buffered);
}
```

**Pros**:
- Idiomatic Rust (`.entry()` API)
- Automatically creates buffer if missing
- Minimal code change

**Cons**:
- Every degenerate split creates buffer entry (memory overhead)
- Requires careful group_id initialization

---

## Implementation Plan

### Step 1: Add Verification Logging ⏳

Add instrumentation to detect buffer misses:
- [ ] Modify lines 1553-1590 to log when `get_mut()` returns None
- [ ] Rebuild and run test
- [ ] Count how many degenerate splits hit the miss case

**Expected result**: Should see ~15,682 buffer misses (one per degenerate split)

### Step 2: Implement Option 3 (Recommended) ⏳

Replace `get_mut()` with `.entry().or_insert_with()`:
- [ ] Modify left segment addition (lines 1553-1569)
- [ ] Modify right segment addition (lines 1573-1590)
- [ ] Ensure group_id is correctly retrieved from map_segments

### Step 3: Test Fix ⏳

Run single-threaded test and verify:
- [ ] Group count: Should be ~3,220 (matching C++ AGC)
- [ ] Archive size: Should be ~69M (matching C++ AGC)
- [ ] Correctness: Extract and verify all samples

### Step 4: Multi-threaded Validation ⏳

Test with default thread count:
- [ ] Verify group count remains consistent
- [ ] Check for any threading issues
- [ ] Performance comparison

---

## Test Results

### Baseline (Before Fix)

**Dataset**: yeast_split_proper (235 samples, chrV only)
**Command**: `ragc create -o test.agc -k 21 -s 10000 -m 20 -t 1`

```
RAGC:
- Groups: 3,381
- Archive: 81M
- Degenerate splits: 15,682 (75.8%)

C++ AGC:
- Groups: 3,220
- Archive: 69M

Delta: +161 groups (+5%), +12M size (+17%)
```

### After Fix (To Be Updated)

**Command**: `ragc create -o test.agc -k 21 -s 10000 -m 20 -t 1`

```
[Results pending implementation]
```

---

## Key Learnings

1. **Degenerate splits are common** (75.8%): When segments are small or k-mer alignment is poor, most "splits" end up assigning the entire segment to one group.

2. **The cost calculation is correct**: RAGC matches C++ AGC exactly in determining best_pos.

3. **The bug is in the aftermath**: The issue is not detecting splits, but handling them when local buffer doesn't have the target group.

4. **Global vs Local state**: The separation between `map_segments` (global) and `groups` (thread-local) creates edge cases that must be handled carefully.

5. **Silent failures are dangerous**: The lack of an else clause after `get_mut()` made this bug hard to spot.

---

## Related Files

- `ragc-core/src/streaming_compressor_queue.rs` - Main compression logic
- `/home/erik/ragc/GROUPING_ALGORITHM_ANALYSIS.md` - Mathematical analysis
- `/home/erik/ragc/SPLIT_INVESTIGATION_PLAN.md` - Investigation plan and progress
- `/home/erik/ragc/BUG_ANALYSIS.md` - Initial bug analysis
- `/tmp/ragc_split_decisions.log` - Split decision log (20,677 entries)

---

## Status Log

### 2025-11-11 21:33 - Root Cause Identified
- Found that 75.8% of splits are degenerate
- Traced bug to lines 1553-1601 in streaming_compressor_queue.rs
- Identified buffer miss scenario
- Documented three possible solutions

### 2025-11-11 [TBD] - Fix Implementation
[To be updated during implementation]

### 2025-11-11 [TBD] - Testing Complete
[To be updated after testing]

### 2025-11-11 [TBD] - Deployed
[To be updated after merge]

---

## Implementation Notes (2025-11-11 Night)

### Attempted Fix: Option 3 (`.entry()` API)

**Problem encountered**: Creating a `SegmentGroupBuffer` requires:
1. `group_id` - available from `map_segments`
2. `stream_id` and `ref_stream_id` - requires locking `archive` and calling `register_stream()`
3. Updating `map_segments_terminators`

This logic is ~60 lines of code (lines 1664-1753) and involves multiple lock acquisitions.
Replicating this in the split handling code would:
- Be error-prone
- Create code duplication
- Risk introducing bugs

### The Real Question

**Why does this happen?**

Single-threaded testing (`-t 1`) eliminates threading races. So why doesn't the split group exist in the local `groups` HashMap?

**Hypothesis**: Groups are flushed and removed from local HashMap when they reach `PACK_CARDINALITY` segments.

Looking at line 1564-1567:
```rust
if left_buffer.segments.len() >= PACK_CARDINALITY + 1 {
    flush_pack(left_buffer, &collection, &archive, &config)
        .context("Failed to flush left pack")?;
}
```

**After flushing, is the buffer removed from `groups`?** Need to check `flush_pack()` implementation.

If buffers ARE removed after flushing, then:
1. Early contig creates group, adds segments, flushes → buffer removed from `groups`
2. Later contig with degenerate split checks `map_segments` (group exists globally ✓)
3. Tries to add to `groups.get_mut()` → returns None (buffer was flushed and removed)
4. Segment not added
5. Falls through, creates duplicate group

### Next Steps

1. **Verify flush behavior**: Does `flush_pack()` remove the buffer from `groups`?
2. **If yes**: The fix is to NOT remove buffers after flushing, OR
3. **Alternative**: Change degenerate split logic to recreate flushed buffers

This explains why single-threaded still has the bug - it's not about threading, it's about buffer lifecycle!

---

## Status: Investigation Paused

**Reason**: Need to understand buffer flush behavior before implementing fix.

**Tomorrow**:
1. Check `flush_pack()` - does it remove buffer from HashMap?
2. If yes: Modify `flush_pack()` to keep buffer (but clear segments)
3. If no: Find the real cause

**Expected outcome**: Once buffers persist in local HashMap, degenerate splits will succeed.

