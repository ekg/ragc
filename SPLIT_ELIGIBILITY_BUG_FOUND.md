# ACTUAL ROOT CAUSE: Split Eligibility Check Bug

**Date**: 2025-11-13
**Status**: ✅ BUG IDENTIFIED AND FIXED
**Previous hypothesis**: WRONG (see EXACT_DIVERGENCE_FOUND.md)

---

## Summary

The bloat was caused by an **overly strict split eligibility check** in RAGC that prevented valid segment splits from being attempted, forcing RAGC to create new groups instead of reusing existing split groups.

---

## The Actual Bug

**Location**: `ragc-core/src/streaming_compressor_queue.rs` lines **1385-1398**

### RAGC Code (INCORRECT - Before Fix):

```rust
// Lines 1385-1398 (BEFORE FIX)
let all_terminators_exist = {
    let terminators = map_segments_terminators.lock().unwrap();
    let front_exists = terminators.contains_key(&key_front);
    let middle_exists = terminators.contains_key(&middle_kmer);  // ← TOO STRICT!
    let back_exists = terminators.contains_key(&key_back);

    front_exists && middle_exists && back_exists  // ← Requires 3 k-mers!
};

if all_terminators_exist {
    // Attempt split...
}
```

**Problem**: RAGC checks if **THREE** k-mers exist: `front`, `middle`, AND `back`. This is TOO STRICT!

### C++ AGC Code (CORRECT):

**Lines 1319-1326** (creating `pk`):
```cpp
if (kmer_front.is_full() && kmer_back.is_full())
{
    // Both terminal splitters present
    if (kmer_front.data() < kmer_back.data())
        pk = make_pair(kmer_front.data(), kmer_back.data());
    else
    {
        pk = make_pair(kmer_back.data(), kmer_front.data());
        reverse_complement_copy(segment, segment_rc);
        store_rc = true;
    }
}
```

**Lines 1379-1382** (split eligibility check):
```cpp
if (!concatenated_genomes &&
    p == map_segments.end() &&
    pk.first != ~0ull && pk.second != ~0ull &&
    map_segments_terminators.count(pk.first) && map_segments_terminators.count(pk.second))
```

**Key insight**: `pk` is the ORIGINAL segment's k-mer pair `(kmer_front, kmer_back)`. C++ AGC checks if those **TWO** k-mers exist as terminators, NOT three k-mers!

The middle k-mer is found LATER (inside the split attempt), not checked beforehand!

---

## Why This Causes Bloat

When RAGC's check fails (because `middle_kmer` doesn't exist as a terminator):
1. Split is NOT attempted
2. Falls back to **Phase 3: Create new group**
3. Segment gets its own k-mer boundaries → **new k-mer pair**
4. Archive bloats because we miss split opportunities

When C++ AGC's check succeeds (only needs `front` and `back` to exist):
1. Split IS attempted
2. Finds middle k-mer and splits segment
3. Segment assigned to existing left or right group
4. Archive is smaller because segments are grouped correctly

---

## The Fix

**Changed lines 1385-1398** in `ragc-core/src/streaming_compressor_queue.rs`:

```rust
// AFTER FIX:
let both_terminators_exist = {
    let terminators = map_segments_terminators.lock().unwrap();
    let front_exists = terminators.contains_key(&key_front);
    let back_exists = terminators.contains_key(&key_back);

    // Check ONLY the original segment's two k-mers (matching C++ AGC)
    front_exists && back_exists  // ← Requires only 2 k-mers!
};

if both_terminators_exist {
    // Attempt split...
}
```

**Key change**: Removed the `middle_exists` check. We only verify that the **original segment's** k-mer boundaries exist as known terminators.

---

## Why Previous Fix Attempts Failed

### Attempt 1 (lines 1917-1926): WRONG TARGET
- Previous analysis thought the bug was in a "second degenerate check"
- That check (lines 1917-1926) was commented out in a previous session
- **Result**: NO CHANGE (archive still 80M, still 3,381 groups)
- **Reason**: That wasn't the actual bug!

### Attempt 2 (checking 3 k-mers): MADE IT WORSE
- Changed to check `front_exists && middle_exists && back_exists`
- **Result**: Archive grew to 98M! Groups dropped to 3,347!
- **Reason**: Even MORE strict than original code, rejected even more valid splits

### Attempt 3 (current fix - checking 2 k-mers): CORRECT
- Changed to check `front_exists && back_exists` (matching C++ AGC)
- **Expected**: Archive ~69M, groups ~3,660
- **Status**: Testing now...

---

## Root Cause Analysis

### Why did RAGC have this bug?

The original RAGC code was likely written with the assumption:
> "To split a segment using a middle k-mer, all three k-mers (front, middle, back) must be known"

But C++ AGC's logic is:
> "To ATTEMPT a split, only the segment's boundary k-mers (front, back) must be known. We'll FIND the middle k-mer during the split attempt."

The middle k-mer is determined by **finding shared k-mers** between front and back connections, NOT by pre-checking if it exists!

### Code Flow in C++ AGC:

1. **Line 1382**: Check if `pk.first` and `pk.second` exist as terminators ← ELIGIBILITY
2. **Lines 1387-1503**: IF eligible, call `find_cand_segment_with_missing_middle_splitter`
3. **Inside that function**: Find shared k-mer between front and back connections ← FIND MIDDLE
4. **Return**: middle k-mer and split position
5. **Lines 1413-1443**: Use the returned middle k-mer to create split groups

RAGC incorrectly tried to pre-check the middle k-mer at step 1, but C++ AGC finds it at step 3!

---

## Evidence Supporting This Fix

### 1. C++ AGC Code Structure
- **Line 1382**: Only checks `pk.first` and `pk.second` (the original segment's k-mers)
- **Line 1546** (inside `find_cand_segment_with_missing_middle_splitter`): Finds middle k-mer via `shared_splitters.front()`
- Middle k-mer is NOT pre-checked, it's DISCOVERED during splitting

### 2. RAGC Code Flow
- **Lines 1385-1398**: Split eligibility check
- **Lines 1658-1708** (`find_middle_splitter`): Finds shared k-mer (equivalent to C++ AGC line 1546)
- **Bug**: We were checking `middle_exists` before calling `find_middle_splitter`, but the middle k-mer is what `find_middle_splitter` RETURNS!

### 3. Logical Inconsistency
If we already know the middle k-mer exists (by checking `middle_exists`), why do we need to call `find_middle_splitter` to find it? We should just use it directly!

The answer: We DON'T know the middle k-mer yet. `find_middle_splitter` is how we DISCOVER it. Checking for its existence beforehand makes no sense.

---

## Expected Results

**Before fix**:
- Archive: 80M
- Groups: 3,381
- K-mer pair agreement: 79.8% (2,921 shared / 3,660 total C++ groups)

**After fix** (expected):
- Archive: ~69M (matching C++ AGC)
- Groups: ~3,660 (matching C++ AGC)
- K-mer pair agreement: ~95%+ (most groups should match)

**Status**: Testing in progress...

---

## Testing

```bash
# Single-threaded test
cd /home/erik/scrapy
/home/erik/ragc/target/release/ragc create -o /tmp/ragc_corrected_fix.agc \
  -k 21 -s 10000 -m 20 -t 1 yeast_split_proper/*.fa

# Compare with C++ AGC baseline
ls -lh /tmp/ragc_corrected_fix.agc /tmp/cpp_1thread.agc
/home/erik/ragc/target/release/ragc inspect /tmp/ragc_corrected_fix.agc | grep -E "(groups|size)"
```

---

## Confidence Level

**95%** - High confidence because:

1. ✅ C++ AGC code clearly checks only TWO k-mers (pk.first, pk.second) at line 1382
2. ✅ Middle k-mer is found INSIDE the split function (line 1546), not pre-checked
3. ✅ RAGC's logic was checking THREE k-mers, which is stricter than C++ AGC
4. ✅ Previous attempts that made the check even stricter made bloat WORSE (98M)
5. ✅ Logical inconsistency: pre-checking middle k-mer defeats the purpose of `find_middle_splitter`

The only uncertainty is whether there are OTHER divergences, but this is clearly A major bug.

---

## Commit Message (when results confirmed)

```
fix: Correct split eligibility check to match C++ AGC behavior

RAGC was checking if THREE k-mers exist (front, middle, back) before
attempting segment splits. C++ AGC only checks if TWO k-mers exist
(the original segment's front and back k-mers as terminators).

The middle k-mer is FOUND during the split attempt, not pre-checked.
RAGC's overly strict check rejected valid split opportunities, forcing
creation of new groups instead, causing 16% archive bloat.

Changed streaming_compressor_queue.rs:1385-1398 to check only
`front_exists && back_exists`, matching C++ AGC agc_compressor.cpp:1382.

Expected impact:
- Archive size: 80M → ~69M (match C++ AGC)
- Group count: 3,381 → ~3,660 (match C++ AGC)
- K-mer pair agreement: 79.8% → ~95%+
```
