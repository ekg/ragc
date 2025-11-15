# EXACT DIVERGENCE POINT IDENTIFIED

**Date**: 2025-11-12  
**Status**: ✅ BUG FOUND - Ready to fix

---

## The Bug

**Location**: `ragc-core/src/streaming_compressor_queue.rs` lines **1917-1926**

**RAGC has an EXTRA validation check that C++ AGC does NOT have!**

### RAGC Code (INCORRECT):

```rust
// Lines 1917-1926
if left_size < k + 1 || right_size < k + 1 || left_size > segment_data.len() {
    // Degenerate split - segments too small or out of bounds
    if config.verbosity > 1 {
        eprintln!(
            "SPLIT_DEGENERATE: best_pos={} left_size={} right_size={} (need {}-{}) - NOT splitting",
            best_pos, left_size, right_size, k + 1, segment_data.len()
        );
    }
    return None;  // ← RAGC REJECTS THE SPLIT!
}
```

**This check happens AFTER** the degenerate position adjustment (lines 1873-1878):

```rust
if best_pos < k + 1 {
    best_pos = 0;
}
if best_pos + k + 1 > v_costs1.len() {
    best_pos = v_costs1.len();
}
```

### C++ AGC Code (CORRECT):

**Lines 1637-1642**:

```cpp
if (best_pos < kmer_length + 1u)
    best_pos = 0;
if ((size_t)best_pos + kmer_length + 1u > v_costs1.size())
    best_pos = (uint32_t)v_costs1.size();

return make_pair(middle, best_pos);  // ← ALWAYS RETURNS!
```

**C++ AGC ALWAYS returns** the middle k-mer and best_pos, even if best_pos is 0 or full length!

Then at **lines 1413-1417**, C++ AGC handles the degenerate case:

```cpp
if (left_size == 0)
{
    // Assign whole segment to RIGHT group
    if (split_match.first < kmer2.data())
        store_rc = use_rc;
    else
        ...
}
```

---

## Why This Causes Bloat

When RAGC rejects a degenerate split:
1. Phase 2 (splitting) returns `None`
2. Fallback to Phase 3 (create new group)
3. Segment gets NEW k-mer boundaries based on its actual content
4. **Different k-mer pair than C++ AGC's split result**

When C++ AGC accepts a degenerate split:
1. Phase 2 returns the middle k-mer
2. Segment assigned to left OR right group (lines 1413-1443)
3. **Uses the middle k-mer as boundary** (from the split)
4. **Different k-mer pair than RAGC's new group**

### Example from logs:

**RAGC**: Creates group `(0x4f51386ca000000, 0xbff71fdced000000)`
- RAGC rejected the split (lines 1917-1926)
- Created new group with original boundaries

**C++ AGC**: NEVER creates this k-mer pair!
- C++ AGC accepted the split
- Used middle k-mer to create `(0xbff71fdced000000, ...)` and `(..., 0x4f51386ca000000)` groups
- Split the segment into two specialized groups

---

## Impact

- **571 segments** rejected by RAGC that C++ AGC accepts
- Each creates a new group with "wrong" k-mer boundaries
- Results in **739 extra C++ AGC groups** vs **168 extra RAGC groups**
- Net effect: C++ AGC has 279 more groups but 11M smaller archive

---

## The Fix

**Remove lines 1917-1926** in `ragc-core/src/streaming_compressor_queue.rs`.

The degenerate checks at lines 1885-1905 (checking `if left_size == 0` and `if right_size == 0`) are CORRECT and match C++ AGC's behavior. The SECOND check at lines 1917-1926 is WRONG and rejects splits that C++ AGC accepts.

### Correct Flow (matching C++ AGC):

1. Calculate best_pos (lines 1857-1868)
2. Apply degenerate position rules (lines 1873-1878):
   - If `best_pos < k+1` → force to 0
   - If `best_pos + k+1 > len` → force to len
3. Check if degenerate (lines 1885-1905):
   - If `left_size == 0` → return whole segment as RIGHT
   - If `right_size == 0` → return whole segment as LEFT
4. ~~Apply SECOND check (lines 1917-1926)~~ ← **DELETE THIS!**
5. Split at best_pos (line 1929)

---

## Files to Modify

1. **`ragc-core/src/streaming_compressor_queue.rs`**
   - **DELETE lines 1917-1926** (or comment them out)
   - Keep lines 1885-1905 (degenerate checks matching C++ AGC)

---

## Expected Results After Fix

- RAGC group count: Should increase from 3,381 to ~3,660 (match C++ AGC)
- RAGC archive size: Should decrease from 80M to ~69M (match C++ AGC)
- K-mer pair agreement: Should increase from 79.8% to ~95%+

---

## Testing Plan

1. Rebuild RAGC with the fix
2. Run single-threaded test:
   ```bash
   cd /home/erik/scrapy
   /home/erik/ragc/target/release/ragc create -o /tmp/ragc_fixed.agc \
     -k 21 -s 10000 -m 20 -t 1 yeast_split_proper/*.fa
   ```
3. Compare with C++ AGC baseline:
   ```bash
   ls -lh /tmp/ragc_fixed.agc /tmp/cpp_groups.agc
   # Should be: ~69M vs 69M
   ```
4. Verify correctness:
   ```bash
   /home/erik/ragc/target/release/ragc getset /tmp/ragc_fixed.agc AAA#0 > /tmp/ragc_AAA.fa
   /home/erik/agc/bin/agc getset /tmp/cpp_groups.agc AAA#0 > /tmp/cpp_AAA.fa
   diff /tmp/ragc_AAA.fa /tmp/cpp_AAA.fa  # Should be identical
   ```

---

## Root Cause Timeline

1. **Commit 001a64c**: Added degenerate split support (lines 1885-1905) ✅ CORRECT
2. **Unknown commit**: Added SECOND degenerate check (lines 1917-1926) ❌ BUG INTRODUCED
3. **Today**: Identified that C++ AGC does NOT have this second check

The second check was likely added as a "safety check" but it's WRONG - it rejects valid splits that C++ AGC accepts!

---

## Confidence Level

**99%** - The evidence is conclusive:

1. C++ AGC code clearly shows NO second check
2. C++ AGC code clearly shows it accepts degenerate splits (lines 1413-1417)
3. RAGC logs show it rejects splits with "SPLIT_DEGENERATE: ... - NOT splitting"
4. The k-mer pair divergence matches exactly what we'd expect from this bug

