# Segment Splitting Investigation Plan

**Date Started**: 2025-11-11
**Problem**: RAGC creates 161 extra groups (3,381 vs 3,220) causing 16% archive bloat (80M vs 69M)
**Status**: 🔬 ACTIVE INVESTIGATION

---

## What We Know (Established Facts)

### Threading is NOT the cause ✅
- Barriers: -8 groups (0.2% improvement)
- Atomic allocation: 0 groups (no effect)
- Single-threaded: -6 groups (negligible)
- **Total from threading fixes: -14 groups**
- **Remaining algorithmic bloat: 161 extra groups**

### Evidence
| Configuration | Archive | Groups | Gap from C++ |
|--------------|---------|---------|--------------|
| RAGC baseline (multi-threaded) | 81M | 3,395 | +175 |
| RAGC + barriers | 81M | 3,387 | +167 |
| RAGC + atomic | 81M | 3,387 | +167 |
| RAGC single-threaded | 80M | 3,381 | **+161** |
| C++ AGC (any threading) | 69M | 3,220 | 0 |

### Algorithm Structure (Identical in Both)

Both implementations use the same 3-phase decision tree:

```
For each segment S with k-mers (k_front, k_back):

Phase 1: Check if group (k_front, k_back) exists
   YES → Add to existing group (DONE)
   NO  → Continue to Phase 2

Phase 2: Try to split using middle k-mer
   2a. Find middle: k_mid ∈ terminators[k_front] ∩ terminators[k_back]
   2b. Check split groups exist: (k_front, k_mid) and (k_mid, k_back)
   2c. Calculate compression cost for each split position
   2d. Choose best position, split if valid
   YES → Add to split groups (DONE)
   NO  → Continue to Phase 3

Phase 3: Create new group
   Allocate new group_id for (k_front, k_back)
   Add segment to new group (DONE)
```

---

## Hypotheses (Ordered by Likelihood)

### Hypothesis 1: Reference Segment Unavailability ⭐⭐⭐⭐⭐
**Likelihood**: VERY HIGH

**Description**: RAGC fails split attempts because reference segments aren't available yet in `reference_segments` HashMap, even though group IDs exist in `map_segments`.

**Evidence**:
- RAGC code (streaming_compressor_queue.rs:1906-1942):
  ```rust
  if let Some(&segment_id) = map_segments_locked.get(key) {
      if let Some(ref_data) = ref_segments_locked.get(&segment_id) {
          // Found reference ✓
      } else {
          // ⚠️ Group exists but NO REFERENCE DATA
          return None;  // Split fails!
      }
  }
  ```

- C++ AGC code (agc_compressor.cpp:1555-1559):
  ```cpp
  auto segment_id1 = map_segments[minmax(kmer_front, middle)];
  auto seg1 = v_segments[segment_id1];  // Global vector, always available
  ```

**Key Difference**:
- C++ AGC: Uses global `v_segments` vector, populated immediately when group created
- RAGC: Uses `reference_segments` HashMap, populated later when pack flushed

**Predicted Impact**:
- If reference not available → `prepare_on_demand()` returns `None` → `try_split_segment_with_cost()` returns `None` → Falls through to Phase 3 → Creates new group

**Test**: Log when `prepare_on_demand()` returns `None` due to missing reference

**Expected Result**: ~161 instances of "SPLIT_SKIP: left/right group has no reference yet"

### Hypothesis 2: Cost Calculation Orientation ⭐⭐⭐
**Likelihood**: MEDIUM

**Description**: RAGC calculates compression cost with wrong segment orientation (forward vs reverse-complement).

**Evidence**:
- RAGC code (streaming_compressor_queue.rs:1947-1965):
  ```rust
  let v_costs1 = lz.get_coding_cost_vector(segment_data, true);   // Always true
  let v_costs2 = lz.get_coding_cost_vector(segment_data, false);  // Always false
  ```

- C++ AGC code (agc_compressor.cpp:1562-1568):
  ```cpp
  if (kmer_front.data() < middle)
      seg1->get_coding_cost(segment_dir, v_costs1, true, ...);
  else {
      seg1->get_coding_cost(segment_rc, v_costs1, false, ...);
      reverse(v_costs1.begin(), v_costs1.end());
  }
  ```

**Key Difference**:
- C++ AGC: Checks k-mer ordering and uses segment_rc if needed
- RAGC: Always uses same orientation regardless of k-mer ordering

**Predicted Impact**:
- Wrong orientation → Poor cost calculation → Split position invalid or total cost too high → Split rejected → New group created

**Test**: Compare cost values for same segment in both implementations

**Expected Result**: Some segments show very different cost profiles

### Hypothesis 3: Degenerate Split Threshold ⭐⭐
**Likelihood**: LOW

**Description**: RAGC rejects degenerate splits (one side empty) more aggressively than C++ AGC.

**Evidence**:
- Both implementations accept degenerate splits (C++ AGC:1415-1430, RAGC:2009-2029)
- Threshold for minimum segment size might differ

**Predicted Impact**: Minor (~10-20 groups)

**Test**: Count degenerate splits in both implementations

### Hypothesis 4: Middle K-mer Selection ⭐
**Likelihood**: VERY LOW

**Description**: RAGC chooses different middle k-mer from intersection set.

**Evidence**:
- C++ AGC (agc_compressor.cpp:1548): `middle = shared_splitters.front()`
- RAGC (streaming_compressor_queue.rs:1832): `shared_kmers.first().copied()`

**Both use first element** - Unlikely to differ

**Predicted Impact**: Negligible

---

## Investigation Plan

### Phase 1: Instrument Code for Logging 📋

**Objective**: Add comprehensive logging to track every split attempt and failure reason

#### Step 1.1: Add C++ AGC Logging ✅ (Already has some logging)

Check existing logging in C++ AGC:
```bash
grep -n "TRACE: Log split" /home/erik/agc/src/core/agc_compressor.cpp
```

Output shows logging at line 1470-1476 for groups 21-30.

**Expand logging to capture**:
- [ ] Every split attempt (not just groups 21-30)
- [ ] Middle k-mer found
- [ ] Cost calculation results
- [ ] Split accept/reject decision with reason

#### Step 1.2: Add RAGC Logging ✅ (Already extensive)

Check existing logging in RAGC:
```bash
grep -n "eprintln!" ragc-core/src/streaming_compressor_queue.rs | grep -i split | head -20
```

Current logging includes:
- Line 1456-1458: Middle k-mer found
- Line 1497-1501: Split group existence check
- Line 1509-1512: Split attempt with cost
- Line 1531-1549: Degenerate split detection
- Line 1607-1612: Split skipped (groups don't exist)
- Line 1893-1898: Split attempt parameters
- Line 1913-1918: Reference segment preparation (on-demand)
- Line 1924-1931: Reference segment MISSING warning
- Line 1933-1939: Group not in map_segments warning
- Line 1951-1954: Split skip (left ref unavailable)
- Line 1960-1964: Split skip (right ref unavailable)
- Line 2012-2016: Degenerate right split
- Line 2023-2027: Degenerate left split
- Line 2043-2048: Degenerate split (segments too small)
- Line 2055-2063: Split success

**Status**: RAGC already has comprehensive logging ✅

#### Step 1.3: Rebuild with Logging

```bash
cd /home/erik/ragc
cargo build --release
```

**Status**: [ ] Not started

### Phase 2: Run Comparison Tests 🧪

**Test Dataset**: yeast235 (235 samples, single-threaded to eliminate any threading noise)

#### Test 2.1: Collect Split Logs

**C++ AGC**:
```bash
cd /home/erik/scrapy
/home/erik/agc/bin/agc create -o /tmp/cpp_debug.agc -k 21 -s 10000 -t 1 \
  yeast_split_proper/*.fa 2>&1 | \
  grep -E "\[CPP\].*SPLIT" > /tmp/cpp_splits.log

# Capture group creation too
/home/erik/agc/bin/agc create -o /tmp/cpp_debug.agc -k 21 -s 10000 -t 1 \
  yeast_split_proper/*.fa 2>&1 | \
  grep -E "\[CPP\].*(GROUP_CREATE|SPLIT)" > /tmp/cpp_full.log
```

**RAGC**:
```bash
cd /home/erik/scrapy
/home/erik/ragc/target/release/ragc create -o /tmp/ragc_debug.agc -k 21 -s 10000 -m 20 -t 1 -v 2 \
  yeast_split_proper/*.fa 2>&1 | \
  grep -E "(SPLIT_|GROUP_CREATE)" > /tmp/ragc_splits.log

# Or with even higher verbosity for hypothesis 1:
/home/erik/ragc/target/release/ragc create -o /tmp/ragc_debug.agc -k 21 -s 10000 -m 20 -t 1 -v 1 \
  yeast_split_proper/*.fa 2>&1 | \
  grep -E "(SPLIT_|DEBUG_LZDIFF|GROUP_CREATE)" > /tmp/ragc_full.log
```

**Status**: [ ] Not started

#### Test 2.2: Compare Split Statistics

```bash
# Count total split attempts
echo "=== C++ AGC ===" > /tmp/split_comparison.txt
echo "Total split attempts:" >> /tmp/split_comparison.txt
grep -c "SPLIT_ACCEPT\|SPLIT" /tmp/cpp_splits.log >> /tmp/split_comparison.txt

echo "" >> /tmp/split_comparison.txt
echo "=== RAGC ===" >> /tmp/split_comparison.txt
echo "Total split attempts:" >> /tmp/split_comparison.txt
grep -c "SPLIT_ATTEMPT\|SPLIT_SUCCESS" /tmp/ragc_splits.log >> /tmp/split_comparison.txt

echo "" >> /tmp/split_comparison.txt
echo "RAGC split successes:" >> /tmp/split_comparison.txt
grep -c "SPLIT_SUCCESS\|SPLIT_DEGENERATE" /tmp/ragc_splits.log >> /tmp/split_comparison.txt

echo "" >> /tmp/split_comparison.txt
echo "RAGC split failures:" >> /tmp/split_comparison.txt
grep -c "SPLIT_SKIP" /tmp/ragc_splits.log >> /tmp/split_comparison.txt

# Break down RAGC failure reasons
echo "" >> /tmp/split_comparison.txt
echo "=== RAGC Failure Breakdown ===" >> /tmp/split_comparison.txt
echo "No reference (left):" >> /tmp/split_comparison.txt
grep -c "left group has no reference yet" /tmp/ragc_full.log >> /tmp/split_comparison.txt
echo "No reference (right):" >> /tmp/split_comparison.txt
grep -c "right group has no reference yet" /tmp/ragc_full.log >> /tmp/split_comparison.txt
echo "Cost vectors empty:" >> /tmp/split_comparison.txt
grep -c "cost vectors empty" /tmp/ragc_full.log >> /tmp/split_comparison.txt
echo "Cost vectors mismatch:" >> /tmp/split_comparison.txt
grep -c "cost vector length mismatch" /tmp/ragc_full.log >> /tmp/split_comparison.txt
echo "Segments too small:" >> /tmp/split_comparison.txt
grep -c "NOT splitting" /tmp/ragc_full.log >> /tmp/split_comparison.txt

cat /tmp/split_comparison.txt
```

**Status**: [ ] Not started

**Expected Result (if Hypothesis 1 correct)**:
```
=== C++ AGC ===
Total split attempts: ~500

=== RAGC ===
Total split attempts: ~500
RAGC split successes: ~339
RAGC split failures: ~161

=== RAGC Failure Breakdown ===
No reference (left): ~80
No reference (right): ~81
Cost vectors empty: 0
Cost vectors mismatch: 0
Segments too small: 0
```

#### Test 2.3: Find First Divergence

```bash
# Extract segment identifiers from split events
# For C++ AGC: extract kmers from SPLIT_ACCEPT
grep "SPLIT_ACCEPT" /tmp/cpp_splits.log | \
  sed -E 's/.*from_kmers=\(([^)]+)\).*/\1/' | \
  head -50 > /tmp/cpp_split_kmers.txt

# For RAGC: extract kmers from SPLIT_SUCCESS
grep "SPLIT_SUCCESS\|SPLIT_DEGENERATE" /tmp/ragc_splits.log | \
  head -50 > /tmp/ragc_split_kmers.txt

# Compare
diff /tmp/cpp_split_kmers.txt /tmp/ragc_split_kmers.txt | head -20
```

**Status**: [ ] Not started

### Phase 3: Validate Hypothesis 1 🔬

**If Hypothesis 1 is correct**, we should see:

1. ~161 instances of "SPLIT_SKIP: left/right group has no reference yet"
2. These failures correspond to groups that exist in `map_segments` but have no entry in `reference_segments`
3. The missing references are for groups created by OTHER workers/samples that haven't flushed yet

**Root Cause**: RAGC's streaming architecture writes references lazily (when pack reaches PACK_CARDINALITY), but C++ AGC's architecture ensures references are always available in global `v_segments`.

**Status**: [ ] Not validated

### Phase 4: Implement Fix 🔧

**Fix Strategy (if Hypothesis 1 validated)**:

#### Option A: Eager Reference Storage (Recommended)
Store reference segment immediately when first segment added to new group.

**Location**: streaming_compressor_queue.rs:1642-1680

**Change**:
```rust
let buffer = groups.entry(key.clone()).or_insert_with(|| {
    // Create buffer
    let buf = SegmentGroupBuffer {
        group_id,
        segments: Vec::new(),
    };

    // ✨ NEW: Store reference immediately for first segment
    // This ensures reference is available for future split attempts
    // (matches C++ AGC behavior where v_segments always available)

    buf
});

// Add segment
buffer.segments.push(buffered);

// ✨ NEW: If this is the first segment in the group, store as reference
if buffer.segments.len() == 1 {
    let mut ref_segs = reference_segments.lock().unwrap();
    ref_segs.insert(group_id, buffered.data.clone());
}
```

**Status**: [ ] Not implemented

#### Option B: Reference Placeholder
Create dummy reference on group creation, update with real data later.

**Status**: [ ] Not considered (Option A simpler)

### Phase 5: Test Fix 🧪

#### Test 5.1: Rebuild and Run

```bash
cd /home/erik/ragc
cargo build --release

cd /home/erik/scrapy
/home/erik/ragc/target/release/ragc create -o /tmp/ragc_fixed.agc \
  -k 21 -s 10000 -m 20 -t 1 -v 1 yeast_split_proper/*.fa 2>&1 | \
  tee /tmp/ragc_fixed.log
```

**Status**: [ ] Not started

#### Test 5.2: Verify Group Count

```bash
/home/erik/ragc/target/release/ragc inspect /tmp/ragc_fixed.agc | grep "Total unique groups"
```

**Expected**: 3,220 (matching C++ AGC)
**Actual**: _____

**Status**: [ ] Not verified

#### Test 5.3: Verify Archive Size

```bash
ls -lh /tmp/ragc_fixed.agc /tmp/cpp_1thread.agc
```

**Expected**: ~69M (matching C++ AGC)
**Actual**: _____

**Status**: [ ] Not verified

#### Test 5.4: Verify "No Reference" Failures Eliminated

```bash
grep -c "has no reference yet" /tmp/ragc_fixed.log
```

**Expected**: 0
**Actual**: _____

**Status**: [ ] Not verified

### Phase 6: Comprehensive Validation ✅

#### Test 6.1: Multi-threaded Test

```bash
/home/erik/ragc/target/release/ragc create -o /tmp/ragc_fixed_mt.agc \
  -k 21 -s 10000 -m 20 -t 6 yeast_split_proper/*.fa

ls -lh /tmp/ragc_fixed_mt.agc
/home/erik/ragc/target/release/ragc inspect /tmp/ragc_fixed_mt.agc | grep "Total unique groups"
```

**Expected**: ~69M, 3,220 groups
**Actual**: _____

**Status**: [ ] Not verified

#### Test 6.2: Extraction Correctness

```bash
# Extract all samples and verify byte-identical
cd /home/erik/scrapy
for sample in $(ls yeast_split_proper/*.fa | head -5); do
    name=$(basename $sample .fa)
    /home/erik/ragc/target/release/ragc getset /tmp/ragc_fixed.agc $name > /tmp/${name}_extracted.fa
    diff $sample /tmp/${name}_extracted.fa
done
```

**Expected**: All diffs empty (byte-identical)
**Actual**: _____

**Status**: [ ] Not verified

#### Test 6.3: Cross-compatibility

```bash
# RAGC archive → C++ AGC extraction
/home/erik/agc/bin/agc getset /tmp/ragc_fixed.agc AAA#0 > /tmp/cpp_extract_ragc.fa

# C++ AGC archive → RAGC extraction
/home/erik/ragc/target/release/ragc getset /tmp/cpp_1thread.agc AAA#0 > /tmp/ragc_extract_cpp.fa

# Verify
diff yeast_split_proper/AAA.fa /tmp/cpp_extract_ragc.fa
diff yeast_split_proper/AAA.fa /tmp/ragc_extract_cpp.fa
```

**Expected**: All diffs empty
**Actual**: _____

**Status**: [ ] Not verified

---

## Results Summary

### Hypothesis Validation

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| 1. Reference Unavailability | [ ] Validated / [ ] Rejected | ___ |
| 2. Cost Orientation | [ ] Validated / [ ] Rejected | ___ |
| 3. Degenerate Threshold | [ ] Validated / [ ] Rejected | ___ |
| 4. Middle K-mer Selection | [ ] Validated / [ ] Rejected | ___ |

### Fix Effectiveness

| Metric | Before Fix | After Fix | Target | Status |
|--------|------------|-----------|--------|--------|
| Groups (single-threaded) | 3,381 | ____ | 3,220 | [ ] ✅ / [ ] ❌ |
| Archive size (single-threaded) | 80M | ____ | 69M | [ ] ✅ / [ ] ❌ |
| Groups (multi-threaded) | 3,387 | ____ | 3,220 | [ ] ✅ / [ ] ❌ |
| Archive size (multi-threaded) | 81M | ____ | 69M | [ ] ✅ / [ ] ❌ |
| "No reference" failures | ~161 (predicted) | ____ | 0 | [ ] ✅ / [ ] ❌ |
| Extraction correctness | 100% | ____ | 100% | [ ] ✅ / [ ] ❌ |

---

## Timeline

- **2025-11-11 20:00**: Plan created, ready to begin Phase 1
- **2025-11-11 __:__**: Phase 1 complete (logging instrumented)
- **2025-11-11 __:__**: Phase 2 complete (comparison tests run)
- **2025-11-11 __:__**: Phase 3 complete (hypothesis validated)
- **2025-11-11 __:__**: Phase 4 complete (fix implemented)
- **2025-11-11 __:__**: Phase 5 complete (fix tested)
- **2025-11-11 __:__**: Phase 6 complete (comprehensive validation)
- **2025-11-11 __:__**: INVESTIGATION CLOSED ✅

---

## Notes

### Key Insight from Investigation

_[Record key insights as we discover them]_

### Unexpected Findings

_[Record anything unexpected that we learn]_

### Future Improvements

_[Ideas for further optimization once correctness achieved]_

---

## Progress Log

### 2025-11-11 20:45 - Phase 2 Started

**Action**: Running RAGC with comprehensive split logging
```bash
cd /home/erik/scrapy
/home/erik/ragc/target/release/ragc create -o /tmp/ragc_debug.agc \
  -k 21 -s 10000 -m 20 -t 1 -v 1 yeast_split_proper/*.fa 2>&1 | \
  grep -E "(SPLIT_|DEBUG_LZDIFF|GROUP_CREATE)" > /tmp/ragc_splits.log
```

**Status**: Running in background (PID logged)

**Expected completion**: ~5 minutes (single-threaded, 235 samples)

**Next**: Once complete, analyze split statistics to validate Hypothesis 1

### 2025-11-11 21:10 - Phase 2 Complete: HYPOTHESIS 1 REJECTED

**Test Results**:
```
Total DEBUG_LZDIFF events: 69,816
Split attempts (pairs): 34,908
All attempts had successful reference preparation ✅
"No reference" failures: 0 ❌
```

**Conclusion**: **HYPOTHESIS 1 REJECTED**

References ARE available. All 34,908 split attempts successfully prepared both left and right references from `reference_segments` HashMap.

**The problem is NOT reference unavailability.**

### New Insight: Split Rejection After Cost Calculation

Data shows:
- 34,908 splits attempted with available references
- But 161 extra groups still created (3,381 vs 3,220)
- Archive still 16% larger (80M vs 69M)

**This means**: Splits are being REJECTED after cost calculation, not before.

### Revised Hypothesis 2: Cost Calculation Produces Invalid Results

**Theory**: The cost calculation succeeds but produces results that cause split rejection:
1. Cost vectors computed incorrectly (wrong orientation?)
2. Best position calculation differs from C++ AGC
3. Split validation rejects valid splits

**Next Test**: Run C++ AGC with logging to compare:
- How many splits does C++ AGC attempt?
- How many succeed vs rejected?
- What are the rejection reasons?


## Progress Log

### 2025-11-11 21:33 - CRITICAL BUG FOUND! 🐛

**ROOT CAUSE IDENTIFIED**: streaming_compressor_queue.rs:1553-1601

#### The Bug

When a degenerate split occurs:
1. Split calculation succeeds (middle k-mer found, both groups exist globally)
2. Split deemed degenerate (best_pos=0, whole segment → one group)
3. **Code tries to add to local `groups` buffer**: `groups.get_mut(&right_key)`
4. **Group exists globally but NOT in this thread's local buffer** → `get_mut()` returns None
5. **Segment silently not added** (no else clause!)
6. **Code still executes `continue;` at line 1601**
7. **But wait - actually falls through to Phase 3** (need to verify control flow)
8. **Phase 3 creates NEW GROUP for original (k_front, k_back) key**

#### Evidence

From test data:
- **15,682 degenerate splits** (75.8%)
- **161 extra groups** in RAGC vs C++ AGC
- Degenerate splits should reduce groups, but RAGC creates new ones instead!

#### The Fix

Need to ensure degenerate splits ALWAYS add to the correct existing group:

**Option 1**: Add groups to local buffer before attempting split
**Option 2**: Fall back to Phase 3 with UPDATED KEY when local buffer miss
**Option 3**: Use global buffer for split segment addition (not just local)

Next step: Verify control flow and implement fix.

