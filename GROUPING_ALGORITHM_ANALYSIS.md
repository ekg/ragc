# Mathematical Analysis: RAGC vs C++ AGC Grouping Algorithm

**Date**: 2025-11-11
**Problem**: RAGC creates 161 extra groups (3,381 vs 3,220) causing 16% archive bloat

---

## Executive Summary

Both implementations use the same **high-level algorithm** but differ in **when and how** they attempt segment splitting. The key difference is in the **splitting preconditions** that allow a segment to be split into existing groups.

---

## Algorithm Overview: Segment Grouping with Dynamic Splitting

### Input
- Segment S with boundary k-mers: `(k_front, k_back)`
- Global registry `map_segments`: `(k1, k2) → group_id`
- Terminator registry `map_segments_terminators`: `k → [k_shared1, k_shared2, ...]`

### Decision Tree

For each segment S:

```
1. Check if group (k_front, k_back) exists in map_segments
   YES → Add S to existing group (DONE)
   NO  → Continue to step 2

2. Try to split S using middle k-mer:
   2a. Find middle k-mer k_mid:
       k_mid ∈ terminators[k_front] ∩ terminators[k_back]

   2b. Check if BOTH split groups exist:
       left_group  = (min(k_front, k_mid), max(k_front, k_mid))
       right_group = (min(k_mid, k_back), max(k_mid, k_back))

       IF left_group ∈ map_segments AND right_group ∈ map_segments:
           → Attempt cost-based splitting
       ELSE:
           → Continue to step 3

   2c. Calculate split cost:
       For each position i in segment:
           cost[i] = LZ_diff_cost(left_ref, S[0..i]) + LZ_diff_cost(right_ref, S[i..end])

       best_pos = argmin(cost[i])

   2d. Validate split position:
       IF best_pos creates valid segments (size ≥ k+1):
           → Split S at best_pos
           → Add left part to left_group
           → Add right part to right_group
           → DONE
       ELSE:
           → Continue to step 3

3. Create new group for (k_front, k_back):
   group_id = next_available_id()
   map_segments[(k_front, k_back)] = group_id
   Add S to new group
```

---

## Critical Difference: Split Attempt Timing

### C++ AGC Logic (agc_compressor.cpp:1380-1484)

```cpp
// Line 1380-1384: Split condition
if (!concatenated_genomes &&
    p == map_segments.end() &&                    // Group doesn't exist
    pk.first != ~0ull && pk.second != ~0ull &&    // Both k-mers present
    map_segments_terminators.count(pk.first) &&   // Terminators exist
    map_segments_terminators.count(pk.second))
{
    // Line 1408: Find middle k-mer and attempt split
    auto split_match = find_cand_segment_with_missing_middle_splitter(...);

    if (split_match.first != ~0ull)
    {
        // Line 1454, 1468: Verify split groups MUST exist
        segment_id = map_segments.at(pk);          // must exists (CRASH if not)
        segment_id2 = map_segments.at(pk2);        // must exists (CRASH if not)

        // Add to existing groups
        buffered_seg_part.add_known(...);
    }
}

// Line 1506-1508: If split failed or impossible, create new group
if (p == map_segments.end())
{
    buffered_seg_part.add_new(...);
}
```

**Key point**: C++ AGC uses `.at()` which **asserts** that both split groups exist. If they don't, it crashes. This means `find_cand_segment_with_missing_middle_splitter` MUST verify group existence internally.

### RAGC Logic (streaming_compressor_queue.rs:1437-1615)

```rust
// Phase 1: Check if group exists
let key_exists = {
    let seg_map = map_segments.lock().unwrap();
    seg_map.contains_key(&key)
};

// Phase 2: Try to split into existing groups
// Lines 1447-1615
if key_front != MISSING_KMER && key_back != MISSING_KMER {
    // Find middle k-mer
    let middle_kmer_opt = {
        let terminators = map_segments_terminators.lock().unwrap();
        find_middle_splitter(key_front, key_back, &terminators)
    };

    if let Some(middle_kmer) = middle_kmer_opt {
        // Check if BOTH split groups exist
        let both_groups_exist = {
            let map = map_segments.lock().unwrap();
            let left_exists = map.contains_key(&left_key);
            let right_exists = map.contains_key(&right_key);
            left_exists && right_exists
        };

        if both_groups_exist {
            // Attempt cost-based split
            let split_result = try_split_segment_with_cost(...);

            if let Some((left_data, right_data, _mid)) = split_result {
                // Add to existing groups
                ...
                continue;  // Skip Phase 3
            }
        }
    }
}

// Phase 3: Create new group if split failed/impossible
let group_id = {
    let mut seg_map = map_segments.lock().unwrap();
    *seg_map.entry(key.clone()).or_insert_with(|| {
        group_counter.fetch_add(1, Ordering::SeqCst)
    })
};
```

**Key point**: RAGC explicitly checks `both_groups_exist` before attempting split. If either group is missing, it falls through to create a new group.

---

## The Mathematical Problem

### Definition: Group Creation Probability

Let's define:
- `G(t)` = set of groups that exist at time t
- `S` = segment with k-mers `(k_front, k_back)`
- `m` = middle k-mer such that `m ∈ terminators[k_front] ∩ terminators[k_back]`

**Splitting succeeds if and only if**:
```
(k_front, m) ∈ G(t) AND (m, k_back) ∈ G(t)
```

**If split fails, new group created**:
```
G(t+1) = G(t) ∪ {(k_front, k_back)}
```

### Problem: Dependency Ordering

Consider three samples processed sequentially:

**Sample 1**: Creates groups G₁ = {(A,B), (B,C), ...}

**Sample 2**: Encounters segment with k-mers (A,C)
- Middle k-mer found: m = B
- Split attempt:
  - Left group: (A,B) ∈ G₁? **YES**
  - Right group: (B,C) ∈ G₁? **YES**
  - **Split succeeds** → Adds to existing groups
  - **No new group created**

**Sample 3**: Encounters segment with k-mers (A,D)
- Middle k-mer found: m = B
- Split attempt:
  - Left group: (A,B) ∈ G₁? **YES**
  - Right group: (B,D) ∈ G₁ ∪ G₂? **DEPENDS**

**If (B,D) was created in Sample 2 → Split succeeds**
**If (B,D) doesn't exist yet → New group (A,D) created**

### The 161 Extra Groups

**Hypothesis**: RAGC is too conservative in split attempts, creating groups that C++ AGC would avoid by splitting.

**Possible causes**:
1. **Missing reference segments**: `try_split_segment_with_cost` returns `None` if reference data unavailable
2. **Cost calculation fails**: `v_costs1` or `v_costs2` empty/mismatched
3. **Degenerate splits rejected**: Split position creates segments < k+1 bytes
4. **Middle k-mer selection differs**: Different choice from intersection set

---

## Next Steps: Systematic Comparison

### Test 1: Log all split attempts

Add detailed logging to both implementations:
- When is middle k-mer found?
- When do split groups both exist?
- When does cost calculation succeed/fail?
- What is the split position chosen?

### Test 2: Compare split attempt statistics

Run both on same dataset:
```bash
# C++ AGC
agc create -o cpp.agc -k 21 -s 10000 -t 1 samples/*.fa 2>&1 | grep SPLIT > cpp_splits.txt

# RAGC
ragc create -o ragc.agc -k 21 -s 10000 -t 1 samples/*.fa 2>&1 | grep SPLIT > ragc_splits.txt

# Compare
diff cpp_splits.txt ragc_splits.txt
```

### Test 3: Identify first divergence

Find the first segment where C++ AGC splits but RAGC creates a new group:
```bash
diff <(grep "SPLIT_ACCEPT" cpp_splits.txt) <(grep "SPLIT_SUCCESS" ragc_splits.txt) | head -1
```

Then investigate WHY that split failed in RAGC but succeeded in C++ AGC.

---

## Specific Areas to Investigate

### Area 1: Reference Segment Availability

**C++ AGC** (agc_compressor.cpp:1555-1556):
```cpp
auto segment_id1 = map_segments[minmax(kmer_front.data(), middle)];
auto segment_id2 = map_segments[minmax(middle, kmer_back.data())];
auto seg1 = v_segments[segment_id1];  // Global vector, always available
auto seg2 = v_segments[segment_id2];
```

**RAGC** (streaming_compressor_queue.rs:1906-1942):
```rust
let prepare_on_demand = |key: &SegmentGroupKey, label: &str| -> Option<LZDiff> {
    let ref_segments_locked = reference_segments.lock().unwrap();

    if let Some(&segment_id) = map_segments_locked.get(key) {
        if let Some(ref_data) = ref_segments_locked.get(&segment_id) {
            // Reference exists!
            ...
        } else {
            // ⚠️ Segment ID exists but no reference data!
            return None;
        }
    }
    None
};
```

**Critical difference**: RAGC might fail split if reference segment hasn't been written yet, even if group ID exists in `map_segments`.

### Area 2: Cost Calculation Orientation

**C++ AGC** (agc_compressor.cpp:1562-1568):
```cpp
if (kmer_front.data() < middle)
    seg1->get_coding_cost(segment_dir, v_costs1, true, zstd_dctx);
else
{
    seg1->get_coding_cost(segment_rc, v_costs1, false, zstd_dctx);
    reverse(v_costs1.begin(), v_costs1.end());
}
```

**RAGC** (streaming_compressor_queue.rs:1947-1965):
```rust
let v_costs1 = if let Some(lz) = prepare_on_demand(left_key, "left") {
    lz.get_coding_cost_vector(segment_data, true)  // Always true for left?
} else {
    return None;
};

let v_costs2 = if let Some(lz) = prepare_on_demand(right_key, "right") {
    lz.get_coding_cost_vector(segment_data, false)  // Always false for right?
} else {
    return None;
};
```

**Question**: Does RAGC handle orientation correctly when `key_front > middle` or `middle > key_back`?

---

## Expected Outcome

Once we identify the specific condition causing splits to fail in RAGC, we can:
1. Fix the condition
2. Verify group count drops to 3,220 (matching C++ AGC)
3. Verify archive size drops to 69M (matching C++ AGC)
4. Confirm 100% byte-identical extraction

The bloat is algorithmic, not concurrent - this is a logic bug, not a race condition.
