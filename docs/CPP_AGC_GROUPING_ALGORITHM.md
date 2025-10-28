# C++ AGC Segment Grouping Algorithm - Complete Analysis

**Date**: 2025-10-27
**Analysis of**: `/home/erik/agc/src/core/agc_compressor.cpp`

## Overview

This document describes the complete segment grouping algorithm used by C++ AGC and identifies critical differences with RAGC's implementation.

---

## C++ AGC Segment Processing Pipeline

For each segment with terminal k-mers (kmer_front, kmer_back):

### Step 1: Direct Group Lookup

```cpp
// Line 1333-1365
pair<uint64_t, uint64_t> pk = minmax(kmer1.data(), kmer2.data());
auto p = map_segments.find(pk);

if (p != map_segments.end()) {
    // Group exists! Use it directly
    segment_id = p->second;
}
```

**If group exists**: Done! Add segment to the group.

**If group doesn't exist**: Proceed to Step 2.

---

### Step 2: Segment Splitting (MISSING IN RAGC!)

**Conditions** (Line 1366-1369):
- `!concatenated_genomes` (multi-file mode)
- No group exists for (k1, k2)
- **Both k-mers are valid** (not MISSING_KMER)
- **Both k-mers exist in terminators map**

**Algorithm** (`find_cand_segment_with_missing_middle_splitter`, lines 1502-1622):

#### 2.1 Find Shared Splitters

```cpp
// Lines 1504-1525
auto p_front = map_segments_terminators.find(kmer_front.data());
auto p_back = map_segments_terminators.find(kmer_back.data());

// Find k-mers that pair with BOTH front and back
vector<uint64_t> shared_splitters;
auto p_shared = set_intersection(
    p_front->second.begin(), p_front->second.end(),
    p_back->second.begin(), p_back->second.end(),
    shared_splitters.begin()
);

// Take FIRST shared k-mer as middle splitter
uint64_t middle = shared_splitters.front();
```

**Key insight**: If k1 pairs with X and k2 pairs with X, then we can split the segment at X!

#### 2.2 Calculate Optimal Split Position

```cpp
// Lines 1532-1621
// Get existing segments for (k1, middle) and (middle, k2)
auto segment_id1 = map_segments[minmax(kmer_front, middle)];
auto segment_id2 = map_segments[minmax(middle, kmer_back)];

auto seg1 = v_segments[segment_id1];
auto seg2 = v_segments[segment_id2];

// Calculate LZ compression cost at EACH position
seg1->get_coding_cost(segment_dir, v_costs1, ...);  // Cost for left part
seg2->get_coding_cost(segment_dir, v_costs2, ...);  // Cost for right part

// Find position that minimizes total cost
uint32_t best_sum = ~0u;
uint32_t best_pos = 0;

for (uint32_t i = 0; i < v_costs1.size(); ++i) {
    uint32_t cs = v_costs1[i] + v_costs2[i];
    if (cs < best_sum) {
        best_sum = cs;
        best_pos = i;
    }
}

return make_pair(middle, best_pos);
```

#### 2.3 Split the Segment

```cpp
// Lines 1397-1454
uint32_t left_size = split_match.second;
uint32_t right_size = segment.size() - split_match.second;

if (left_size > 0 && right_size > 0) {
    // Split segment into 2 parts with OVERLAP
    uint32_t seg2_start_pos = left_size - kmer_length / 2;
    segment2.assign(segment.begin() + seg2_start_pos, segment.end());
    segment.resize(seg2_start_pos + kmer_length);

    // Part 1: (kmer_front, middle)
    pk = make_pair(kmer_front.data(), middle);
    segment_id = map_segments.at(pk);  // MUST exist

    // Part 2: (middle, kmer_back)
    pk2 = make_pair(middle, kmer_back.data());
    segment_id2 = map_segments.at(pk2);  // MUST exist

    // Store BOTH parts
}
```

**Critical**: The overlap size is `kmer_length / 2` to ensure the middle k-mer is present in both parts.

**RAGC Status**: ❌ **COMPLETELY MISSING**

---

### Step 3: One-Kmer Candidate Search

**Conditions**: Segment has only one k-mer (boundary segment).

**Algorithm** (`find_cand_segment_with_one_splitter`, lines 1625-1728):

#### 3.1 Find Candidate Groups

```cpp
// Lines 1634-1673
auto p = map_segments_terminators.find(kmer.data());
if (p == map_segments_terminators.end()) {
    // No candidates - create new MISSING_KMER group
    return make_pair((kmer.data(), ~0ull), ...);
}

// Collect ALL candidate groups
vector<tuple<k1, k2, is_rc, segment_ptr>> v_candidates;
for (auto cand_kmer : p->second) {
    pair<uint64_t, uint64_t> cand_pk = minmax(cand_kmer, kmer.data());
    v_candidates.push_back({
        cand_pk.first,
        cand_pk.second,
        is_rc,
        v_segments[map_segments[cand_pk]]  // Reference segment
    });
}
```

#### 3.2 Sort by Reference Size Similarity

```cpp
// Lines 1676-1685
int64_t segment_size = segment_dir.size();
stable_sort(v_candidates.begin(), v_candidates.end(),
    [segment_size](const auto& x, const auto& y) {
        int64_t x_size = get<3>(x)->get_ref_size();
        int64_t y_size = get<3>(y)->get_ref_size();

        // Sort by distance to segment size
        return abs(segment_size - x_size) < abs(segment_size - y_size);
    }
);
```

**Insight**: Segments of similar size compress better together (LZ matches are more likely).

#### 3.3 Test Compression with Each Candidate

```cpp
// Lines 1698-1711
uint64_t best_estim_size = segment_dir.size() - 16u;
pair<uint64_t, uint64_t> best_pk = (~0ull, ~0ull);

for (auto& candidate : v_candidates) {
    auto estim_size = get<3>(candidate)->estimate(
        get<2>(candidate) ? segment_rc : segment_dir,
        best_estim_size,
        zstd_dctx
    );

    if (estim_size < best_estim_size) {
        best_estim_size = estim_size;
        best_pk = make_pair(get<0>(candidate), get<1>(candidate));
        is_best_rc = get<2>(candidate);
    }
}

return make_pair(best_pk, is_best_rc);
```

**Key**: Uses LZ compression estimation to find the **best matching** candidate, not just the first one!

**RAGC Status**: ⚠️ **PARTIALLY IMPLEMENTED**
- RAGC tracks candidates via `group_terminators` ✓
- RAGC uses FIRST candidate instead of testing compression ✗
- RAGC doesn't sort by size similarity ✗

---

### Step 4: Fallback Minimizer Procedure

**Conditions** (Line 1461): No group found after Steps 1-3.

```cpp
if (p == map_segments.end() && fallback_filter) {
    // Use minimizer-based grouping
    // (Not analyzed in detail - appears to be a last resort)
}
```

---

## Complete C++ AGC Decision Tree

```
For each segment:
├─ Both k-mers valid?
│  ├─ YES: Group (k1, k2) exists?
│  │  ├─ YES: Add to group ✓ DONE
│  │  ├─ NO: Can split segment? (both k-mers in terminators)
│  │  │  ├─ YES: Find shared splitters
│  │  │  │  ├─ Found: Split segment at optimal position → 2 segments
│  │  │  │  └─ Not found: Try fallback minimizers
│  │  │  └─ NO: Try fallback minimizers
│  └─ NO: Only one k-mer?
│     ├─ YES: Find candidate groups via terminators
│     │  ├─ Found candidates:
│     │  │  ├─ Sort by reference size similarity
│     │  │  ├─ Test compression with each candidate
│     │  │  └─ Use best match
│     │  └─ No candidates: Create MISSING_KMER group
│     └─ NO: Neither k-mer valid (error case)
```

---

## RAGC vs C++ AGC Comparison

| Feature | C++ AGC | RAGC | Impact |
|---------|---------|------|--------|
| **Direct group lookup** | ✓ | ✓ | None |
| **Segment splitting** | ✓ Lines 1366-1456 | ❌ **MISSING** | **CRITICAL** |
| **One-kmer terminators tracking** | ✓ | ✓ (commit 1f2b580) | Fixed |
| **One-kmer compression testing** | ✓ All candidates | ❌ First candidate only | Moderate |
| **One-kmer size-based sorting** | ✓ | ❌ | Minor |
| **Fallback minimizers** | ✓ | ? | Unknown |

---

## Root Cause of Group Explosion

**Problem**: RAGC creates 803 groups vs ~240 optimal for yeast10.

**Analysis**:
- 5206 segments ÷ 50 (PACK_CARDINALITY) = ~104 ideal packs
- 803 groups = 6.5 segments/pack (should be 50)
- This means segments are creating **new unique (k1, k2) pairs** instead of joining existing groups

**Root Cause**: Missing segment splitting feature!

### Example Scenario

Imagine we have these existing groups:
- Group A: (k1, k_middle)
- Group B: (k_middle, k2)

When a segment arrives with (k1, k2):

**C++ AGC**:
1. No group for (k1, k2) exists
2. Finds k_middle as shared splitter
3. Splits segment into two parts:
   - Part 1 → Group A: (k1, k_middle)
   - Part 2 → Group B: (k_middle, k2)
4. **Result**: No new group created! ✓

**RAGC**:
1. No group for (k1, k2) exists
2. ❌ Doesn't try to split
3. Creates **new group** for (k1, k2)
4. **Result**: Unnecessary new group created! ✗

### Impact

For a dataset with 240 splitters creating a network of k-mer pairs:
- **C++ AGC**: Segments split to match existing groups → ~240 groups
- **RAGC**: Each unique (k1, k2) pair creates a new group → 803 groups

The 3.3x group explosion comes from segments that **could be split** but aren't.

---

## Compression Quality Impact

### Group Count → Pack Count → Compression Ratio

```
More groups → Fewer segments/pack → Less LZ compression → Larger files
```

**Measurements**:

| Implementation | Groups | Packs | Segs/Pack | File Size | vs C++ AGC |
|----------------|--------|-------|-----------|-----------|------------|
| **C++ AGC** | ~240 | ~104 | ~50 | 9.6M | — |
| **RAGC (before fix)** | 1537 | 1537 | 3.4 | 23M | +140% |
| **RAGC (one-kmer fix)** | 803 | 803 | 6.5 | 15M | +56% |
| **RAGC (with splitting)** | ~240? | ~104? | ~50? | ~9.6M? | ~0% |

**Key insight**: The 56% file size gap is NOT a compression parameter issue - it's a **grouping algorithm gap**.

---

## Implementation Priority

### Critical: Segment Splitting

**Without this feature**, RAGC cannot match C++ AGC's compression quality.

**Implementation requirements**:
1. Detect when (k1, k2) group doesn't exist but both k-mers are valid
2. Find shared splitters using set intersection on terminators
3. Calculate LZ compression cost at each position
4. Find optimal split position
5. Split segment with overlap (kmer_length/2)
6. Store as two separate segments

**Complexity**: High
- Requires LZ compression cost calculation (`get_coding_cost`)
- Requires segment splitting with overlap
- Requires storing multiple segments per original segment
- May impact overall architecture

### Important: Compression-Based Candidate Selection

**Current RAGC**: Uses first candidate
**C++ AGC**: Tests all candidates, selects best compression

**Impact**: Moderate (only affects boundary segments)

**Implementation requirements**:
1. Sort candidates by reference size similarity
2. Implement `estimate()` function for LZ compression estimation
3. Test each candidate, track best

**Complexity**: Medium
- Need LZ compression estimation (already exists in RAGC)
- Need to iterate all candidates (may be slow)
- Could optimize by testing only top N candidates

---

## Testing Strategy

### Verify Segment Splitting

After implementing splitting:

```bash
# Should match C++ AGC exactly
/home/erik/ragc/target/release/ragc create -o ragc_split.agc -k 31 -v 2 samples/*.fa
/home/erik/agc/bin/agc create -o cpp.agc -k 31 samples/*.fa

# Compare
ls -lh ragc_split.agc cpp.agc

# Should see:
# - Similar group count (~240)
# - Similar pack count (~104)
# - Similar file size (~9.6M)
```

### Debug Output

Add logging to show:
- When splitting is attempted
- Shared splitters found
- Split position chosen
- Resulting segment pairs

---

## Questions for User

1. **Priority**: Should we implement segment splitting immediately, or explore other optimizations first?

2. **Architecture**: Segment splitting requires storing multiple segments per contig. Should we:
   - Modify existing segment storage?
   - Create a separate "split segments" collection?
   - Restructure the entire grouping phase?

3. **Testing**: Do we have any datasets where C++ AGC and RAGC produce identical file sizes? This would help verify our understanding.

4. **Fallback minimizers**: Should we investigate the fallback minimizer procedure as well, or is segment splitting sufficient?

---

## Summary

**C++ AGC uses a sophisticated 4-tier grouping strategy**:
1. Direct group lookup (fast path)
2. **Segment splitting** when (k1, k2) doesn't exist ← **MISSING IN RAGC**
3. One-kmer candidate search with compression testing ← **SIMPLIFIED IN RAGC**
4. Fallback minimizer grouping ← **UNKNOWN IN RAGC**

**The 56% file size gap is explained by**:
- Missing segment splitting: Creates 3.3x more groups (803 vs ~240)
- Simplified one-kmer matching: Uses first candidate instead of best compression

**To achieve C++ AGC parity**: Segment splitting is **essential**, compression testing is **beneficial**.
