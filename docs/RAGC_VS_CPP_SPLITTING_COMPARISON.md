# RAGC vs C++ AGC Splitting Implementation Comparison

## Executive Summary

This document compares RAGC's segment splitting implementation against C++ AGC to identify differences that may explain the group fragmentation issue (RAGC creates 1713 groups vs C++ AGC's 1589 groups, +7.8% fragmentation).

**Status**: RAGC implementation closely matches C++ AGC algorithm with one critical difference in the LZDiff fallback path.

## Side-by-Side Comparison

### Step 1: Splitting Conditions

**C++ AGC** (`agc_compressor.cpp:1394-1398`):
```cpp
if (!concatenated_genomes &&
    p == map_segments.end() &&
    pk.first != ~0ull && pk.second != ~0ull &&
    map_segments_terminators.count(pk.first) &&
    map_segments_terminators.count(pk.second))
```

**RAGC** (`streaming_compressor_queue.rs:1292`):
```rust
if !key_exists && key_front != MISSING_KMER && key_back != MISSING_KMER {
```

**Analysis**:
- ✅ Both check segment is new (`!key_exists` / `p == end()`)
- ✅ Both check k-mers are valid (not MISSING)
- ⚠️ **DIFFERENCE**: C++ AGC explicitly checks k-mers exist in terminators
- ⚠️ **DIFFERENCE**: RAGC checks terminators INSIDE `find_middle_splitter()` (returns None if not found)

**Verdict**: Functionally equivalent. RAGC defers terminator check to `find_middle_splitter()`.

### Step 2: Find Middle Splitter

**C++ AGC** (`agc_compressor.cpp:1531-1554`):
```cpp
auto p_front = map_segments_terminators.find(kmer_front.data());
auto p_back = map_segments_terminators.find(kmer_back.data());

if (p_front == map_segments_terminators.end() ||
    p_back == map_segments_terminators.end())
    return make_pair(~0ull, 0);

// Find intersection
auto p_shared = set_intersection(
    p_front->second.begin(), p_front->second.end(),
    p_back->second.begin(), p_back->second.end(),
    shared_splitters.begin());

shared_splitters.erase(remove(shared_splitters.begin(),
    shared_splitters.end(), ~0ull), shared_splitters.end());

uint64_t middle = shared_splitters.front(); // Take 1st shared k-mer
```

**RAGC** (`streaming_compressor_queue.rs:1533-1584`):
```rust
let front_connections = terminators.get(&front_kmer)?;
let back_connections = terminators.get(&back_kmer)?;

// Find intersection (manual set_intersection implementation)
let mut i = 0;
let mut j = 0;
let mut shared_kmers = Vec::new();

while i < front_connections.len() && j < back_connections.len() {
    if front_connections[i] == back_connections[j] {
        shared_kmers.push(front_connections[i]);
        i += 1;
        j += 1;
    } else if front_connections[i] < back_connections[j] {
        i += 1;
    } else {
        j += 1;
    }
}

// Filter out MISSING_KMER
shared_kmers.retain(|&k| k != MISSING_KMER);

shared_kmers.first().copied()
```

**Analysis**:
- ✅ Both return None/empty if terminators not found
- ✅ Both use set intersection algorithm
- ✅ Both filter out MISSING_KMER
- ✅ Both take FIRST shared k-mer
- ✅ **IDENTICAL LOGIC**

**Verdict**: Perfect match.

### Step 3: Check Split Groups Exist

**C++ AGC** (`agc_compressor.cpp:1561-1562`):
```cpp
auto segment_id1 = map_segments[minmax(kmer_front.data(), middle)];
auto segment_id2 = map_segments[minmax(middle, kmer_back.data())];
```

**RAGC** (`streaming_compressor_queue.rs:1329-1342`):
```rust
let both_groups_exist = {
    let map = map_segments.lock().unwrap();
    let left_exists = map.contains_key(&left_key);
    let right_exists = map.contains_key(&right_key);

    left_exists && right_exists
};

if both_groups_exist {
    // Proceed with split
}
```

**Analysis**:
- ❌ **CRITICAL DIFFERENCE**: C++ uses `operator[]` which ASSUMES keys exist
- ✅ **RAGC IS MORE DEFENSIVE**: Explicitly checks both groups exist
- ⚠️ **CONCERN**: C++ AGC might crash if groups don't exist (using `.at()` later at lines 1468, 1482)

**Verdict**: RAGC is safer, but this could be a source of divergence if C++ AGC's assumption is wrong.

### Step 4: Calculate Compression Costs

**C++ AGC** (`agc_compressor.cpp:1567-1630`):
```cpp
auto seg1 = v_segments[segment_id1];
auto seg2 = v_segments[segment_id2];

// Calculate costs
seg1->get_coding_cost(segment_dir, v_costs1, true, zstd_dctx);
partial_sum(v_costs1.begin(), v_costs1.end(), v_costs1.begin());

seg2->get_coding_cost(segment_dir, v_costs2, false, nullptr);
partial_sum(v_costs2.rbegin(), v_costs2.rend(), v_costs2.rbegin());
```

**RAGC** (`streaming_compressor_queue.rs:1650-1706`):
```rust
// Get left and right group buffers
let left_buffer = groups.get(left_key)?;
let right_buffer = groups.get(right_key)?;

// Check if both groups have LZDiff prepared
let (left_lz, right_lz) = match (&left_buffer.lz_diff, &right_buffer.lz_diff) {
    (Some(left), Some(right)) => (left, right),
    _ => {
        // 🚨 FALLBACK TO MIDPOINT 🚨
        let best_pos = segment_len / 2;
        // ... validate and split at midpoint ...
        return Some((left_data, right_data, middle_kmer));
    }
};

// Calculate compression costs
let v_costs1 = left_lz.get_coding_cost_vector(segment_data, true);
let v_costs2 = right_lz.get_coding_cost_vector(segment_data, false);
```

**Analysis**:
- ❌ **CRITICAL DIFFERENCE**: RAGC has midpoint fallback when LZDiff not ready
- ❌ **C++ AGC ASSUMPTION**: Reference is ALWAYS written before split attempt
- ⚠️ **POTENTIAL BUG**: RAGC's fallback could split at suboptimal positions

**Question**: Does C++ AGC guarantee that both split groups have references written before attempting split?

**Investigation needed**: Check when `lz_diff` is prepared in RAGC vs when C++ AGC writes references.

### Step 5: Find Best Split Position

**C++ AGC** (`agc_compressor.cpp:1632-1643`):
```cpp
uint32_t best_sum = ~0u;
uint32_t best_pos = 0;

for (uint32_t i = 0; i < v_costs1.size(); ++i) {
    uint32_t cs = v_costs1[i] + v_costs2[i];
    if (cs < best_sum) {
        best_sum = cs;
        best_pos = i;
    }
}
```

**RAGC** (`streaming_compressor_queue.rs:1723-1732`):
```rust
let mut best_sum = u32::MAX;
let mut best_pos = 0usize;

for i in 0..v_costs1.len() {
    let cs = v_costs1[i].saturating_add(v_costs2[i]);
    if cs < best_sum {
        best_sum = cs;
        best_pos = i;
    }
}
```

**Analysis**:
- ✅ Both iterate all positions
- ✅ Both find minimum combined cost
- ✅ **IDENTICAL LOGIC** (saturating_add is just safe arithmetic)

**Verdict**: Perfect match.

### Step 6: Validate Split Position

**C++ AGC** (`agc_compressor.cpp:1645-1648`):
```cpp
if (best_pos < kmer_length + 1u)
    best_pos = 0;
if ((size_t)best_pos + kmer_length + 1u > v_costs1.size())
    best_pos = (uint32_t)v_costs1.size();
```

**RAGC** (`streaming_compressor_queue.rs:1734-1742`):
```rust
let k = config.k;
if best_pos < k + 1 {
    best_pos = 0; // Too close to start
}
if best_pos + k + 1 > v_costs1.len() {
    best_pos = v_costs1.len(); // Too close to end
}
```

**Analysis**:
- ✅ Both check minimum distance from start (k+1)
- ✅ Both check minimum distance from end (k+1)
- ✅ Both adjust to boundaries if too close
- ✅ **IDENTICAL LOGIC**

**Verdict**: Perfect match.

### Step 7: Check for Degenerate Split

**C++ AGC** (`agc_compressor.cpp:1442-1459` - happens AFTER validation):
```cpp
// If split position is at boundaries, don't split
if (best_pos == 0 || best_pos >= segment.size())
    // Don't split - add as singleton
```

**RAGC** (`streaming_compressor_queue.rs:1744-1755`):
```rust
if best_pos == 0 || best_pos >= segment_data.len() {
    // Degenerate split - position at or beyond segment boundaries
    return None;
}
```

**Analysis**:
- ✅ Both check if position is at boundaries
- ✅ Both skip split if degenerate
- ✅ **IDENTICAL LOGIC**

**Verdict**: Perfect match.

### Step 8: Split Segment with Overlap

**C++ AGC** (`agc_compressor.cpp:1450-1454`):
```cpp
uint32_t seg2_start_pos = left_size - kmer_length / 2;
segment2.assign(segment.begin() + seg2_start_pos, segment.end());
segment.resize((size_t)seg2_start_pos + kmer_length);
```

**RAGC** (`streaming_compressor_queue.rs:1607-1627`):
```rust
fn split_segment_at_position(
    segment_data: &[u8],
    split_pos: usize,
    k: usize,
) -> (Vec<u8>, Vec<u8>) {
    let seg2_start_pos = split_pos.saturating_sub(k / 2);

    // Right segment: [seg2_start_pos .. end]
    let right = segment_data[seg2_start_pos..].to_vec();

    // Left segment: [0 .. seg2_start_pos + k]
    let left_end = seg2_start_pos + k;
    let left = segment_data[..left_end].to_vec();

    (left, right)
}
```

**Analysis**:
- ✅ Both calculate `seg2_start_pos = split_pos - k/2`
- ✅ Both create right segment from `[seg2_start_pos..end]`
- ✅ Both create left segment from `[0..seg2_start_pos + k]`
- ✅ Both create k-byte overlap
- ✅ **IDENTICAL LOGIC**

**Verdict**: Perfect match.

## Critical Differences Summary

### 1. LZDiff Fallback (MAJOR)

**Location**: `streaming_compressor_queue.rs:1658-1699`

**Issue**: RAGC falls back to midpoint splitting when LZDiff is not ready.

**C++ AGC behavior**: ASSUMES reference is always written before split attempt.

**Impact**: Could cause suboptimal splits if references aren't ready, leading to different segment boundaries and group assignments.

**Investigation needed**:
- When is `lz_diff` set in RAGC vs when C++ AGC writes references?
- Does C++ AGC guarantee references are written before attempting splits?
- Are there cases where RAGC attempts splits before references are ready?

### 2. Defensive Group Existence Check (MINOR)

**Location**: `streaming_compressor_queue.rs:1329-1342`

**Issue**: RAGC explicitly checks both split groups exist before attempting split.

**C++ AGC behavior**: Uses `operator[]` which assumes groups exist (could insert defaults if not found).

**Impact**: Probably none - if C++ AGC's assumption is correct, this is just safer programming.

## Hypotheses for Group Fragmentation

Based on this analysis, here are potential root causes for RAGC's +124 extra groups:

### Hypothesis 1: LZDiff Not Ready (MOST LIKELY)

**Theory**: RAGC's midpoint fallback creates different segment boundaries than C++ AGC's cost-based splits.

**Evidence**:
- RAGC has 444 fewer segments (11,208 vs 11,647)
- Sample AIF#2 chrI: 15 segments (RAGC) vs 16 (C++)
- Midpoint splits could create different terminal k-mers

**Test**: Add logging to see how often midpoint fallback is used.

**Fix**: Ensure LZDiff is written BEFORE attempting splits (match C++ AGC's processing order).

### Hypothesis 2: Processing Order Difference

**Theory**: Priority queue changes when segments are processed, affecting which groups exist when split is attempted.

**Evidence**:
- Priority queue implementation is correct
- But processing order could affect which groups are created first
- If segment A creates group X, then segment B tries to split into X+Y, but Y doesn't exist yet...

**Test**: Log the order segments are processed in RAGC vs C++ AGC.

**Fix**: Ensure references are written before segments that might split into those groups.

### Hypothesis 3: Terminator Map Population Timing

**Theory**: `map_segments_terminators` might not be fully populated when splits are attempted.

**Evidence**:
- RAGC populates terminators during segment processing
- C++ AGC might populate terminators in a different phase

**Test**: Verify terminator map contents match C++ AGC at split time.

**Fix**: Ensure all terminal k-mers are registered before attempting splits.

## Next Steps

1. **Add verbose logging** to RAGC splitting logic:
   - Log when midpoint fallback is used
   - Log which splits succeed vs fail (groups not ready)
   - Log segment processing order

2. **Compare with C++ AGC**:
   - Run both with same dataset
   - Compare split decisions for specific segments
   - Identify first divergence point

3. **Investigate LZDiff preparation timing**:
   - When is `lz_diff` set in RAGC?
   - When does C++ AGC write references?
   - Are there ordering guarantees?

4. **Fix root cause**:
   - If LZDiff timing issue: Ensure references written before splits attempted
   - If processing order issue: Match C++ AGC's segment processing order
   - If terminator timing issue: Pre-populate terminator map

## Testing Strategy

To identify the root cause:

```bash
# Test with verbose logging
cd /home/erik/scrapy
/home/erik/ragc/target/release/ragc create \
    -o /tmp/test_verbose.agc \
    -k 21 -s 10000 -m 20 -v 2 \
    yeast_split_proper/*.fa 2>&1 | tee /tmp/ragc_verbose.log

# Search for fallback usage
grep "SPLIT_FALLBACK" /tmp/ragc_verbose.log | wc -l

# Search for split failures
grep "SPLIT_SKIP" /tmp/ragc_verbose.log | wc -l

# Search for successful splits
grep "SPLIT_SUCCESS" /tmp/ragc_verbose.log | wc -l
```

Compare segment counts for specific samples:

```bash
# Check AIF#2 chrI (known to have 15 vs 16 segments)
/home/erik/ragc/target/release/ragc inspect /tmp/test_verbose.agc --segments | \
    grep -A 20 "AIF#2/AIF#2#chrI"
```

## Conclusion

RAGC's splitting implementation is **very close** to C++ AGC, with one critical difference: the **midpoint fallback when LZDiff is not ready**. This is the most likely cause of the group fragmentation issue.

**Recommended fix**: Ensure `lz_diff` is prepared (reference written) BEFORE attempting splits on segments that could use those groups. This may require adjusting the processing order or adding a two-phase approach (write references first, then attempt splits).
