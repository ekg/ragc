# C++ AGC Segment Splitting Algorithm

## Overview

C++ AGC splits segments with missing middle k-mers to improve compression by merging them with existing groups. This document comprehensively explains the splitting algorithm based on analysis of `agc_compressor.cpp`.

## When Splitting is Attempted

**Location**: `agc_compressor.cpp` lines 1394-1398

**Conditions** (ALL must be true):
1. NOT in `concatenated_genomes` mode
2. Group key NOT already in `map_segments` (`p == map_segments.end()`)
3. Front k-mer is valid (not `~0ull` / MISSING_KMER)
4. Back k-mer is valid (not `~0ull` / MISSING_KMER)
5. Front k-mer exists in `map_segments_terminators`
6. Back k-mer exists in `map_segments_terminators`

```cpp
if (!concatenated_genomes &&
    p == map_segments.end() &&
    pk.first != ~0ull && pk.second != ~0ull &&
    map_segments_terminators.count(pk.first) && map_segments_terminators.count(pk.second))
```

**Critical insight**: Splitting only happens for NEW segments (not in map_segments yet).

## Step 1: Find Middle Splitter

**Location**: `agc_compressor.cpp` lines 1531-1554

**Function**: `find_cand_segment_with_missing_middle_splitter()`

**Algorithm**:
1. Look up front k-mer in `map_segments_terminators` → get vector of connected k-mers
2. Look up back k-mer in `map_segments_terminators` → get vector of connected k-mers
3. Find intersection using `set_intersection()` (both vectors are sorted)
4. Remove `~0ull` (MISSING_KMER) from intersection
5. If intersection is empty → return `(~0ull, 0)` (no middle splitter found)
6. Otherwise, take FIRST shared k-mer as middle splitter

```cpp
auto p_front = map_segments_terminators.find(kmer_front.data());
auto p_back = map_segments_terminators.find(kmer_back.data());

// Find intersection
auto p_shared = set_intersection(
    p_front->second.begin(), p_front->second.end(),
    p_back->second.begin(), p_back->second.end(),
    shared_splitters.begin());

// Take first shared k-mer
uint64_t middle = shared_splitters.front();
```

**Example**:
- Segment has front k-mer `AAA` and back k-mer `TTT`
- `map_segments_terminators[AAA]` = `[GGG, CCC, TTT, ...]`
- `map_segments_terminators[TTT]` = `[AAA, GGG, CCC, ...]`
- Intersection = `[GGG, CCC]`
- Middle splitter = `GGG` (first one)

## Step 2: Get Existing Segment Groups

**Location**: `agc_compressor.cpp` lines 1561-1565

**Algorithm**:
1. Create left group key: `minmax(kmer_front, middle)`
2. Create right group key: `minmax(middle, kmer_back)`
3. Look up left segment ID: `map_segments[left_key]`
4. Look up right segment ID: `map_segments[right_key]`
5. Get actual segment objects from `v_segments`

```cpp
auto segment_id1 = map_segments[minmax(kmer_front.data(), middle)];
auto segment_id2 = map_segments[minmax(middle, kmer_back.data())];

auto seg1 = v_segments[segment_id1];
auto seg2 = v_segments[segment_id2];
```

**CRITICAL**: Uses `operator[]` which assumes keys exist. C++ AGC expects both split groups to already exist in `map_segments` before attempting split.

## Step 3: Calculate Compression Costs

**Location**: `agc_compressor.cpp` lines 1567-1630

**Algorithm**:
1. Calculate prefix costs using left segment's `get_coding_cost()`:
   - For each position in segment, calculate LZ diff cost
   - Apply `partial_sum()` to get cumulative prefix costs
   - Result: `v_costs1[i]` = cost of compressing `segment[0..i]` against left reference

2. Calculate suffix costs using right segment's `get_coding_cost()`:
   - For each position in segment, calculate LZ diff cost
   - Apply `partial_sum()` in reverse to get cumulative suffix costs
   - Result: `v_costs2[i]` = cost of compressing `segment[i..end]` against right reference

```cpp
// Left segment costs (prefix)
if (kmer_front.data() < middle)
    seg1->get_coding_cost(segment_dir, v_costs1, true, zstd_dctx);
else {
    seg1->get_coding_cost(segment_rc, v_costs1, false, zstd_dctx);
    reverse(v_costs1.begin(), v_costs1.end());
}
partial_sum(v_costs1.begin(), v_costs1.end(), v_costs1.begin());

// Right segment costs (suffix)
if (middle < kmer_back.data()) {
    seg2->get_coding_cost(segment_dir, v_costs2, false, nullptr);
    partial_sum(v_costs2.rbegin(), v_costs2.rend(), v_costs2.rbegin());
} else {
    seg2->get_coding_cost(segment_rc, v_costs2, true, nullptr);
    partial_sum(v_costs2.begin(), v_costs2.end(), v_costs2.begin());
    reverse(v_costs2.begin(), v_costs2.end());
}
```

**Note**: The orientation (forward vs reverse-complement) depends on k-mer ordering to ensure proper alignment with reference segments.

## Step 4: Find Best Split Position

**Location**: `agc_compressor.cpp` lines 1632-1643

**Algorithm**:
1. For each position `i` from 0 to segment length:
   - Calculate combined cost: `cs = v_costs1[i] + v_costs2[i]`
   - Track position with minimum combined cost
2. Return position with best (lowest) cost

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

## Step 5: Validate Split Position

**Location**: `agc_compressor.cpp` lines 1645-1648

**Algorithm**:
1. If `best_pos < kmer_length + 1` → set `best_pos = 0` (too close to start)
2. If `best_pos + kmer_length + 1 > segment.size()` → set `best_pos = segment.size()` (too close to end)

```cpp
if (best_pos < kmer_length + 1u)
    best_pos = 0;
if ((size_t)best_pos + kmer_length + 1u > v_costs1.size())
    best_pos = (uint32_t)v_costs1.size();
```

**Degenerate Check** (lines 1442-1459): After validation, if `best_pos` is at boundaries (0 or segment.size()), the split is considered degenerate and NOT performed.

## Step 6: Split Segment with K-mer Overlap

**Location**: `agc_compressor.cpp` lines 1450-1454

**Algorithm**:
1. Calculate `seg2_start_pos = best_pos - kmer_length / 2`
2. Create right segment: `segment2 = segment[seg2_start_pos..end]`
3. Create left segment: `segment = segment[0..seg2_start_pos + kmer_length]`

```cpp
// Split segment into 2 parts (with overlap of size kmer_length)
uint32_t seg2_start_pos = left_size - kmer_length / 2;
segment2.assign(segment.begin() + seg2_start_pos, segment.end());
segment.resize((size_t)seg2_start_pos + kmer_length);
```

**Overlap Calculation**:
- Left segment ends at: `seg2_start_pos + k = (best_pos - k/2) + k = best_pos + k/2`
- Right segment starts at: `best_pos - k/2`
- Overlap region: `[best_pos - k/2 .. best_pos + k/2]` = **k bytes**

**Example** (k=21, best_pos=100):
- `seg2_start_pos = 100 - 10 = 90`
- Left: `[0..90+21] = [0..111]` (111 bytes)
- Right: `[90..end]` (end - 90 bytes)
- Overlap: `[90..111]` = 21 bytes ✓

## Step 7: Add Split Segments to Groups

**Location**: `agc_compressor.cpp` lines 1468-1488

**Algorithm**:
1. Look up left group ID: `segment_id = map_segments.at(left_key)` (must exist)
2. Look up right group ID: `segment_id2 = map_segments.at(right_key)` (must exist)
3. Add left segment to left group
4. Add right segment to right group

```cpp
segment_id = map_segments.at(pk);          // must exists
// ... add left segment to segment_id group ...

segment_id2 = map_segments.at(pk2);         // must exists
// ... add right segment to segment_id2 group ...
```

**CRITICAL**: Uses `.at()` which throws exception if key not found. This confirms that C++ AGC REQUIRES both split groups to exist before splitting.

## Key Insights

### 1. Groups Must Exist Before Splitting
C++ AGC does NOT create new groups when splitting. It only splits segments to merge them into EXISTING groups. This is why:
- Step 2 uses `operator[]` to look up segment IDs
- Step 7 uses `.at()` which throws if groups don't exist
- The algorithm assumes both split groups are already in `map_segments`

### 2. First Shared K-mer Only
Only the FIRST shared k-mer is used as middle splitter. Other potential middle k-mers are ignored (see comment at line 1554).

### 3. Compression Cost Based
Split position is NOT just the midpoint - it's determined by minimizing the combined compression cost against both reference segments.

### 4. K-byte Overlap Required
Segments overlap by exactly K bytes to ensure both segments have complete terminal k-mers for future matching.

### 5. Degenerate Splits Prevented
If the best split position would create segments smaller than `k+1` bytes, the split is NOT performed.

## Common Questions

**Q: Why check if group exists in map_segments vs just checking terminators?**

A: `map_segments_terminators` tells us which k-mers have been SEEN in segments. `map_segments` tells us which k-mer PAIRS have actual segment groups. A k-mer can appear in terminators without having a group (if it was only seen as a middle k-mer, not a terminal pair).

**Q: Why use operator[] in Step 2 but .at() in Step 7?**

A: This appears to be an inconsistency in the C++ code. Both should use `.at()` to ensure groups exist. Using `operator[]` would insert default values if keys don't exist, leading to invalid segment IDs.

**Q: What happens if no middle splitter is found?**

A: The segment is NOT split and is added as a new singleton group to `map_segments`.

**Q: Can a segment be split into more than 2 parts?**

A: No, C++ AGC only performs binary splits (2 parts). Recursive splitting could theoretically happen if a split segment itself has a middle splitter, but this would occur in a future iteration.

**Q: Why are cost vectors reversed for RC segments?**

A: When using reverse-complement orientation, the segment is conceptually reversed. Cost calculations must account for this to ensure proper alignment with the reference segment.

## Testing Strategy

To verify RAGC matches C++ AGC:

1. **Check terminator map population**: Ensure `map_segments_terminators` contains same k-mers as C++
2. **Check middle splitter selection**: Log which middle k-mers are found and selected
3. **Check group existence**: Verify both split groups exist before attempting split
4. **Check compression costs**: Compare cost vectors with C++ AGC (if possible)
5. **Check split positions**: Log best_pos for sample segments
6. **Check segment sizes**: Verify split segments have correct sizes with k-byte overlap

## References

- `agc_compressor.cpp` lines 1289-1528: `add_segment()` function
- `agc_compressor.cpp` lines 1531-1651: `find_cand_segment_with_missing_middle_splitter()` function
- `segment.cpp` lines 101-116: `get_coding_cost()` function
