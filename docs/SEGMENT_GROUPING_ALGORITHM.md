# C++ AGC Segment Grouping Algorithm Analysis

**Goal**: Exactly match C++ AGC's segment grouping logic to achieve identical compression ratios.

**Current Status**: RAGC creates 70% more segment groups (1226 vs 719), causing 27% archive bloat.

---

## Overview: add_segment() Function (Lines 1279-1815)

The core logic for determining which segment group a segment belongs to.

### Input Parameters
- `segment` - The actual sequence data (contig_t, which is vector<uint8_t>)
- `kmer_front` - The front terminal k-mer (CKmer object)
- `kmer_back` - The back terminal k-mer (CKmer object)

### Output
- `pk` - The segment group key pair<uint64_t, uint64_t> (front_kmer, back_kmer)
- `store_rc` - Boolean indicating if segment should be stored as reverse complement

---

## Phase 1: Determine Segment Group Key Based on K-mer Presence

C++ AGC has 4 distinct cases based on which terminal k-mers are present:

### Case 1: Neither k-mer present (!front.is_full() && !back.is_full()) [Lines 1290-1305]
**Status**: ❓ Need to verify RAGC behavior

```cpp
if (!kmer_front.is_full() && !kmer_back.is_full())
{
    // No terminal splitters present

    if (fallback_filter)       // Try fallback minimizers procedure
    {
        tie(pk, store_rc) = find_cand_segment_using_fallback_minimizers(segment, 1);

        if(pk != pk_empty && store_rc)
            reverse_complement_copy(segment, segment_rc);
    }
    else
        pk = pk_empty;
}
```

**Logic**:
- If fallback mode enabled: Use minimizer-based grouping
- Otherwise: Mark as empty (pk = ~0ull, ~0ull)

**RAGC Implementation**:
```rust
// TODO: Check what RAGC does with MISSING_KMER, MISSING_KMER segments
// Current behavior: ?
```

**Testing**:
- [ ] Create test case with segment having no splitters
- [ ] Verify RAGC creates same group key
- [ ] Verify fallback mode behavior (if enabled)

---

### Case 2: Both k-mers present (front.is_full() && back.is_full()) [Lines 1306-1318]
**Status**: ⚠️ **CRITICAL DIFFERENCE IDENTIFIED**

```cpp
else if (kmer_front.is_full() && kmer_back.is_full())
{
    // Both terminal splitters present

    if (kmer_front.data() < kmer_back.data())
        pk = make_pair(kmer_front.data(), kmer_back.data());
    else
    {
        pk = make_pair(kmer_back.data(), kmer_front.data());  // SWAP!
        reverse_complement_copy(segment, segment_rc);
        store_rc = true;
    }
}
```

**Logic**:
- **NORMALIZE**: Ensure front < back by swapping if needed
- If swapped: Store segment as reverse complement
- **Effect**: Reduces unique group keys by ~50% (approximately)

**RAGC Current Behavior**:
```rust
// segment.rs lines 115-130
let seg_front = front_kmer;
let seg_back = kmer_value;
// NO NORMALIZATION - uses raw values
```

**Why This Matters**:
- Segment (A, B) and (B, A) are treated as same group in C++ AGC
- RAGC treats them as different groups → 70% more groups!

**RAGC Implementation Needed**:
```rust
// In streaming_compressor_queue.rs, when creating SegmentGroupKey:
let (key_front, key_back, should_reverse) = if front != MISSING_KMER && back != MISSING_KMER {
    if front < back {
        (front, back, false)
    } else {
        (back, front, true)  // SWAP + mark for RC
    }
} else {
    (front, back, false)
};

// If should_reverse:
//   1. Reverse complement the segment data
//   2. Set is_rev_comp = true in BufferedSegment
```

**Testing**:
- [ ] Create test with segment having both k-mers where front > back
- [ ] Verify RAGC swaps and stores RC
- [ ] Verify extraction still correct
- [ ] Count unique groups - should drop significantly

---

### Case 3: Only front k-mer present (front.is_full()) [Lines 1319-1340]
**Status**: ❓ Need to analyze

```cpp
else if (kmer_front.is_full())
{
    CKmer kmer = kmer_front;
    reverse_complement_copy(segment, segment_rc);

    tie(pk, store_rc) = find_cand_segment_with_one_splitter(
        kmer, segment, segment_rc, zstd_dctx, bar);

    if (pk.first == ~0ull || pk.second == ~0ull)
    {
        auto pk_alt = pk;
        bool store_rc_alt = false;

        tie(pk_alt, store_rc_alt) = find_cand_segment_using_fallback_minimizers(segment, 5);

        if (pk_alt != pk_empty)
        {
            pk = pk_alt;
            store_rc = store_rc_alt;
        }
    }
}
```

**Logic**:
1. Call `find_cand_segment_with_one_splitter()` to find matching group
2. If no match found (pk contains ~0ull): Try fallback minimizers
3. Returns (pk, store_rc) pair

**Key Function**: `find_cand_segment_with_one_splitter()` - **THIS IS THE CRITICAL MATCHING LOGIC**

**RAGC Current Behavior**:
```rust
// segment.rs lines 117-129
let (seg_front, seg_back) = if front_kmer == MISSING_KMER {
    // Apply orientation logic
    if is_dir {
        (MISSING_KMER, kmer_value)
    } else {
        (kmer_value, MISSING_KMER)
    }
} else {
    (front_kmer, kmer_value)
};
```

**Difference**:
- C++ AGC: Searches existing segments to find best match
- RAGC: Just uses the literal k-mer values

**Testing**:
- [ ] Analyze find_cand_segment_with_one_splitter() behavior
- [ ] Implement equivalent matching logic in RAGC
- [ ] Verify groups match C++ AGC

---

### Case 4: Only back k-mer present (back.is_full()) [Lines 1341-1365]
**Status**: ❓ Need to analyze

```cpp
else if (kmer_back.is_full())
{
    CKmer kmer = kmer_back;
    kmer.swap_dir_rc();  // IMPORTANT: Swap orientation!
    reverse_complement_copy(segment, segment_rc);
    bool store_dir;

    tie(pk, store_dir) = find_cand_segment_with_one_splitter(
        kmer, segment_rc, segment, zstd_dctx, bar);
    store_rc = !store_dir;

    if (pk.first == ~0ull || pk.second == ~0ull)
    {
        auto pk_alt = pk;
        bool store_dir_alt = false;

        tie(pk_alt, store_dir_alt) = find_cand_segment_using_fallback_minimizers(segment_rc, 5);

        if (pk_alt != pk_empty)
        {
            pk = pk_alt;
            store_rc = !store_dir_alt;
        }
    }
}
```

**Logic**:
1. **Swap k-mer orientation** (swap_dir_rc)
2. Call `find_cand_segment_with_one_splitter()` with **segment_rc, segment** (reversed order!)
3. Invert the store_rc result (!store_dir)
4. Fallback to minimizers if needed

**RAGC Current Behavior**:
```rust
// Similar orientation logic but no matching attempt
```

**Testing**:
- [ ] Understand swap_dir_rc() behavior
- [ ] Implement equivalent in RAGC
- [ ] Verify groups match

---

## Phase 2: find_cand_segment_with_one_splitter() [Lines 1637-1815]

**Purpose**: When segment has only one terminal k-mer, find the best existing segment group to match against.

**This is the core matching algorithm that RAGC is missing!**

### Step 1: Check if k-mer is a known terminator [Lines 1646-1658]

```cpp
auto p = map_segments_terminators.find(kmer.data());
if (p == map_segments_terminators.end())
{
    // This k-mer doesn't terminate any existing segment
    // Return default based on k-mer orientation
    if (kmer.is_dir_oriented())
        best_pk = make_pair(kmer.data(), ~0ull);
    else
    {
        best_pk = make_pair(~0ull, kmer.data());
        is_best_rc = true;
    }

    return make_pair(best_pk, is_best_rc);
}
```

**Logic**:
- `map_segments_terminators`: Maps k-mer → set of other k-mers it pairs with
- If k-mer not in map: No existing segments to match → return (kmer, MISSING) or (MISSING, kmer)

**RAGC Needs**:
```rust
// Build map_segments_terminators as segments are added
// HashMap<u64, HashSet<u64>> mapping kmer → set of paired kmers
```

**Testing**:
- [ ] Track which k-mers appear at segment boundaries
- [ ] Build terminator map incrementally
- [ ] Use for matching

---

### Step 2: Normalize candidate pair keys [Lines 1662-1685]

```cpp
for (auto cand_kmer : p->second)
{
    pair<uint64_t, uint64_t> cand_pk;

    v_candidates.emplace_back();
    auto& ck = v_candidates.back();

    if (cand_kmer < kmer.data())
    {
        cand_pk = make_pair(cand_kmer, kmer.data());
        get<0>(ck) = cand_kmer;
        get<1>(ck) = kmer.data();
        get<2>(ck) = true;  // Will need RC
    }
    else
    {
        cand_pk = make_pair(kmer.data(), cand_kmer);
        get<0>(ck) = kmer.data();
        get<1>(ck) = cand_kmer;
        get<2>(ck) = false; // No RC needed
    }

    get<3>(ck) = v_segments[map_segments[cand_pk]];  // Get segment object
}
```

**Logic**:
- For each candidate pairing, **normalize to front < back**
- Store whether RC is needed (get<2>)
- Retrieve the actual segment object for estimation

**Candidate Tuple Structure**:
- `get<0>`: front k-mer (normalized)
- `get<1>`: back k-mer (normalized)
- `get<2>`: bool - needs reverse complement
- `get<3>`: shared_ptr<CSegment> - the reference segment

---

### Step 3: Sort candidates by size similarity [Lines 1688-1697]

```cpp
int64_t segment_size = (int64_t)segment_dir.size();
stable_sort(v_candidates.begin(), v_candidates.end(),
    [segment_size](const auto& x, const auto& y) {
        int64_t x_size = get<3>(x)->get_ref_size();
        int64_t y_size = get<3>(y)->get_ref_size();

        // Primary: Closest size to current segment
        if (abs(segment_size - x_size) < abs(segment_size - y_size))
            return true;
        if (abs(segment_size - x_size) > abs(segment_size - y_size))
            return false;

        // Tiebreaker: Smaller size first
        return x_size < y_size;
    });
```

**Logic**:
- **Prioritize segments with similar size** to current segment
- Hypothesis: Similar size → better compression
- Tiebreaker: Prefer smaller reference

**RAGC Needs**:
```rust
// Track reference segment size for each group
// Sort candidates by |ref_size - current_size|
```

---

### Step 4: Estimate compression for each candidate [Lines 1710-1800]

```cpp
for (auto& candidate : v_candidates)
{
    auto estim_size = get<3>(candidate)->estimate(
        get<2>(candidate) ? segment_rc : segment_dir,
        (uint32_t)best_estim_size,
        zstd_dctx
    );
    auto cand_pk = make_pair(get<0>(candidate), get<1>(candidate));

    if (estim_size < best_estim_size ||
        (estim_size == best_estim_size && cand_pk < best_pk) ||
        (estim_size == best_estim_size && cand_pk == best_pk && !get<2>(candidate)))
    {
        best_estim_size = estim_size;
        best_pk = cand_pk;
        is_best_rc = get<2>(candidate);
    }
}
```

**Logic**:
- Call `segment->estimate()` to predict compressed size using this reference
- Choose candidate with smallest estimated size
- **Tiebreaker 1**: Lexicographically smaller pk
- **Tiebreaker 2**: Prefer non-RC version

**Key Method**: `CSegment::estimate()` - Performs LZ compression simulation

**RAGC Needs**:
```rust
// Implement estimate() that:
//   1. Applies LZ differential encoding against reference
//   2. Estimates ZSTD compressed size
//   3. Returns estimated bytes
```

---

### Step 5: Fallback if no match [Lines 1803-1812]

```cpp
if (best_pk == empty_pk)
{
    if (kmer.is_dir_oriented())
        best_pk = make_pair(kmer.data(), ~0ull);
    else
    {
        best_pk = make_pair(~0ull, kmer.data());
        is_best_rc = true;
    }
}
```

**Logic**:
- If all candidates failed: Use default based on orientation
- This creates a new group with one MISSING k-mer

---

## Implementation Plan

### Milestone 1: Add Case 2 Normalization (Both k-mers present)
**Goal**: Reduce groups by ~40-50% by normalizing (front, back) pairs

**Changes**:
1. In `streaming_compressor_queue.rs`, modify segment grouping logic:
   ```rust
   let (key_front, key_back, should_reverse) =
       if front != MISSING && back != MISSING {
           if front < back {
               (front, back, false)
           } else {
               (back, front, true)
           }
       } else {
           (front, back, false)
       };
   ```

2. Implement reverse complement when should_reverse = true

**Expected Result**: Groups drop from ~1226 to ~600-700
**Success Criteria**:
- [ ] Archive size reduces by ~10-15%
- [ ] Extraction still byte-for-byte correct
- [ ] MD5 hashes match

---

### Milestone 2: Build Terminator Map
**Goal**: Track which k-mers pair together

**Changes**:
1. Add to StreamingQueueCompressor:
   ```rust
   terminator_map: Arc<Mutex<HashMap<u64, HashSet<u64>>>>
   ```

2. When adding segment with both k-mers:
   ```rust
   terminator_map.entry(front).or_default().insert(back);
   terminator_map.entry(back).or_default().insert(front);
   ```

**Expected Result**: Foundation for matching logic
**Success Criteria**:
- [ ] Map populated correctly
- [ ] No performance regression

---

### Milestone 3: Implement find_cand_segment_with_one_splitter()
**Goal**: Match segments with one k-mer to existing groups

**Changes**:
1. Create `find_matching_group()` function
2. Implement size-based candidate sorting
3. Implement estimate() for compression prediction
4. Use for Cases 3 & 4 (one k-mer present)

**Expected Result**: Groups drop to ~719 (matching C++ AGC)
**Success Criteria**:
- [ ] Unique groups match C++ AGC count
- [ ] Archive size within 5% of C++ AGC
- [ ] Extraction correct

---

### Milestone 4: Implement Orientation Logic for Case 4
**Goal**: Handle back-k-mer-only case correctly

**Changes**:
1. Implement swap_dir_rc() equivalent
2. Handle reversed argument order
3. Invert RC flag

**Expected Result**: Complete algorithmic parity
**Success Criteria**:
- [ ] Archive size within 2% of C++ AGC
- [ ] All test cases pass

---

## Testing Strategy

After each milestone:

```bash
# 1. Create archives
./target/release/ragc create -o /tmp/ragc.agc -k 21 -s 10000 -m 20 [input files]
/home/erik/agc/bin/agc create -k 21 -s 10000 -l 20 -o /tmp/cpp.agc [input files]

# 2. Compare sizes
ls -lh /tmp/ragc.agc /tmp/cpp.agc

# 3. Verify correctness
./target/release/ragc getset /tmp/ragc.agc [sample] > /tmp/ragc_extract.fa
/home/erik/agc/bin/agc getset /tmp/cpp.agc [sample] > /tmp/cpp_extract.fa
md5sum /tmp/ragc_extract.fa /tmp/cpp_extract.fa

# 4. Count unique groups
grep "Flushed.*segment groups" in output

# 5. Performance check
time ./target/release/ragc create ...
time /home/erik/agc/bin/agc create ...
```

---

## Current Implementation Gaps

### ❌ Missing in RAGC:
1. **Case 2 normalization** (front < back swap)
2. **Terminator map** for tracking k-mer pairings
3. **find_cand_segment_with_one_splitter()** matching logic
4. **Compression estimation** (estimate() method)
5. **Size-based candidate sorting**
6. **Case 4 orientation handling** (swap_dir_rc)

### ✓ Present in RAGC:
1. Segmentation (creates segments correctly)
2. Basic grouping by (front, back) keys
3. LZ differential encoding
4. ZSTD compression
5. Extraction/decompression

---

## Key Data Structures Needed

### C++ AGC Has:
```cpp
map<uint64_t, set<uint64_t>> map_segments_terminators;
map<pair<uint64_t, uint64_t>, int> map_segments;  // pk → segment_id
vector<shared_ptr<CSegment>> v_segments;          // segment_id → CSegment
```

### RAGC Needs:
```rust
// In StreamingQueueCompressor:
terminator_map: Arc<Mutex<HashMap<u64, HashSet<u64>>>>,
segment_groups_by_key: Arc<Mutex<HashMap<SegmentGroupKey, SegmentGroupInfo>>>,

struct SegmentGroupInfo {
    group_id: u32,
    reference_size: usize,
    reference_data: Vec<u8>,
    // ... other fields for estimation
}
```

---

## Compression Estimation Algorithm

**C++ AGC's estimate() method** (approximate logic):

1. Apply LZ differential encoding:
   - Find longest match between segment and reference
   - Generate copy/literal commands

2. Serialize commands to byte buffer

3. Compress with ZSTD (limited iterations for speed)

4. Return compressed size

**RAGC Implementation**:
- Reuse existing LZ encoder from segment_compression.rs
- Add quick ZSTD estimation
- Return size without actually storing

---

## Progress Tracking

### Milestone 1: Case 2 Normalization
- [ ] Implement swap logic
- [ ] Implement reverse complement
- [ ] Test correctness
- [ ] Measure group count reduction
- [ ] Measure archive size reduction

### Milestone 2: Terminator Map
- [ ] Add data structure
- [ ] Populate during compression
- [ ] Test map contents
- [ ] Verify no performance regression

### Milestone 3: Matching Logic
- [ ] Implement find_matching_group()
- [ ] Implement candidate sorting
- [ ] Implement compression estimation
- [ ] Test group count
- [ ] Measure archive size

### Milestone 4: Complete Parity
- [ ] Handle all 4 cases
- [ ] Pass all correctness tests
- [ ] Achieve <2% size difference
- [ ] Document remaining differences

---

## Notes & Observations

- **Key Insight**: C++ AGC doesn't just group by k-mers, it actively **matches** segments to find best compression
- **Size matters**: Sorting by size similarity suggests segments of similar length compress better together
- **Estimation is critical**: The estimate() method predicts compression, avoiding trial-and-error
- **Normalization everywhere**: front < back normalization appears in multiple places
- **Orientation is complex**: swap_dir_rc() and reversed argument orders show orientation handling is subtle

---

## Questions to Investigate

1. How does CSegment::estimate() work exactly?
2. What is swap_dir_rc() doing with k-mer orientation?
3. Why does fallback use different parameters (1 vs 5)?
4. How is map_segments_terminators built initially?
5. What is the performance impact of compression estimation?

---

**Last Updated**: 2025-11-06
**Status**: Analysis complete, ready for incremental implementation
