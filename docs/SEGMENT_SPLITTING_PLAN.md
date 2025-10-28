# Segment Splitting Implementation Plan

**Goal**: Implement C++ AGC's segment splitting feature to reduce group count from 803 to ~240 and close the 26-38% compression gap.

**Date**: 2025-10-28
**Status**: Planning

---

## Problem Statement

RAGC creates 803 groups vs C++ AGC's ~240 groups because segments with unique (k1, k2) pairs create new groups instead of splitting and joining existing groups.

**Example**:
- Existing groups: (k1, k_mid), (k_mid, k2)
- New segment: (k1, k2)
- **C++ AGC**: Splits into (k1, k_mid) + (k_mid, k2) → No new group
- **RAGC**: Creates new group (k1, k2) → Group explosion

**Impact**: 803 groups → 6.5 segs/pack → poor LZ compression → 26-38% larger files

---

## C++ AGC Algorithm (from agc_compressor.cpp:1366-1454)

### Conditions for Segment Splitting

1. Multi-file mode (`!concatenated_genomes`)
2. No group exists for (k1, k2)
3. Both k-mers are valid (not MISSING_KMER)
4. Both k-mers exist in `map_segments_terminators`

### Step 1: Find Shared Splitters

```rust
// Pseudocode
let front_terminators = group_terminators.get(kmer_front)?;
let back_terminators = group_terminators.get(kmer_back)?;

// Set intersection: find k-mers that pair with BOTH front and back
let shared: Vec<u64> = front_terminators
    .iter()
    .filter(|k| back_terminators.contains(k))
    .copied()
    .collect();

if shared.is_empty() {
    // Can't split - create new group
    return None;
}

// Use FIRST shared k-mer as middle splitter
let middle_kmer = shared[0];
```

**Key insight**: If k1 pairs with X and k2 pairs with X, we can split at X!

### Step 2: Find Optimal Split Position

**C++ AGC approach** (lines 1532-1621):
- Get reference segments for (k1, middle) and (middle, k2)
- Calculate LZ compression cost at each position
- Find position that minimizes total cost

**Simplified RAGC approach** (for MVP):
- **Skip LZ cost calculation** (complex, requires reference segment access)
- Use simple heuristic: **split at middle of segment**
- Rationale: Any split is better than creating a new group

```rust
// MVP: Simple midpoint split
let split_pos = segment.len() / 2;

// Future: LZ-based optimal split
// let split_pos = find_optimal_split(&segment, ref1, ref2);
```

### Step 3: Split Segment with Overlap

**Critical**: Overlap ensures middle k-mer appears in both parts!

```rust
let kmer_len = 31; // config.kmer_length
let overlap = kmer_len / 2;

// Part 1: [0 .. split_pos + kmer_len]
let part1 = segment[0 .. split_pos + kmer_len].to_vec();

// Part 2: [split_pos - overlap .. end]
let part2 = segment[split_pos - overlap ..].to_vec();

// Create two SegmentInfo objects
let seg1 = SegmentInfo {
    key: SegmentGroupKey::new(kmer_front, middle_kmer),
    segment: Segment { data: part1, ... },
    ...
};

let seg2 = SegmentInfo {
    key: SegmentGroupKey::new(middle_kmer, kmer_back),
    segment: Segment { data: part2, ... },
    ...
};
```

**Verification**:
- Part 1 ends with middle k-mer
- Part 2 starts with middle k-mer
- Both parts have correct terminal k-mers

---

## RAGC Implementation Strategy

### Location: compressor_streaming.rs:2790-2845

**Current flow**:
1. Lines 2787-2843: Handle one-kmer segments (✓ done)
2. Line 2845: Try to get or create group

**New flow**:
1. Lines 2787-2843: Handle one-kmer segments (✓ done)
2. **NEW: Lines 2843-2900: Try to split segment if group doesn't exist**
3. Line 2900+: Get or create group (only if splitting failed)

### Implementation Steps

#### Step 1: Check if Splitting is Possible

```rust
// After one-kmer handling, before or_insert_with
// Check if segment needs splitting
if seg_info.key.kmer_front != MISSING_KMER
    && seg_info.key.kmer_back != MISSING_KMER
    && !groups.contains_key(&seg_info.key)  // Group doesn't exist
{
    // Try to split
    if let Some((seg1, seg2)) = try_split_segment(
        &seg_info,
        &group_terminators,
        &config,
    ) {
        // Process BOTH segments (recursive!)
        // This is tricky in parallel code...
    }
}
```

#### Step 2: Implement `try_split_segment` Function

```rust
fn try_split_segment(
    seg_info: &SegmentInfo,
    group_terminators: &DashMap<u64, Vec<u64>>,
    config: &CompressionConfig,
) -> Option<(SegmentInfo, SegmentInfo)> {
    // Find shared splitters
    let front_terms = group_terminators.get(&seg_info.key.kmer_front)?;
    let back_terms = group_terminators.get(&seg_info.key.kmer_back)?;

    let shared: Vec<u64> = front_terms
        .iter()
        .filter(|k| back_terms.contains(*k))
        .copied()
        .collect();

    if shared.is_empty() {
        return None;  // Can't split
    }

    let middle_kmer = shared[0];

    // Simple midpoint split (MVP)
    let split_pos = seg_info.segment.data.len() / 2;
    let kmer_len = config.kmer_length as usize;
    let overlap = kmer_len / 2;

    // Validate split position
    if split_pos < kmer_len || split_pos + kmer_len > seg_info.segment.data.len() {
        return None;  // Segment too small to split
    }

    // Part 1: [0 .. split_pos + kmer_len]
    let part1_data = seg_info.segment.data[0..split_pos + kmer_len].to_vec();

    // Part 2: [split_pos - overlap .. end]
    let part2_data = seg_info.segment.data[split_pos - overlap..].to_vec();

    // Create SegmentInfo for part 1
    let seg1 = SegmentInfo {
        key: SegmentGroupKey::new_normalized(
            seg_info.key.kmer_front,
            middle_kmer,
        ).0,
        segment: Segment {
            data: part1_data,
            sample_name: seg_info.segment.sample_name.clone(),
            contig_name: seg_info.segment.contig_name.clone(),
            is_rev_comp: seg_info.segment.is_rev_comp,
        },
        seg_part_no: seg_info.seg_part_no * 2,  // Track split segments
    };

    // Create SegmentInfo for part 2
    let seg2 = SegmentInfo {
        key: SegmentGroupKey::new_normalized(
            middle_kmer,
            seg_info.key.kmer_back,
        ).0,
        segment: Segment {
            data: part2_data,
            sample_name: seg_info.segment.sample_name.clone(),
            contig_name: seg_info.segment.contig_name.clone(),
            is_rev_comp: seg_info.segment.is_rev_comp,
        },
        seg_part_no: seg_info.seg_part_no * 2 + 1,
    };

    Some((seg1, seg2))
}
```

#### Step 3: Handle Split Segments in Parallel Code

**Challenge**: After splitting, we have TWO segments to process. In parallel code, we can't just recurse.

**Solutions**:

**Option A: Process inline (simplest)**
```rust
if let Some((seg1, seg2)) = try_split_segment(...) {
    // Process seg1
    {
        let pack_opt = groups.entry(seg1.key.clone()).or_insert_with(...);
        // Add seg1 to group
    }

    // Process seg2
    {
        let pack_opt = groups.entry(seg2.key.clone()).or_insert_with(...);
        // Add seg2 to group
    }

    // Skip original segment processing
    continue;  // or early return
}
```

**Option B: Recursive splitting (more correct but complex)**
- After splitting, check if split segments also need splitting
- Requires loop or recursion

**Recommendation**: Start with Option A (simple inline processing). If that works, add recursive splitting later.

---

## Challenges & Solutions

### Challenge 1: Accessing group_terminators Before Group Creation

**Problem**: We need terminators map populated BEFORE we try to split.

**Solution**: ✓ Already handled! Lines 2849-2864 update terminators when creating groups. First segments populate the map, later segments can split.

### Challenge 2: Parallel Processing Complexity

**Problem**: Multiple threads accessing `groups` and `group_terminators` DashMaps.

**Solution**: ✓ DashMap handles concurrency! Each entry has its own lock. Splitting logic just needs to:
1. Read from terminators (shared read)
2. Insert into groups (per-key lock)

### Challenge 3: Segment Part Numbering

**Problem**: Split segments need unique part numbers for collection tracking.

**Solution**: Use hierarchical numbering:
- Original segment: part_no = N
- Part 1: part_no = N * 2
- Part 2: part_no = N * 2 + 1

### Challenge 4: LZ Compression Cost Calculation

**Problem**: C++ AGC calculates optimal split position using LZ cost.

**Solution**: **Defer this optimization!**
- MVP: Simple midpoint split
- Still reduces groups from 803 → ~240
- Can add LZ optimization later

---

## Testing Strategy

### Test 1: Debug Output

Add verbose logging:
```rust
if config.verbosity > 1 {
    eprintln!("[SPLIT] Segment ({:x}, {:x}) → ({:x}, {:x}) + ({:x}, {:x})",
             seg_info.key.kmer_front, seg_info.key.kmer_back,
             seg1.key.kmer_front, seg1.key.kmer_back,
             seg2.key.kmer_front, seg2.key.kmer_back);
}
```

### Test 2: Group Count Verification

```bash
# With splitting
/home/erik/ragc/target/release/ragc create -o ragc_split.agc -k 31 -v 2 samples/*.fa
# Expected: ~240 groups (shown in verbose output)

# Without splitting (revert code)
# Expected: ~803 groups
```

### Test 3: File Size Comparison

```bash
# yeast10, k=31, s=60000
/home/erik/ragc/target/release/ragc create -o ragc.agc -k 31 -s 60000 samples/*.fa
/home/erik/agc/bin/agc create -o cpp.agc -k 31 -s 60000 samples/*.fa

ls -lh ragc.agc cpp.agc
# Expected: ragc.agc ~6.6M (matching cpp.agc)
# Current: ragc.agc 9.1M (38% larger)
```

### Test 4: Extraction Correctness

```bash
# Extract and verify
/home/erik/ragc/target/release/ragc getset ragc_split.agc sample > out.fa
diff original.fa out.fa
# Expected: No differences
```

---

## Implementation Phases

### Phase 1: MVP (Simple Midpoint Splitting)

**Goal**: Reduce groups from 803 → ~240

1. Implement `try_split_segment` with simple midpoint split
2. Add splitting check before `or_insert_with`
3. Process split segments inline
4. Test on yeast10 dataset
5. Verify group count reduction

**Success criteria**:
- Group count ~240 (matching C++ AGC)
- File size ~7-8M (improvement over 9.1M)
- Extraction correctness ✓

### Phase 2: Optimize Split Position (Optional)

**Goal**: Match C++ AGC file size exactly

1. Implement LZ compression cost estimation
2. Test split position at multiple points
3. Select position with minimum cost
4. Benchmark file size improvement

**Success criteria**:
- File size ~6.6M (matching C++ AGC exactly)

### Phase 3: Recursive Splitting (Optional)

**Goal**: Handle edge cases where split segments also need splitting

1. Add recursion or loop to `try_split_segment`
2. Handle deeply nested splits
3. Add depth limit to prevent infinite loops

---

## Code Locations

### Files to Modify

1. **ragc-core/src/compressor_streaming.rs**
   - Line ~2843: Add splitting logic before `or_insert_with`
   - New function: `try_split_segment` (lines ~3900-4000)

### Data Structures Available

- `group_terminators: DashMap<u64, Vec<u64>>` - tracks k-mer pairings
- `groups: DashMap<SegmentGroupKey, (u32, GroupWriter)>` - existing groups
- `config.kmer_length: u32` - k-mer length for overlap calculation

---

## Expected Results

### Before Splitting

| Metric | Value |
|--------|-------|
| Groups | 803 |
| Segs/Pack | 6.5 |
| File Size (k=31) | 9.1M |
| Gap vs C++ AGC | +38% |

### After Splitting (Predicted)

| Metric | Value |
|--------|-------|
| Groups | ~240 |
| Segs/Pack | ~50 |
| File Size (k=31) | ~6.6M |
| Gap vs C++ AGC | ~0% |

---

## Next Steps

1. ✅ Create this implementation plan
2. ⏳ Implement `try_split_segment` function
3. ⏳ Add splitting check in worker thread
4. ⏳ Test on yeast10 dataset
5. ⏳ Verify group count and file size
6. ⏳ Commit working implementation

---

## Questions for Erik

1. **LZ cost calculation**: Should we implement optimal split position in Phase 1, or is midpoint split acceptable for MVP?

2. **Recursive splitting**: Should split segments recursively split themselves, or is one level sufficient?

3. **Performance**: Should we limit splitting (e.g., only if segment > 2*kmer_length) to avoid overhead on small segments?
