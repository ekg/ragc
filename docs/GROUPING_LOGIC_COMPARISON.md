# RAGC vs C++ AGC Grouping Logic Comparison

This document provides detailed guidance on the grouping logic differences between RAGC and C++ AGC. It is intended to help future developers and agents identify and fix grouping discrepancies.

## Summary of Key Differences

| Aspect | C++ AGC | RAGC |
|--------|---------|------|
| Terminator map update | Immediate (per-segment) | Batch-deferred (per-sample) |
| Intra-batch visibility | Groups visible immediately | Groups not visible until batch ends |
| Cost estimation | CLZDiff_V2::Estimate | FFI call to same (when ffi_cost enabled) |
| Candidate sorting | stable_sort by size difference | sort_by with tie-breaker |
| Selection logic | Two-pass with incrementing barriers | Matches C++ AGC |

## Critical Finding: Terminator Map Timing

### The Core Issue

**C++ AGC updates terminators IMMEDIATELY after each group is created.**
**RAGC updates terminators only at BATCH BOUNDARIES (sample changes).**

This means:
- In C++ AGC: Segment N can match with groups created by segments 1..N-1 in same batch
- In RAGC: Segment N **cannot** match with groups created by segments 1..N-1 in same batch

### C++ AGC Behavior (agc_compressor.cpp)

```cpp
// Lines 1017-1024: Update terminators IMMEDIATELY after creating a group
if (kmer1 != ~0ull && kmer2 != ~0ull)
{
    map_segments_terminators[kmer1].push_back(kmer2);
    sort(map_segments_terminators[kmer1].begin(), map_segments_terminators[kmer1].end());

    if (kmer1 != kmer2)
    {
        map_segments_terminators[kmer2].push_back(kmer1);
        sort(map_segments_terminators[kmer2].begin(), map_segments_terminators[kmer2].end());
    }
}
```

### RAGC Behavior (streaming_compressor_queue.rs)

```rust
// Lines 2457-2482: Update BATCH-LOCAL terminators
if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
    let mut term_map = batch_local_terminators.lock().unwrap();
    term_map.entry(key.kmer_front).or_insert_with(Vec::new).push(key.kmer_back);
    // ... (batch-local, not global)
}

// Lines 1700-1710: Merged to global ONLY at batch boundary (sample change)
let mut global_terms = map_segments_terminators.lock().unwrap();
for (kmer, connections) in batch_terms.iter() {
    let entry = global_terms.entry(*kmer).or_insert_with(Vec::new);
    entry.extend(connections.iter().cloned());
    // ...
}
```

### Impact

This timing difference causes different grouping decisions:

1. **First sample**: Both implementations behave identically (no prior groups)
2. **Subsequent samples**: May differ if:
   - A segment would match with a group created earlier in the same batch
   - C++ AGC finds the match, RAGC doesn't (creates new group instead)

## Critical Code Sections

### C++ AGC Key Files and Lines

| File | Lines | Function | Purpose |
|------|-------|----------|---------|
| `agc_compressor.cpp` | 1275-1500 | `add_segment` | Main segment classification entry point |
| `agc_compressor.cpp` | 1625-1802 | `find_cand_segment_with_one_splitter` | Case 3a/3b: one terminator present |
| `agc_compressor.cpp` | 1807-1950 | `find_cand_segment_using_fallback_minimizers` | Fallback grouping strategy |
| `agc_compressor.cpp` | 1017-1024 | (inside add_segment) | **Immediate terminator update** |
| `agc_compressor.cpp` | 1676-1685 | `stable_sort` lambda | Candidate sorting by size difference |
| `agc_compressor.cpp` | 1726-1732 | Pass 1 | Estimate computation with incrementing barriers |
| `agc_compressor.cpp` | 1775-1788 | Pass 2 | Selection with tie-breaking |

### RAGC Key Files and Lines

| File | Lines | Function | Purpose |
|------|-------|----------|---------|
| `streaming_compressor_queue.rs` | 1161-1427 | `find_group_with_one_kmer` | Case 3a/3b equivalent |
| `streaming_compressor_queue.rs` | 1434-1620 | `find_cand_segment_using_fallback_minimizers` | Fallback grouping equivalent |
| `streaming_compressor_queue.rs` | 2454-2482 | (inside `or_insert_with`) | **Batch-local terminator update** |
| `streaming_compressor_queue.rs` | 1700-1710 | `flush_batch` | **Merge batch-local to global** |
| `streaming_compressor_queue.rs` | 1275-1289 | `candidates.sort_by` | Candidate sorting |
| `streaming_compressor_queue.rs` | 1307-1367 | Pass 1 | Estimate computation |
| `streaming_compressor_queue.rs` | 1377-1396 | Pass 2 | Selection with tie-breaking |

## Potential Fixes

### Option 1: Immediate Terminator Update (Recommended)

Update global terminators immediately after creating a new group, matching C++ AGC:

```rust
// In the or_insert_with closure for new groups:
if key.kmer_front != MISSING_KMER && key.kmer_back != MISSING_KMER {
    // Update GLOBAL terminators immediately (like C++ AGC)
    let mut global_terms = map_segments_terminators.lock().unwrap();
    global_terms.entry(key.kmer_front).or_insert_with(Vec::new).push(key.kmer_back);
    if key.kmer_front != key.kmer_back {
        global_terms.entry(key.kmer_back).or_insert_with(Vec::new).push(key.kmer_front);
    }
    // Sort after insertion
    if let Some(v) = global_terms.get_mut(&key.kmer_front) { v.sort_unstable(); v.dedup(); }
    if key.kmer_front != key.kmer_back {
        if let Some(v) = global_terms.get_mut(&key.kmer_back) { v.sort_unstable(); v.dedup(); }
    }
}
```

### Option 2: Look Up Both Maps

Modify `find_group_with_one_kmer` to look in both global AND batch-local terminators:

```rust
let connected_kmers = {
    let global_terms = map_segments_terminators.lock().unwrap();
    let batch_terms = batch_local_terminators.lock().unwrap();

    let mut connections: Vec<u64> = Vec::new();
    if let Some(vec) = global_terms.get(&kmer) {
        connections.extend(vec.iter().cloned());
    }
    if let Some(vec) = batch_terms.get(&kmer) {
        connections.extend(vec.iter().cloned());
    }
    connections.sort_unstable();
    connections.dedup();
    connections
};
```

## Debugging Checklist

When investigating grouping discrepancies:

1. **Enable verbose logging**:
   ```bash
   ./ragc create -v 2 ... 2>&1 | grep -E "(NEW_GROUP|REUSE_GROUP|CASE3)"
   ```

2. **Export segment layouts**:
   ```bash
   ./ragc inspect archive.agc --segment-layout > layout.csv
   ```

3. **Compare CSV layouts**:
   ```bash
   diff cpp_layout.csv ragc_layout.csv | head -50
   ```

4. **Check for intra-batch group creation**:
   Look for patterns where C++ AGC reuses a group but RAGC creates a new one, especially for segments within the same sample.

5. **Trace specific segments**:
   ```bash
   RAGC_TRACE_SAMPLE="sample_name" RAGC_TRACE_CONTIG="chrI" RAGC_TRACE_INDEX=5 ./ragc create ...
   ```

## Key Invariants (Must Always Match)

1. **Cost estimation**: With `ffi_cost` feature enabled, RAGC uses exact same cost calculation
2. **Selection threshold**: `segment_len - 16` (both implementations)
3. **Tie-breaking order**: estimate < best_estimate, then pk < best_pk, then !needs_rc
4. **MISSING_KMER sentinel**: `u64::MAX` / `~0ull`
5. **Group key normalization**: smaller k-mer first, larger second

## Edge Cases to Watch

1. **Empty candidate lists**: When no terminators found, both create new terminator group
2. **All candidates fail threshold**: Both fall back to terminator group with MISSING
3. **Same cost, different pk**: Lexicographic comparison should match
4. **Same cost, same pk, different RC**: Prefer forward orientation
5. **First segment in sample**: Should have identical behavior (no prior batch-local state)

## Test Verification Protocol

After any fix:

```bash
# 1. Build
cargo build --release

# 2. Create archives
./target/release/ragc create -o ragc_test.agc -k 21 -s 10000 -m 20 -t 1 samples/*.fa
/home/erik/agc/bin/agc create -o cpp_test.agc -k 21 -s 10000 -l 20 -t 1 samples/*.fa

# 3. Compare segment layouts
./target/release/ragc inspect ragc_test.agc --segment-layout > ragc_layout.csv
./target/release/ragc inspect cpp_test.agc --segment-layout > cpp_layout.csv
diff ragc_layout.csv cpp_layout.csv

# 4. Compare archive sizes
ls -la ragc_test.agc cpp_test.agc

# 5. Verify round-trip correctness
./target/release/ragc getset ragc_test.agc sample_name > extracted.fa
diff original.fa extracted.fa
```

## Related Documentation

- `CLAUDE.md`: Project-level instructions and verification protocol
- `docs/PIPELINE_COMPARISON.md`: Architecture comparison
- `docs/MEMORY_PROFILING.md`: Performance analysis

## Testing Results

### Fix 1: Immediate Terminator Updates (2025-11-25)

Changed RAGC to update global `map_segments_terminators` immediately when creating new groups (matching C++ AGC's behavior at lines 1017-1024).

**Result: No improvement** - Archive still 11.09% larger.

This indicates the terminator timing is NOT the root cause.

### Fix 2: Use Actual is_dir Orientation (2025-11-25)

**ROOT CAUSE FOUND**: RAGC was hardcoding `kmer_is_dir` values instead of using actual orientation from segment detection.

**The Bug**:
- Case 3a: Always used `kmer_is_dir = true`
- Case 3b: Always used `kmer_is_dir = false`

**The Fix**:
- Case 3a: Use `segment.front_kmer_is_dir` directly
- Case 3b: Use `!segment.back_kmer_is_dir` (inverted due to `swap_dir_rc()`)

**Why This Matters**:
C++ AGC's `swap_dir_rc()` swaps the k-mer's internal `kmer_dir` and `kmer_rc` fields, which inverts `is_dir_oriented()` (which checks `kmer_dir <= kmer_rc`). This determines whether the MISSING key is `(kmer, MISSING)` or `(MISSING, kmer)`.

**Result**:
- AAB#0#chrIX,0 now correctly reuses group 40 instead of creating new group
- Archive size: 3,791,829 bytes (3.4% larger than C++ AGC's 3,667,525)
- Group count: 1352 vs 1349 (only 3 extra groups)

## Key Code Changes

In `streaming_compressor_queue.rs` lines 1568-1606:

```rust
// Case 3a: Use actual orientation from segment detection
find_group_with_one_kmer(
    segment.front_kmer,
    segment.front_kmer_is_dir,  // Was: hardcoded true
    ...
)

// Case 3b: Invert orientation (matches swap_dir_rc behavior)
let kmer_is_dir_after_swap = !segment.back_kmer_is_dir;
find_group_with_one_kmer(
    segment.back_kmer,
    kmer_is_dir_after_swap,  // Was: hardcoded false
    ...
)
```

## Remaining Differences

The 3.4% size difference (3 extra groups) is due to:
1. Some candidate selection differences where RAGC doesn't find a winner
2. Different LZ encoding cost estimates causing different grouping choices

## Next Investigation Areas

1. **Cost estimation accuracy**: Verify FFI cost function produces identical results
2. **Candidate sorting**: Ensure stable_sort behavior matches C++ AGC
3. **Threshold handling**: Check edge cases in threshold comparison

## Conclusion

The major grouping bug has been fixed. RAGC now produces archives within 3.4% of C++ AGC's size (down from 11%) with nearly identical group counts (1352 vs 1349).
