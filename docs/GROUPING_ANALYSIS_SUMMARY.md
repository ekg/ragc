# Segment Grouping Analysis Summary

**Date**: 2025-12-04
**Analysis Duration**: ~2 hours
**Files Created**:
1. `CPP_AGC_SEGMENT_GROUPING.md` - Complete C++ AGC specification
2. `SEGMENT_GROUPING_COMPARISON.md` - Side-by-side comparison with RAGC

---

## Key Findings

### 1. Group ID Assignment Order: CRITICAL DIFFERENCE

**C++ AGC** (agc_compressor.h, lines 391-397):
- Iterates `s_seg_part` (a `set<kk_seg_part_t>`)
- Ordering: Pure lexicographic by `(sample_name, contig_name, seg_part_no)`
- **No sample priority** - samples within same batch are interleaved alphabetically

**RAGC** (agc_compressor.rs, lines 2541-2582):
- Iterates sorted `pending_batch_segments` (Vec)
- Ordering: `(sample_priority DESC, sample_name ASC, contig_name ASC, place ASC)`
- **Sample priority first** - earlier samples get lower group IDs regardless of name

**Impact**: If samples are processed in non-alphabetical order, group IDs will differ.

**Example**:
```
Files: zzz.fa, aaa.fa (in that order)

C++ AGC sorting:
  AAA#0/chr1/seg0 → group_id = 16
  ZZZ#0/chr1/seg0 → group_id = 17

RAGC sorting (priority 0 > priority -1):
  ZZZ#0/chr1/seg0 → group_id = 16
  AAA#0/chr1/seg0 → group_id = 17
```

---

### 2. Raw Group Count: CONFIGURATION MISMATCH

**C++ AGC** (agc_compressor.cpp, line 1048):
```cpp
if (group_id < (int)no_raw_groups)  // no_raw_groups = 16
    in_group_id = v_segments[group_id]->add_raw(seg_data, zstd_cctx, zstd_dctx);
else
    in_group_id = v_segments[group_id]->add(seg_data, zstd_cctx, zstd_dctx);
```
- **16 raw groups** (groups 0-15) use simple encoding without LZ
- Group 0 is for orphans (both k-mers MISSING)
- Groups 1-15 are for first k-mer pairs encountered

**RAGC** (agc_compressor.rs, line 2659):
```rust
let is_raw_group = group_id == 0;  // Only group 0 is raw
```
- **1 raw group** (only group 0) for orphan segments
- All other groups use LZ encoding

**Impact**: Groups 1-15 will have different encoding, causing size differences.

---

### 3. Segment Classification: Subtle Timing Difference

**C++ AGC**:
- Lookup in `map_segments` during `add_segment()`
- Only checks **global map** (updated during previous batch's `store_segments()`)
- NEW segments go to `s_seg_part` set
- EXISTING segments go to `vl_seg_part[group_id]` vector

**RAGC**:
- Lookup in `map_segments` AND `batch_local_groups`
- Checks **both global and batch-local** maps
- NEW segments go to `pending_batch_segments` Vec
- EXISTING segments processed immediately

**Impact**: If two segments with same k-mer pair arrive in same batch:
- C++ AGC: Both classified as NEW initially, grouped in `process_new()`
- RAGC: First classified as NEW, second may see it in batch-local map (if processed after first's group assigned)

However, in current RAGC implementation, **both are deferred to `pending_batch_segments`** and sorted together, so this should match C++ AGC.

---

### 4. Priority Queue Order: FASTA Order Tie-breaker

**C++ AGC** (queue.h, line 155):
- `multimap<pair<size_t, size_t>, T>` where key is `(priority, cost)`
- Tie-breaker for same (priority, cost): **multimap insertion order** (implementation-defined)

**RAGC** (agc_compressor.rs, lines 194-212):
- `BinaryHeap<ContigTask>` with custom `cmp()`
- Tie-breaker for same (priority, cost): **sequence number** (FASTA order)

**Impact**: For contigs with same size in same sample, order is:
- C++ AGC: Non-deterministic (insertion order)
- RAGC: Deterministic (FASTA order)

This was **fixed in RAGC** to be deterministic, but may not match C++ AGC's actual order.

---

### 5. Map Update Timing: Implementation Detail

**C++ AGC**:
- `map_segments` updated in `store_segments()` (AFTER group ID assignment)
- Multi-threaded writing with "keep lower group ID" logic

**RAGC**:
- `map_segments` updated in `flush_batch()` (DURING group ID assignment)
- Single-threaded sequential assignment, no conflict resolution needed

**Impact**: None if assignment order is same, but may interact with other differences.

---

## Root Cause Analysis

The **1% size difference** between RAGC (566KB) and C++ AGC (572KB) on yeast chr5 test is likely due to:

1. **Primary**: Sample priority vs lexicographic ordering difference
   - Causes different group ID assignments
   - Cascades to different segment splitting decisions
   - Results in slightly different grouping efficiency

2. **Secondary**: Raw group count difference (1 vs 16)
   - Groups 1-15 get LZ encoding in RAGC, raw in C++ AGC
   - Raw encoding less efficient than LZ for similar segments
   - But RAGC's LZ may be better, offsetting some bloat

3. **Tertiary**: FASTA order determinism
   - May cause different contig processing order
   - Different order → different reference segments chosen
   - Different references → different LZ match quality

---

## Recommended Fixes (Priority Order)

### FIX 1: Remove Sample Priority from Sorting (HIGH PRIORITY)

**File**: `ragc-core/src/agc_compressor.rs`, lines 243-266

**Change**:
```rust
impl Ord for PendingSegment {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // REMOVE sample_priority comparison entirely
        // Sort purely lexicographically like C++ AGC

        // By sample_name
        match self.sample_name.cmp(&other.sample_name) {
            std::cmp::Ordering::Equal => {
                // Then by contig_name
                match self.contig_name.cmp(&other.contig_name) {
                    std::cmp::Ordering::Equal => {
                        // Finally by place (seg_part_no)
                        self.place.cmp(&other.place)
                    }
                    other => other,
                }
            }
            other => other,
        }
    }
}
```

**Do the same** for `BufferedSegment::cmp()` (lines 287-310)

**Rationale**: Match C++ AGC's pure lexicographic ordering within batches.

### FIX 2: Use 16 Raw Groups (MEDIUM PRIORITY)

**File**: `ragc-core/src/agc_compressor.rs`

**Change line 2659**:
```rust
let is_raw_group = group_id < NO_RAW_GROUPS;  // Use 16 raw groups like C++ AGC
```

**Add check elsewhere** for groups 0-15 to use raw encoding.

**Rationale**: Match C++ AGC's raw group configuration.

### FIX 3: Verify Priority Queue Tie-breaker (LOW PRIORITY)

**Test**: Log actual contig processing order in both implementations

**If mismatch**: Investigate C++ AGC's multimap insertion order behavior

**Rationale**: Ensure deterministic order matches between implementations.

---

## Testing Protocol

After each fix:

1. **Build both implementations**:
   ```bash
   cd /home/erik/agc && make clean && make -j8
   cd /home/erik/ragc && cargo build --release
   ```

2. **Create archives**:
   ```bash
   /home/erik/agc/bin/agc create -o /tmp/cpp_test.agc -k 21 -s 10000 -l 20 -t 1 samples/*.fa
   ./target/release/ragc create -o /tmp/ragc_test.agc -k 21 -s 10000 -m 20 -t 1 samples/*.fa
   ```

3. **Compare layouts**:
   ```bash
   ./target/release/ragc inspect /tmp/cpp_test.agc --segment-layout > /tmp/cpp_layout.csv
   ./target/release/ragc inspect /tmp/ragc_test.agc --segment-layout > /tmp/ragc_layout.csv
   diff /tmp/cpp_layout.csv /tmp/ragc_layout.csv
   ```

4. **Compare sizes**:
   ```bash
   ls -lh /tmp/cpp_test.agc /tmp/ragc_test.agc
   ```

5. **If still different**: Enable debug logging and compare group ID assignments:
   ```bash
   AGC_DEBUG_CASE3=1 /home/erik/agc/bin/agc create ... 2> cpp_groups.log
   RAGC_GROUP_LOG=1 ./target/release/ragc create ... 2> ragc_groups.log
   ```

---

## Expected Outcome

After applying fixes:

- **Segment layouts should match** (same group_ids for same segments)
- **Archive sizes should match** (within compression variance)
- **Different run results should be identical** (determinism)

If still different:

- **Investigate segment splitting** (lines 1419-1528 in C++ AGC, 4082-4499 in RAGC)
- **Check terminator update timing** (visible before vs after group assignment)
- **Compare LZ encoding parameters** (match length, cost calculation)

---

## Files Reference

### C++ AGC Key Files
- `src/core/agc_compressor.cpp` - Main compression logic
  - Lines 972-989: `register_segments()`
  - Lines 992-1068: `store_segments()`
  - Lines 1298-1568: `add_segment()`
  - Lines 2114-2168: `compress_contig()`
- `src/core/agc_compressor.h` - Data structures
  - Lines 384-415: `CBufferedSegPart::process_new()`
  - Lines 157-164: Segment sorting order
- `src/common/queue.h` - Priority queue
  - Lines 153-346: `CBoundedPQueue` implementation

### RAGC Key Files
- `ragc-core/src/agc_compressor.rs` - Main compression logic
  - Lines 2505-2720: `flush_batch()`
  - Lines 2743-2895: Main compression loop
  - Lines 3645-3920: Segment processing
  - Lines 243-310: Sorting order definitions

---

## Conclusion

The grouping mechanisms are **fundamentally the same** but have **critical ordering differences**:

1. **Sample priority** (RAGC) vs **pure lexicographic** (C++ AGC)
2. **16 raw groups** (C++ AGC) vs **1 raw group** (RAGC)
3. **Deterministic FASTA order** (RAGC) vs **insertion order** (C++ AGC)

Fixing issue #1 should eliminate most divergence.
Fixing issue #2 will eliminate encoding differences.
Issue #3 may require C++ AGC investigation to match exactly.

**Next steps**: Apply FIX 1, test, and iterate.
