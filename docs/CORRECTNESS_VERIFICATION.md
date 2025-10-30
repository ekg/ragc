# RAGC vs C++ AGC Correctness Verification

**Date**: 2025-10-30
**Status**: CRITICAL BUG FOUND - Systematic verification in progress
**Bug**: Samples after first are corrupted (30-66% data missing)

---

## Verification Protocol

For each component, we verify:
1. **Logic Match**: RAGC implements the same algorithm as C++ AGC
2. **Output Match**: RAGC produces identical intermediate results as C++ AGC
3. **Integration Match**: Component works correctly with other components

---

## Current Bug Symptoms

| Sample | C++ AGC Size | RAGC Size | Data Loss |
|--------|--------------|-----------|-----------|
| AAA#0 | 12,309,305 | 12,309,305 | **0%** ✓ |
| AAB#0 | 12,299,861 | 8,463,610 | **31.1%** ✗ |
| AAC#0 | 12,208,304 | 8,633,675 | **29.2%** ✗ |
| AAC#1 | 11,993,248 | 8,115,025 | **32.3%** ✗ |
| AAC#2 | 11,767,847 | 7,913,733 | **32.7%** ✗ |
| AAR#0 | 12,155,284 | 7,805,245 | **35.7%** ✗ |
| ABA#0 | 12,350,989 | 7,799,508 | **36.8%** ✗ |
| ABH#0 | 12,141,176 | 7,973,761 | **34.3%** ✗ |
| ACA#0 | 12,532,526 | 4,187,512 | **66.5%** ✗ |
| ACH#0 | 12,397,704 | 7,295,631 | **41.1%** ✗ |

**Pattern**: First sample perfect, all others increasingly corrupted.

**Split Statistics**:
- C++ AGC: 7143 attempts → 6431 executed (89.9% success)
- RAGC: 6710 eligible → 5392 executed (80.4% success)
- Missing: ~1039 splits

---

## Component Verification Checklist

### ✓ Phase 1: K-mer Extraction and Splitter Finding

**C++ AGC Logic** (from `CPP_AGC_EXACT_ARCHITECTURE.md`):
- Enumerate k-mers from first sample (AAA#0)
- Find duplicate k-mers
- Find singleton k-mers
- Find actually-used splitters

**RAGC Implementation**: `ragc-core/src/compressor_streaming.rs:add_fasta_files_with_splitters()`

**Verification Status**: ✓ **VERIFIED**
- Test: AAA#0 decompresses identically
- Splitter count: 11,771 splitters (matches C++ AGC logs)
- Singleton k-mers: 11,302,379 (matches C++ AGC logs)
- Duplicate k-mers: 202,849 (matches C++ AGC logs)

**Conclusion**: Splitter finding is correct. Bug is NOT here.

---

### ⚠️ Phase 2: Contig Segmentation

**C++ AGC Logic**:
```cpp
for (each k-mer in contig) {
    if (bloom_splitters.check(d) && hs_splitters.check(d)) {
        add_segment(sample_name, id, seg_part_no, segment_data,
                   split_kmer, kmer, ...);
    }
}
add_segment(...); // Final segment
```

**RAGC Implementation**: `ragc-core/src/contig_compression.rs:compress_contig_with_context()`

**Verification Needed**:
- [ ] Do we create the same segments for AAB#0 as C++ AGC?
- [ ] Do we split at the same positions?
- [ ] Do we call add_segment the same number of times?

**Test**: Add instrumentation to count segments per sample.

---

### ⚠️ Phase 3: Segment Addition and Split Checking

**C++ AGC Logic** (from `agc_compressor.cpp:1275-1504`):
```cpp
pair_segment_desc_t CAGCCompressor::add_segment(...) {
    pk = make_pair(kmer_front.canonical(), kmer_back.canonical());

    auto p = map_segments.find(pk);

    // SPLIT CHECKING
    if (!concatenated_genomes &&
        p == map_segments.end() &&
        pk.first != ~0ull && pk.second != ~0ull &&
        map_segments_terminators.count(pk.first) &&
        map_segments_terminators.count(pk.second)) {

        auto split_match = find_cand_segment_with_missing_middle_splitter(...);

        if (split_match.first != ~0ull) {
            // SPLIT INTO TWO SEGMENTS
            segment_id = map_segments.at(make_pair(kmer_front, middle_kmer));
            segment_id2 = map_segments.at(make_pair(middle_kmer, kmer_back));

            // Buffer both segments as KNOWN
            buffered_seg_part.add_known(segment_id, ...);
            buffered_seg_part.add_known(segment_id2, ...);
            return;
        }
    }

    // No split - buffer as NEW or KNOWN
    if (p == map_segments.end()) {
        buffered_seg_part.add_new(pk.first, pk.second, ...);
    } else {
        buffered_seg_part.add_known(segment_id, ...);
    }
}
```

**RAGC Implementation**: `ragc-core/src/compressor_streaming.rs` (Phase 4: Process KNOWN segments)

**Critical Differences to Check**:
- [ ] Do we check splits for EVERY segment or only NEW segments?
- [ ] Do we require BOTH k-mers in terminators?
- [ ] Do we correctly add BOTH split halves to archive?
- [ ] Do we handle the case where split succeeds?

**Known Issue**:
```
RAGC stats:
  Both target groups exist: 6710
  Splits executed: 5392
  Missing: 1318 splits (19.6% failure rate)
```

**Hypothesis**: **Splits are being detected but not executed correctly**, causing segments to be dropped.

---

### ⚠️ Phase 4: Split Execution Logic

**What happens when a split executes?**

**C++ AGC**:
1. Find middle k-mer
2. Split segment into two parts at middle k-mer position
3. Look up group_id1 for (kmer_front, middle_kmer)
4. Look up group_id2 for (middle_kmer, kmer_back)
5. Buffer BOTH parts as KNOWN segments with their group_ids
6. **Both parts are added to archive**

**RAGC** (`compressor_streaming.rs` lines ~6000-6500):
Need to check:
- [ ] After split, do we add BOTH segments?
- [ ] Do we correctly look up both group IDs?
- [ ] Do we return both segments or only one?
- [ ] Do we handle split failure correctly?

**Critical Code to Review**: `split_and_add_segment_known()` or equivalent

---

### ⚠️ Phase 5: Buffered Segment Handling

**C++ AGC**:
- Segments buffered during processing phase
- Flushed during `store_segments()` at synchronization barrier
- ALL buffered segments are added to archive

**RAGC**:
- Need to verify buffering and flushing logic

**Verification Needed**:
- [ ] Are all buffered segments flushed?
- [ ] Are split segments properly buffered?
- [ ] Is there any path where segments get dropped?

---

### ⚠️ Phase 6: Terminator Addition

**C++ AGC Logic** (`store_segments()` lines 995-1049):
```cpp
if (kmer1 != ~0ull && kmer2 != ~0ull) {
    map_segments_terminators[kmer1].push_back(kmer2);
    sort(map_segments_terminators[kmer1].begin(),
         map_segments_terminators[kmer1].end());

    if (kmer1 != kmer2) {
        map_segments_terminators[kmer2].push_back(kmer1);
        sort(map_segments_terminators[kmer2].begin(),
             map_segments_terminators[kmer2].end());
    }
}
```

**RAGC Implementation**: Need to find equivalent

**Verification Needed**:
- [ ] Do we add terminators for ALL new groups?
- [ ] Do we add bidirectional connections (k1→k2 AND k2→k1)?
- [ ] Do we add terminators at the right time (before next sample)?
- [ ] Are terminators from AAA#0 available when processing AAB#0?

---

## Debugging Strategy

### Step 1: Instrument Split Execution ✓ DONE

Already have logs showing:
- Split attempts
- Split successes
- Split failures

### Step 2: Find Where Segments Are Dropped

Add counters for:
```rust
// In segment processing
total_segments_input: usize,
segments_added_as_new: usize,
segments_added_as_known: usize,
segments_split: usize,
segments_after_split: usize, // Should be segments_split * 2

// Check invariant:
assert_eq!(total_segments_input,
           segments_added_as_new +
           segments_added_as_known +
           segments_after_split);
```

If this fails, we're dropping segments somewhere.

### Step 3: Compare Group Structure

Extract from both archives:
- Total number of groups
- Segments per group
- K-mer pairs for each group

Compare:
```bash
cpp_agc_groups=$(agc info cpp.agc | parse groups)
ragc_groups=$(agc info ragc.agc | parse groups)

diff cpp_agc_groups ragc_groups
```

### Step 4: Compare Segment Sequences

For AAB#0, extract:
- C++ AGC: List of (start_pos, end_pos, group_id) for each segment
- RAGC: Same list
- Compare: Are segments assigned to different groups?

### Step 5: Binary Search for Divergence Point

Process samples one at a time:
1. Create archive with just AAA#0 → verify perfect
2. Create archive with AAA#0 + AAB#0 → check if AAB#0 corrupts
3. If yes, binary search within AAB#0 to find which contig/segment causes divergence

---

## Next Actions

1. **Find split execution code** in RAGC
2. **Add segment accounting** to track every segment from input to output
3. **Compare first divergent segment** between RAGC and C++ AGC
4. **Identify exact line** where behavior differs

---

## Expected Fix

Once we find where segments are dropped:
- Fix the split execution logic
- Ensure all split halves are added
- Verify all samples decompress identically
- Create regression test with all 10 samples
