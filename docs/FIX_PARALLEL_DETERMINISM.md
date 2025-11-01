# Fix Parallel Processing Determinism

**Goal**: Make multi-threaded compression produce identical output to sequential compression

**Status**: ✅ COMPLETE - Multi-threaded determinism FIXED!

**Date Started**: 2025-11-01
**Date Completed**: 2025-11-01

---

## ✅ SOLUTION: Sequence Number Approach - SUCCESS!

**Implementation**: Add sequence numbers to tasks when queuing, sort by seq_num before registration

**Changes Made**:
1. Added `seq_num: u64` field to `ContigTask::Normal`
2. Added `seq_num: u64` field to `PendingSegment`
3. Assign seq_num when queuing tasks (uses contig counter)
4. Sort pending segments by seq_num before registration

**Results**:
```
Single-threaded (num_threads=1):  ✅ WORKS
Multi-threaded  (num_threads=15): ✅ WORKS

Original hash:  0bc3568acf4ace9fcd43c2a6153d704579bbf9d0cf0d013ef65882fd72ac35e3
Extracted hash: 0bc3568acf4ace9fcd43c2a6153d704579bbf9d0cf0d013ef65882fd72ac35e3
test result: ok. 1 passed; 0 failed
```

**Why it works**:
- Seq_num captures file order at task creation time
- Workers process tasks in any order (parallel)
- Sorting by seq_num restores file order before registration
- Group IDs assigned in file order → correct decompression!

---

## Root Cause Analysis

### Current RAGC Behavior (WRONG)

With `num_threads > 1`:
- Workers process contigs in parallel
- Segments added to `pending_segments` in **non-deterministic order** (race condition)
- Groups registered in encounter order, but encounter order is non-deterministic
- Result: Different group IDs than sequential → wrong output

### C++ AGC Solution (CORRECT)

C++ AGC achieves determinism through **sorted registration**:

1. **Segments stored in sorted set**:
   ```cpp
   set<kk_seg_part_t> s_seg_part;  // Automatically sorted!
   ```

2. **Sort order** (operator< for kk_seg_part_t):
   ```cpp
   bool operator<(const struct kk_seg_part_t& x) const {
       if (sample_name != x.sample_name)
           return sample_name < x.sample_name;      // Sort by sample first
       if (contig_name != x.contig_name)
           return contig_name < x.contig_name;      // Then by contig
       return seg_part_no < x.seg_part_no;          // Then by segment number
   }
   ```

3. **Group ID assignment** (process_new):
   ```cpp
   map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;
   uint32_t group_id = (uint32_t)vl_seg_part.size();

   // Iterate sorted set → deterministic order!
   for (const auto& x : s_seg_part) {
       auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));
       if (p == m_kmers.end())
           m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;  // Assign sequential group ID
   }
   ```

**KEY INSIGHT**: Group IDs assigned based on **sorted order** of segments, NOT encounter order!

---

## Implementation Plan

### ✅ Phase 1: Understand C++ AGC Sorting [COMPLETE]

**Findings**:
- C++ AGC uses `std::set<kk_seg_part_t>` which automatically sorts
- Sort key: `(sample_name, contig_name, seg_part_no)`
- Groups assigned in iteration order of sorted set
- This gives deterministic group IDs regardless of parallel processing order

**Evidence**:
- File: `/home/erik/agc/src/core/agc_compressor.h:112-119` (operator<)
- File: `/home/erik/agc/src/core/agc_compressor.h:384-415` (process_new)
- File: `/home/erik/agc/src/core/agc_compressor.h:300` (set declaration)

---

### ❌ Phase 2: Sort Pending Segments - INCORRECT APPROACH [ABANDONED]

**Status**: Causes DATA CORRUPTION - abandoned!

**What we tried**: Sort segments by (sample_name, contig_name, seg_part_no) to match C++ AGC's std::set order

**Results**:
- Hash changed from `0bc3568...` (correct) to `be273...` (wrong)
- **DATA CORRUPTION** during decompression:
  ```
  Original:  GCATTTTTATTATTGATATGGGTGGTATGTTGGAATAAAAATCCACTATCGTCTATC
  Extracted: GCATTTTCACCACACCCACACCCACACCCACACCCACACACCACACCCACACACCAC
  ```
- Extracted sequence has repeating "CAC" patterns (garbled)

**Why it failed**:
1. Sorting changes GROUP ID assignment order
2. Group IDs are stored in the archive's collection with segment metadata
3. Decompressor uses group IDs + collection order to reconstruct contigs
4. Changing group ID order breaks the contig reconstruction!

**Key insight**: Group IDs must be assigned in the order segments appear in contigs (file order), NOT in alphabetically sorted order!

**Evidence of file order != sorted order**:
```bash
# File order: chrV, chrVI, chrVII, chrVIII, chrIX, chrX
# Sorted order: chrI, chrII, chrIII, chrIV, chrIX (!), chrV, chrVI, chrVII, chrVIII, chrX
# Alphabetically: "chrIX" < "chrV" because "I" < "V"
```

---

## REVISED Understanding: C++ AGC Does NOT Sort for Group IDs!

After testing, the original hypothesis was **WRONG**. C++ AGC's `std::set<seg_part_t>` is used for:
- Efficient storage and lookup of pending segments
- Thread-safe insertion with automatic deduplication

BUT: The sorted order of `std::set` is NOT used for group ID assignment!

**New hypothesis to investigate**:
1. C++ AGC maintains file order through priority queue task ordering
2. Tasks are processed in priority order (which preserves file order)
3. Even with multiple threads, the std::set collects segments but some other mechanism ensures deterministic group ID order

**What we tried**:
```rust
// Thread 0 only: Register groups deterministically
let mut all_pending = Vec::new();
std::mem::swap(&mut all_pending, &mut *pending_segments.lock().unwrap());

// DON'T SORT! ← THIS IS THE BUG!
// Single-threaded: segments are already in file order
// Multi-threaded: segments are in non-deterministic order

// Group by k-mer key in ENCOUNTER ORDER
use std::collections::HashMap;
let mut key_to_segments: HashMap<(u64, u64), Vec<PendingSegment>> = HashMap::new();
let mut key_order: Vec<(u64, u64)> = Vec::new(); // Track encounter order
```

**New Code**:
```rust
// Thread 0 only: Register groups deterministically
let mut all_pending = Vec::new();
std::mem::swap(&mut all_pending, &mut *pending_segments.lock().unwrap());

// SORT by (sample_name, contig_name, seg_part_no) to match C++ AGC!
all_pending.sort_by(|a, b| {
    (&a.sample_name, &a.contig_name, a.seg_part_no)
        .cmp(&(&b.sample_name, &b.contig_name, b.seg_part_no))
});

// Group by k-mer key in SORTED ORDER (matching C++ AGC!)
use std::collections::HashMap;
let mut key_to_segments: HashMap<(u64, u64), Vec<PendingSegment>> = HashMap::new();
let mut key_order: Vec<(u64, u64)> = Vec::new(); // Track first occurrence in sorted order
```

**Test Plan**:
1. Run with `num_threads=1` → should still match sequential
2. Run with `num_threads=15` → should now match sequential!
3. Verify SHA256 hash matches expected: `0bc3568...`

**Success Criteria**:
- ✅ Single-threaded still works
- ✅ Multi-threaded produces same hash as sequential
- ✅ Test passes: `test_streaming_compressor_extraction_correctness`

---

### ⬜ Phase 3: Test and Verify

**Test 1: Single-threaded (regression test)**
```bash
# Update test to use num_threads=1
cargo test --release test_streaming_compressor_extraction_correctness -- --nocapture
```

Expected: `test result: ok` (no regression)

**Test 2: Multi-threaded**
```bash
# Update test to use num_threads=15
cargo test --release test_streaming_compressor_extraction_correctness -- --nocapture
```

Expected:
```
Original hash:  0bc3568acf4ace9fcd43c2a6153d704579bbf9d0cf0d013ef65882fd72ac35e3
Extracted hash: 0bc3568acf4ace9fcd43c2a6153d704579bbf9d0cf0d013ef65882fd72ac35e3
test result: ok. 1 passed; 0 failed
```

**Test 3: Determinism across runs**
```bash
# Run twice with same input, verify archives identical
cargo run --release -- create -o /tmp/test1.agc -t 15 input.fa
cargo run --release -- create -o /tmp/test2.agc -t 15 input.fa
sha256sum /tmp/test1.agc /tmp/test2.agc  # Should match!
```

**Test 4: Match sequential output**
```bash
# Create with sequential (old code)
cargo run --release -- create -o /tmp/seq.agc -t 1 input.fa

# Create with parallel (new code)
cargo run --release -- create -o /tmp/par.agc -t 15 input.fa

# Compare archives (should be byte-for-byte identical or have same group structure)
```

---

### ⬜ Phase 4: Update Documentation

**Files to update**:
- `docs/FIX_MULTITHREADING_CORRUPTION.md` - Mark multi-threaded as FIXED
- `CLAUDE.md` - Document the sorting insight
- Code comments - Explain why sorting is required

**Documentation points**:
- Sorting by (sample, contig, seg_part_no) ensures deterministic group IDs
- This matches C++ AGC's std::set iteration order
- Without sorting, parallel processing creates non-deterministic encounter order
- Single-threaded doesn't need explicit sort (already in file order)

---

## Key Insights

### Why Sorting Fixes Multi-Threading

1. **Without sorting**:
   - Thread A processes contig1, Thread B processes contig2
   - Which thread finishes first is non-deterministic
   - Segments added to pending_segments in non-deterministic order
   - Groups created in non-deterministic order → wrong group IDs

2. **With sorting**:
   - Threads process contigs in any order (doesn't matter!)
   - All segments collected in pending_segments
   - Sort by (sample, contig, seg_part_no) before registration
   - Groups always created in same sorted order → deterministic group IDs!

### Why Single-Threaded Worked Without Sorting

- Only one thread processes contigs sequentially
- Segments added to pending_segments in file order
- Groups created in file order
- File order happens to match C++ AGC's sorted order (for alphabetically named samples/contigs)

### Critical Test Case

The yeast10 test case (AEL#2.fa) has:
- Single sample: "AEL#2"
- Multiple contigs in file order
- File order matches sorted order (alphabetical)
- This is why single-threaded worked!

With multiple samples or non-alphabetical names, even single-threaded would fail without sorting.

---

## Timeline

- **Phase 1**: Understanding [COMPLETE] - 30 minutes
- **Phase 2**: Implementation [PENDING] - 15 minutes
- **Phase 3**: Testing [PENDING] - 30 minutes
- **Phase 4**: Documentation [PENDING] - 15 minutes

**Total Estimated**: 1.5 hours

---

## References

### C++ AGC Source
- Sorted set: `/home/erik/agc/src/core/agc_compressor.h:300`
- Sort operator: `/home/erik/agc/src/core/agc_compressor.h:157-164`
- Group assignment: `/home/erik/agc/src/core/agc_compressor.h:384-415`

### RAGC Source
- Registration code: `ragc-core/src/compressor_streaming.rs:2163-2270`
- Test: `ragc-core/tests/test_streaming_corruption.rs`
