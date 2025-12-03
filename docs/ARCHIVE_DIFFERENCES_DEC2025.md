# Archive Differences Investigation - December 2025

**Date**: 2025-12-03
**Status**: In Progress
**Goal**: Achieve byte-for-byte identical archives with C++ AGC

---

## Current Situation

**Archives**:
- RAGC: `/tmp/ragc_clean.agc` - 15M
- C++ AGC: `/tmp/cpp_clean.agc` - 12M
- **Gap: 25% larger**

**Test Dataset**: 10 yeast samples (chrV only)

**`ragc inspect --compare` Results**:
- RAGC: 12,980 segments
- C++ AGC: 12,987 segments
- Four major differences identified (see below)

---

## Four Root Causes (Priority Order)

### Issue #2: Group_ID Assignment (PRIORITY 1 - MAJOR)

**Impact**: 5,650 segments in different groups (43.5% of segments!)

**Description**:
- C++ AGC collects NEW segments in a sorted BTreeSet (by sample→contig→place order)
- C++ AGC assigns group_ids sequentially when iterating this sorted set
- RAGC was assigning group_ids immediately via atomic counter (random order)
- **Result**: Identical k-mer pairs get different group_ids → poor LZ compression

**Root Cause**:
- RAGC: `group_counter.fetch_add(1)` at registration time (immediate, unsorted)
- C++ AGC: `process_new()` iterates sorted `s_seg_part` and assigns IDs sequentially

**Previous Fix Attempt** (Commit ff00e87):
- Implemented sorted group_id assignment in `flush_batch()`
- **BROKE ARCHIVES**: Created archives with "Unknown frame descriptor" error
- Pushed to origin/main but reverted locally with `git reset --hard HEAD~1`
- Result: 14M archive (vs 15M without fix) - improvement but broken

**Investigation Needed**:
1. Why did commit ff00e87 break archives?
2. What went wrong in the implementation?
3. How to implement sorted assignment WITHOUT breaking compatibility?

**Investigation Results** (2025-12-03):

Analyzed commit ff00e87. The approach was:
- EXISTING groups: Segments added to `segment_groups` buffer immediately in worker_thread
- NEW groups: Segments deferred to `pending_batch_segments` vector
- In `flush_batch()`: Pending segments sorted by (sample, contig, place), group_ids assigned, segments added to buffers

**SUSPECTED BUG**: In the new `flush_batch()` implementation, after adding pending segments to `segment_groups` buffers, **NO WRITE/FLUSH OPERATIONS OCCUR**.

The old code (worker_thread) would:
1. Add segment to buffer
2. Call `flush_pack()` if buffer was full (writes to archive)

The new code (flush_batch) does:
1. Sort pending segments
2. Assign group_ids
3. Add segments to buffers
4. **RETURN** ← Missing flush_pack calls!

This means pending segments are buffered but **never written to the archive**, causing incomplete/corrupted archives that fail with "Unknown frame descriptor" when reading.

**Root Cause**: Missing `flush_pack()` calls in `flush_batch()` for newly created segments.

**Status**: ❌ Not Fixed (broken implementation reverted)

---

### Issue #1: Segment Splitting Logic (PRIORITY 2)

**Impact**: 905 segment index mismatches

**Description**:
- After initial segmentation by splitters, segments can be split further based on compression cost heuristics (`try_split_segment_with_cost()`)
- C++ AGC and RAGC use different heuristics or calculations
- **Example**: C++ AGC keeps 20KB segment whole, RAGC splits it into two 10KB segments
- This causes segment index misalignment in comparisons

**Root Cause**: Unknown - need to compare cost calculation between implementations

**Note**: `ragc inspect --compare` uses segment INDEX as key, which becomes misleading when one archive splits differently

**Investigation Needed**:
1. Compare compression cost calculation formula
2. Compare split decision threshold
3. Match RAGC's split logic to C++ AGC exactly

**Status**: ❌ Not Investigated

---

### Issue #3: Reference Selection (PRIORITY 3)

**Impact**: 5 groups have different reference selection

**Description**:
- Within a segment group, one segment is chosen as the LZ encoding reference
- All other segments are encoded as diffs from this reference
- C++ AGC and RAGC select different segments as reference in 5 groups
- **Result**: Different reference → more literals, fewer matches → larger archive

**Root Cause**: Unknown - need to identify reference selection algorithm in C++ AGC

**Investigation Needed**:
1. Find C++ AGC's reference selection code
2. Document the algorithm (first segment? longest? most common k-mers?)
3. Match RAGC's selection logic

**Status**: ❌ Not Investigated

---

### Issue #4: Pack Structure (PRIORITY 4)

**Impact**: 286 groups have different pack counts

**Description**:
- Segments are organized into ZSTD packs of ~50 segments each
- C++ AGC and RAGC organize packs differently in 286 out of ~1300 groups
- **Result**: Different pack boundaries → different ZSTD compression contexts → different sizes

**Root Cause**: Unknown - need to compare pack organization logic

**Investigation Needed**:
1. Understand how C++ AGC decides pack boundaries
2. Compare with RAGC's packing logic
3. Match the packing decisions exactly

**Status**: ❌ Not Investigated

---

## Investigation Protocol

**From CLAUDE.md**:
1. **Use `ragc inspect` FIRST** - extract as much info as possible from archives themselves
2. **Only instrument code** when absolutely no other way to see the difference
3. **Fix bugs in `ragc inspect`** along the way to improve investigation tools
4. **Be systematic** - fix one issue at a time, verify impact on archive size

**Order of Investigation**:
1. **#2: Group_ID Assignment** (5650 segments affected - MAJOR)
2. **#1: Segment Splitting Logic** (905 mismatches)
3. **#3: Reference Selection** (5 groups)
4. **#4: Pack Structure** (286 groups)

---

## Investigation History

### 2025-12-03 (Morning)

**Commit ff00e87** - Attempted sorted group_id assignment:
- Modified `agc_compressor.rs` to collect NEW segments in `pending_batch_segments`
- Modified `flush_batch()` to sort by (sample_name, contig_name, place) and assign group_ids
- **Result**: 14M archive (improvement from 15M)
- **Problem**: Archives broken - "Unknown frame descriptor" error on read
- **Action**: Reverted with `git reset --hard HEAD~1`

**Commit 71c53d1** - Fixed `ragc inspect --compare` output:
- Changed "DIFFERENT segment boundaries" → "DIFFERENT segment splitting"
- Added clarifying note about index comparison vs genomic position
- Made comparison output less misleading

### 2025-12-03 (Afternoon)

Current state:
- Many background bash processes running (C++ AGC tests with instrumentation)
- Testing various configurations and approaches
- Documents created: `SIZE_DIFFERENCE_INVESTIGATION.md` (outdated, from earlier phase)

---

## Next Actions

**Immediate** (#2 Group_ID Assignment):
1. Understand WHY commit ff00e87 broke archives
2. Review the implementation in `agc_compressor.rs` and `worker.rs`
3. Identify the specific bug that caused "Unknown frame descriptor"
4. Design a fix that implements sorted assignment WITHOUT breaking format

**After #2 is fixed**:
1. Verify archive size improvement (expect 14M, same as ff00e87)
2. Verify archives can be read by both RAGC and C++ AGC
3. Run full roundtrip tests
4. Move to #1 (Segment Splitting Logic)

---

## Success Criteria

- [ ] **#2 Fixed**: Group_ids assigned in sorted order, archives readable
- [ ] **#1 Fixed**: Segment splitting matches C++ AGC
- [ ] **#3 Fixed**: Reference selection matches C++ AGC
- [ ] **#4 Fixed**: Pack structure matches C++ AGC
- [ ] **Final**: Archives are byte-for-byte identical (or < 1% difference with documented reasons)

---

## Key Files

**Investigation Tools**:
- `ragc-cli/src/inspect.rs` - Archive inspection and comparison
- `docs/SIZE_DIFFERENCE_INVESTIGATION.md` - Earlier phase (outdated)
- `docs/ARCHIVE_DIFFERENCES_DEC2025.md` - This document

**Implementation**:
- `ragc-core/src/agc_compressor.rs` - Group_id assignment logic
- `ragc-core/src/worker.rs` - Segment packing and processing
- `ragc-core/src/segment.rs` - Segment splitting logic

**Commits**:
- `ff00e87` (origin/main, BROKEN) - Sorted group_id assignment attempt
- `71c53d1` (local main) - Fixed inspect --compare output
- `8f8dd19` (parent) - Working baseline before ff00e87

---

## Notes

- **Ship of Theseus Requirement**: Archives must be BYTE-FOR-BYTE IDENTICAL
- **Not "works correctly"** - IDENTICAL
- **Not "extracts the same data"** - IDENTICAL
- **Not "compatible format"** - IDENTICAL

This is like reimplementing zlib - not approximately correct, but **EXACT**.
