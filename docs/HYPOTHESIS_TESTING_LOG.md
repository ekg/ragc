# RAGC Group Fragmentation - Hypothesis Testing Log

**Problem**: RAGC creates 1713 groups vs C++ AGC's 1589 groups (+7.8% fragmentation, +124 groups)

**Data points**:
- RAGC: 11,208 segments (fewer)
- C++ AGC: 11,647 segments (more)
- Fewer segments + more groups = segments not merging into existing groups

---

## Hypothesis 1: operator[] default insertion behavior

**Theory**: C++ AGC's `map_segments[key]` inserts default value 0 if key doesn't exist, allowing splits even when groups don't exist. RAGC's defensive check prevents these splits.

**Test**: Remove RAGC's defensive group existence check and use `.get()` with fallback instead of early return.

**Results**:
- Before: 1713 groups (+124 vs target)
- After: **3437 groups** (+1848, +116% fragmentation!)

**Conclusion**: ❌ **DISPROVEN**
- Removing defensive check made fragmentation MUCH worse
- The defensive check is CORRECT and necessary
- C++ AGC must guarantee groups exist through other means

---

## Hypothesis 2: Priority queue ordering differences

**Theory**: RAGC's BinaryHeap ordering might not match C++ AGC's multimap + rbegin() behavior, causing segments to be processed in different order.

**Investigation**:

### C++ AGC Priority Queue
- Type: `multimap<pair<size_t, size_t>, T>`
- Key: `(priority, cost)`
- Ordering: Ascending by key (pair comparison)
- Pop: `rbegin()` → largest key (highest priority, highest cost)
- Priority assignment: starts at `~0ull`, decrements after each sample

### RAGC Priority Queue
- Type: `BinaryHeap<ContigTask>`
- Ordering: Custom `Ord` implementation
  ```rust
  fn cmp(&self, other: &Self) -> Ordering {
      match self.sample_priority.cmp(&other.sample_priority) {
          Equal => self.cost.cmp(&other.cost),
          other_ord => other_ord,
      }
  }
  ```
- Pop: `pop()` → largest element (highest priority, highest cost)
- Priority assignment: starts at `i32::MAX`, decrements when new sample seen

### Comparison
| Aspect | C++ AGC | RAGC | Match? |
|--------|---------|------|--------|
| Priority for first sample | ~0ull (max u64) | i32::MAX | ✓ Different types but same semantics |
| Priority decrement | After each sample | When new sample seen | ✓ Equivalent |
| Ordering: priority | Descending (via rbegin) | Descending (max-heap) | ✓ |
| Ordering: cost | Descending (via rbegin) | Descending (max-heap) | ✓ |

**Status**: ⚠️ **LIKELY CORRECT** but needs verification

**Potential issues**:
1. File reading order might differ
2. Multi-threaded priority assignment could cause non-deterministic ordering
3. Need to verify actual processing order matches

**Next test**: Add logging to compare sample processing order between implementations.

---

## Hypothesis 3: Terminal k-mer extraction differences

**Theory**: The way terminal k-mers are extracted from segments might differ, leading to different grouping keys.

**Areas to investigate**:
1. Forward vs reverse-complement handling
2. Dir-oriented vs RC-oriented k-mer selection
3. Edge cases (segments smaller than k+1 bytes)

**Status**: 🔍 **PENDING**

---

## Hypothesis 4: Terminator map population timing

**Theory**: `map_segments_terminators` might not be fully populated when split attempts occur, causing middle splitter lookups to fail.

**Evidence needed**:
1. When are terminators added to the map?
2. Are all terminal k-mers registered before segments with those k-mers attempt splits?
3. Does processing order affect terminator availability?

**Status**: 🔍 **PENDING**

---

## Hypothesis 5: Raw group (0-15) assignment differences

**Theory**: Round-robin assignment to groups 0-15 (for segments without terminal splitters) might differ due to threading.

**Evidence needed**:
1. How does C++ AGC assign raw group IDs?
2. Does RAGC's atomic fetch_add produce same assignment?
3. Are raw groups merged differently?

**Status**: 🔍 **PENDING**

---

## Hypothesis 6: Segment splitting position calculation

**Theory**: Despite algorithm appearing identical, the actual split positions might differ due to floating point or rounding differences.

**Evidence needed**:
1. Compare split positions for specific segments
2. Check if cost calculation produces identical results
3. Verify k-mer overlap calculation matches

**Status**: 🔍 **PENDING**

---

## Hypothesis 7: End-of-contig segmentation bug

**Theory**: The end-of-contig backward k-mer search has a size restriction (`remaining_after > k`) that prevents splitters from being added at contig ends, causing RAGC to miss segments.

**Investigation**:
- Analyzed `segment.rs:196` which checks `if remaining_after > k` before splitting at splitters
- This appears to skip splitters when less than k bytes remain after the position
- C++ AGC splits at ALL splitters regardless of remaining data size

**Test**: Modified `segment.rs:196` to remove the size restriction:
```rust
// Match C++ AGC: split at ALL splitters regardless of remaining data size
// C++ AGC lines 2042-2045, 2070-2072: no size restriction on splits
if true {  // was: if remaining_after > k
```

**Results**:
- Before: 11,208 segments, 1713 groups
- After: **11,208 segments, 1713 groups** (NO CHANGE)

**Discovery**: Added logging `RAGC_END_SPLIT_ATTEMPT` and found **0 executions**!

**Conclusion**: ❌ **WRONG CODE PATH**
- The end-of-contig backward search code (lines 184-236) is **NEVER EXECUTED**
- All segments are created in the main loop (lines 124-172), which IS correct
- The bug is NOT in segmentation logic
- **The problem is in SPLITTER SELECTION**, not segmentation

**Key insight**: The main segmentation loop correctly splits at EVERY splitter occurrence. But if the wrong splitters are selected in the first place, the segmentation will be wrong.

---

## Hypothesis 8: Splitter selection differences

**Theory**: RAGC and C++ AGC select different sets of splitters from the reference genome, causing different segmentation.

**Evidence**:
- ALL 444 missing segments belong to sample AIF#2
- AIF#2 chrI: 15 segments (RAGC) vs 16 (C++) - missing exactly 21 bytes (k)
- Main segmentation loop is CORRECT - splits at every splitter
- But if certain splitters are missing from RAGC's set, segments won't be created

**Investigation**:

### Debug Output Setup

Added debug output to both implementations:

**C++ AGC** (`agc_compressor.cpp`):
- Line 778: `DEBUG_CPP_FUNCTION_CALLED` - Function entry tracking
- Line 813: `DEBUG_CPP_SPLITTER` - Normal splitter selection
- Line 838: `DEBUG_CPP_SPLITTER_RIGHTMOST` - End-of-contig splitter

**RAGC** (`splitters.rs`):
- Line 265: `DEBUG_RAGC_SPLITTER_USED` - Splitter selection (already present)
- Line 281: `DEBUG_RAGC_SPLITTER_USED_RIGHTMOST` - End-of-contig splitter (already present)

### Results

Ran both implementations on reference genome (AAA#0.fa) with debug output:

**C++ AGC**:
```bash
cd /home/erik/scrapy && /home/erik/agc/bin/agc create -k 21 -s 10000 \
  /tmp/test_cpp_splitters.agc "yeast_split_proper/AAA#0.fa" 2>&1 | \
  grep DEBUG_CPP > /tmp/cpp_full_output.log
```
- **Total splitters**: 1081
- Function calls: 17 (one per contig)
- chrI example: 3 splitters (output interleaved due to threading)

**RAGC** (Pass 2 only):
```bash
cd /home/erik/scrapy && /home/erik/ragc/target/release/ragc create \
  -o /tmp/test_ragc.agc -k 21 -s 10000 -m 20 yeast_split_proper/*.fa 2>&1 | \
  grep "DEBUG_RAGC_SPLITTER_USED" > /tmp/ragc_splitters.log
```
- **Total splitters**: 1226 (filtered to Pass 2 / AAA#0 only)

### ✅ **ROOT CAUSE IDENTIFIED**

**RAGC selects 1226 splitters vs C++ AGC's 1081 (+13.4% more splitters)**

This excess of splitters creates too many split points, leading to:
1. **Different segment boundaries** - segments split at different positions
2. **Mismatched terminal k-mers** - segment ends don't align with groups
3. **Failed group merging** - segments become singleton groups
4. **+7.8% fragmentation** (+124 groups)

### Verification

First three splitters match EXACTLY between implementations:
```
Position: 62,    k-mer: 1460648714153492480, length: 10062
Position: 10062, k-mer: 5634860872301019136, length: 10000
Position: 20102, k-mer: 6077903639865720832, length: 10040
```

But RAGC continues selecting many more splitters where C++ AGC stops.

**Conclusion**: ✅ **ROOT CAUSE FOUND**

The group fragmentation issue is caused by RAGC selecting 13.4% more splitters from the reference genome than C++ AGC. Both implementations claim to use the same algorithm:
- Start with `current_len = segment_size`
- When `current_len >= segment_size` and k-mer is a candidate, mark it as splitter
- Reset `current_len = 0` after each splitter

The discrepancy suggests a subtle difference in:
1. Candidate k-mer filtering (Pass 1)
2. Splitter selection criteria (Pass 2)
3. K-mer reset logic or length tracking

**Conclusion**: ❌ **DISPROVEN** - Splitter selection is IDENTICAL!

### Final Testing (yeast10 dataset)

Ran both implementations on same 10 samples:

**C++ AGC**:
- Singletons: 11,302,379
- Splitters: 1226
- Reference: AAA#0

**RAGC**:
- Singletons: 11,302,379
- Splitters: 1226
- Reference: AAA#0

**Archives created with identical splitters, but:**
- C++ AGC: 5.4M, 1622 groups
- RAGC: 6.9M, 1851 groups (+229 groups, +14% fragmentation)

### Root Cause Analysis using `ragc inspect`

Compared archives with `ragc inspect`:

```bash
ragc inspect /tmp/test_cpp_yeast10.agc > /tmp/cpp_inspect.txt
ragc inspect /tmp/test_ragc_yeast10.agc > /tmp/ragc_inspect.txt
```

**Key Finding**: RAGC creates hundreds of **singleton groups**:
- C++ AGC Group 22: 10 segments (1 ref + 9 deltas) ✓ Merged
- RAGC Group 22: 1 segment ✗ Singleton
- RAGC Groups 1714-1768: Many singletons that should have merged

**Conclusion**: The bug is NOT in splitter selection. Segments fail to merge into existing groups despite identical splitters. Must be:
1. Different segment boundaries (despite same splitters)
2. Wrong terminal k-mer extraction
3. Failed group lookup by terminal k-mers

---

## Hypothesis 9: Segment overlap calculation inconsistency

**Theory**: RAGC uses inconsistent overlap (k bytes vs k-1 bytes) between main loop and end-of-contig handler, causing different segment boundaries.

**Investigation**:

### C++ AGC overlap behavior
- Line 2062: `split_pos = pos + 1 - kmer_length`
- Creates **k-byte overlap** consistently everywhere
- No special case for end-of-contig

### RAGC overlap behavior (BEFORE fix)
- Main loop (line 163): `let new_start = (pos + 1).saturating_sub(k)` → **k-byte overlap** ✓
- End-of-contig (line 205): `let new_start = (pos + 1).saturating_sub(k - 1)` → **(k-1)-byte overlap** ✗

**Fix Applied**: Changed line 205 to use k-byte overlap consistently

**Test Results** (yeast10 dataset):
- Before fix: 1851 groups, 6.9M
- After fix: 1852 groups, 6.9M (+1 group, no size change)
- C++ AGC: 1622 groups, 5.4M (target)

**Conclusion**: ❌ **MINIMAL IMPACT**
- The end-of-contig code path is rarely executed (0 executions in earlier tests)
- Fixing overlap had negligible effect on fragmentation
- **Root cause is elsewhere** - likely in group key calculation or terminal k-mer extraction

---

## Next Steps

1. **Test Hypothesis 10**: Compare terminal k-mer extraction between RAGC and C++ AGC
2. **Test Hypothesis 11**: Verify group key normalization logic matches C++ AGC
3. **Test Hypothesis 12**: Add debug logging to compare actual group keys used

---

## Testing Methodology

For each hypothesis:
1. **State theory clearly** with expected behavior
2. **Identify evidence needed** to prove/disprove
3. **Design minimal test** that isolates the variable
4. **Measure impact** (group count, segment count, archive size)
5. **Document conclusion** with data

**Success criteria**: Find root cause that explains:
- Why RAGC has 124 extra groups (+7.8%)
- Why RAGC has 439 fewer segments (-3.8%)
- Why fragments are created instead of merging into existing groups
