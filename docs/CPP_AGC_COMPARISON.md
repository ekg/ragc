# C++ AGC vs RAGC Algorithm Comparison

**Goal**: Identify all algorithmic differences causing RAGC to produce 68 missing groups (4.2%) and 20% larger archives.

**Status**:
- Groups: 1,554 (RAGC) vs 1,622 (C++ AGC) = 68 missing
- Archive size: 6.5M vs 5.4M = 20% larger
- Extraction: ✅ Works correctly

---

## Phase 1: Input Processing & Contig Ordering

### 1.1 Contig Processing Order
- **C++ AGC**: Priority queue sorted by contig size (largest first)
  - Location: `agc_compressor.cpp:2277` - `pq_contigs_desc->Emplace(..., sample_priority, cost)` where `cost = contig.size()`
- **RAGC**: Now matches! Sorted at CLI level by size descending
  - Location: `ragc-cli/src/main.rs:363-372`
- **Status**: ✅ **IDENTICAL** (as of 2025-11-06)

### 1.2 Sample Priority
- **C++ AGC**: `sample_priority` decreases after synchronization points
  - Location: `agc_compressor.cpp:2266` - `--sample_priority`
- **RAGC**: Processes samples in input file order
- **Status**: ⚠️ **NEEDS INVESTIGATION**
  - Question: Does sample priority affect compression or just processing order?

### 1.3 K-mer Extraction from Contigs
- **C++ AGC**: Location: `agc_compressor.cpp:2064-2091` (segment_data function)
- **RAGC**: Location: `ragc-core/src/segment.rs:split_at_splitters_with_size`
- **Status**: ⏳ **NOT YET CHECKED**

---

## Phase 2: Segment Grouping (K-mer Terminators)

### 2.1 Segment Group Key Creation
- **C++ AGC**: `minmax(kmer_front, kmer_back)` to normalize key
  - Location: `agc_compressor.cpp:1361`
- **RAGC**: `SegmentGroupKey::new(front, back)` - always `(min, max)`
  - Location: `ragc-core/src/streaming_compressor_queue.rs:961-968`
- **Status**: ⏳ **NOT YET CHECKED** - need to verify normalization matches exactly

### 2.2 MISSING_KMER Sentinel Value
- **C++ AGC**: `kmer_t(-1)` = `u64::MAX` (18446744073709551615)
  - Location: `common/agc_basic.h`
- **RAGC**: `MISSING_KMER = u64::MAX`
  - Location: `ragc-core/src/segment.rs:15`
- **Status**: ✅ **IDENTICAL** (fixed 2025-11-06)

### 2.3 Reverse Complement Handling
- **C++ AGC**: Compares `kmer` vs `kmer.get_rev_comp()`, uses canonical (smaller) value
  - Location: `agc_compressor.cpp:1353-1377` (CASE1, CASE2, CASE3)
- **RAGC**: Uses `should_reverse` based on canonical k-mer
  - Location: `ragc-core/src/streaming_compressor_queue.rs:987-1011`
- **Status**: ⏳ **NOT YET CHECKED** - this is complex and could have subtle differences

### 2.4 Terminator Map Population
- **C++ AGC**: Bidirectional edges added to `map_segments_terminators`
  - Location: `agc_compressor.cpp:1319-1334`
- **RAGC**: Bidirectional connections in `map_segments_terminators`
  - Location: `ragc-core/src/streaming_compressor_queue.rs:1179-1201`
- **Status**: ⏳ **NOT YET CHECKED**

---

## Phase 3: Segment Splitting

### 3.1 Split Eligibility Check
- **C++ AGC**: Checks if group key doesn't exist AND both k-mers in terminators map
  - Location: `agc_compressor.cpp:1398-1403`
- **RAGC**: Same logic - `!key_exists` AND both k-mers not MISSING_KMER
  - Location: `ragc-core/src/streaming_compressor_queue.rs:1018`
- **Status**: ✅ **IDENTICAL**

### 3.2 Finding Middle Splitter
- **C++ AGC**: `set_intersection` of terminator lists, removes `~0ull`, returns first
  - Location: `agc_compressor.cpp:1545-1566`
- **RAGC**: Set intersection, filters MISSING_KMER, returns first
  - Location: `ragc-core/src/streaming_compressor_queue.rs:1237-1276`
- **Status**: ✅ **IDENTICAL** (verified 2025-11-06)

### 3.3 Split Position Calculation
- **C++ AGC**: **Uses compression cost heuristic** - tries encoding with both split groups, finds minimum cost position
  - Location: `agc_compressor.cpp:1573-1626`
  - Calls `seg1->get_coding_cost()` and `seg2->get_coding_cost()`
  - Finds `best_pos` with minimum sum of costs
  - Can return `best_pos = 0` or `segment.size()` if split is not beneficial
- **RAGC**: **Uses midpoint heuristic** - splits at `segment_len / 2`
  - Location: `ragc-core/src/streaming_compressor_queue.rs:1288-1290` (find_split_position)
  - Comment says: "C++ AGC uses compression cost (lines 1556-1625), but for a simplified approach we split at the midpoint"
- **Status**: ❌ **DIFFERENT - CRITICAL ISSUE**
  - **Impact**: This could explain BOTH missing groups AND size bloat
  - If C++ AGC's cost calculation returns degenerate positions (0 or segment.size()), it **doesn't split**
  - RAGC always splits at midpoint, may split segments that shouldn't be split

### 3.4 Degenerate Split Handling
- **C++ AGC**: When `best_pos < kmer_length + 1` → sets `best_pos = 0`
  - When `best_pos + kmer_length + 1 > size` → sets `best_pos = segment.size()`
  - These degenerate cases result in `left_size == 0` or `right_size == 0`
  - Location: `agc_compressor.cpp:1628-1631, 1442-1459`
- **RAGC**: No equivalent handling - always splits if middle k-mer found
  - Location: None
- **Status**: ❌ **DIFFERENT - LIKELY ROOT CAUSE**

---

## Phase 4: LZ Encoding

### 4.1 Reference Selection
- **C++ AGC**: First segment in sorted buffer
  - Location: `segment.cpp:42`
- **RAGC**: First segment in sorted buffer
  - Location: `ragc-core/src/streaming_compressor_queue.rs:637`
- **Status**: ⏳ **NOT YET CHECKED**

### 4.2 LZ Diff Algorithm
- **C++ AGC**: Uses `LZDiff` class
  - Location: `common/lz_diff.cpp`
- **RAGC**: Rust port of LZDiff
  - Location: `ragc-core/src/lz_diff.rs`
- **Status**: ⏳ **NOT YET CHECKED**

---

## Phase 5: Compression & Serialization

### 5.1 ZSTD Compression Level
- **C++ AGC**: Default level 17
- **RAGC**: Default level 17
- **Status**: ✅ **IDENTICAL**

### 5.2 Pack Cardinality (segments per pack)
- **C++ AGC**: `PACK_CARDINALITY = 128`
  - Location: `common/agc_basic.h:22`
- **RAGC**: `PACK_CARDINALITY = 128`
  - Location: `ragc-core/src/streaming_compressor_queue.rs:18`
- **Status**: ✅ **IDENTICAL**

---

## Summary of Findings

### ✅ Confirmed Identical
1. Contig processing order (size-sorted, largest first)
2. MISSING_KMER sentinel value (u64::MAX)
3. Split eligibility check logic
4. Middle splitter finding algorithm
5. ZSTD compression level
6. Pack cardinality

### ❌ Confirmed Different
1. **Split position calculation** - C++ uses compression cost, RAGC uses midpoint
2. **Degenerate split handling** - C++ can decide NOT to split, RAGC always splits

### ⏳ Not Yet Checked
1. Sample priority impact
2. K-mer extraction details
3. Segment group key normalization details
4. Reverse complement handling details
5. Terminator map population details
6. LZ encoding algorithm
7. Reference selection

---

## Next Steps

### Immediate Action: Fix Split Position Calculation

**Problem identified**:
- C++ AGC: 940 splits executed (out of 1,356 attempts = 69% success)
- RAGC: 1,044 splits executed (104 extra splits, ~11% more)

**Root cause**: Different split position algorithms
- C++ AGC uses compression cost heuristic + degenerate position handling
- RAGC uses midpoint heuristic + always splits

**Proposed solutions** (in order of complexity):

1. **Simple fix**: Add minimum segment size check
   - Don't split if it would create segments < k+1 bytes
   - Similar to C++ AGC's degenerate position rules
   - Easy to implement, may not be complete

2. **Proper fix**: Implement compression cost heuristic
   - Use `LZDiff::get_coding_cost_vector()` (already implemented in RAGC)
   - Access split groups' LZDiff structures
   - Find position with minimum cost
   - Apply degenerate position rules
   - Complex but correct

### Other Items

2. Verify reverse complement handling (complex, easy to get wrong)

3. Check LZ encoding implementation for differences

4. Verify k-mer extraction matches exactly

---

## Notes

- **Testing methodology**: 10-sample yeast dataset, k=21, s=10000
- **Current divergence**: 68 missing groups (4.2%), 20% size bloat
- **Extraction correctness**: ✅ Verified working
