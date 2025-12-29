# RAGC Compression Quality Investigation - Status

## Current Status (2025-12-29)

**Archive Size Comparison (yeast235):**
- C++ AGC (1 thread): 72.0 MB, 307,708 segments, 3,678 groups
- RAGC (1 thread): 76.8 MB, 299,297 segments, 3,444 groups
- **Gap: 6.55% larger** (down from 7.2%)

## Recent Improvements

### Cost-Based Split Optimization Implemented
- Added `find_split_by_cost()` function matching C++ AGC's `get_coding_cost` approach
- Computes LZ encoding cost at every position using `get_coding_cost_vector()`
- Uses `partial_sum` for cumulative cost from both left and right
- Added `SplitDecision` enum with three outcomes:
  - `SplitAt(pos)`: Actually split the segment
  - `AssignToLeft`: Assign entire segment to left group (when split would create tiny left piece)
  - `AssignToRight`: Assign entire segment to right group (when split would create tiny right piece)

### Results
- Segment count reduced from 317,693 to 299,297 (fewer unnecessary splits)
- Tiny segments (<100 bytes) reduced from 21,857 to 3,362 (fewer than C++ AGC's 4,356!)
- Archive size reduced from 77.2 MB to 76.8 MB

## Key Findings

### 1. Barrier Splitting is Beneficial
- RAGC's barrier splitting feature helps compression
- Disabling it makes archives ~13% larger (81.6 MB)
- This is functionality C++ AGC also has (via `find_cand_segment_with_missing_middle_splitter`)

### 2. Single-Threaded Produces Better Compression
- 1 thread: 76.8 MB
- 8 threads: 80.2 MB
- Parallelism introduces non-determinism in group creation order

### 3. Cost-Based Split Logic (IMPLEMENTED)
**C++ AGC approach (`get_coding_cost`):**
- Computes LZ encoding cost at EVERY position in the segment
- Uses `partial_sum` to accumulate costs from both ends
- Finds the position with minimum combined cost (left + right)
- If best_pos is near edges, assigns entire segment to one group instead of splitting

**RAGC now matches this** (see `find_split_by_cost()` in agc_compressor.rs lines 5885-6051)

## Remaining Gap Analysis

The 6.55% gap (~4.7 MB) may be due to:
1. Different grouping decisions elsewhere in the pipeline
2. LZ encoding efficiency differences
3. Reference segment selection differences
4. Different handling in `find_cand_segment_with_one_splitter` (one k-mer missing)

## Verification Commands

```bash
# Build RAGC
cargo build --release

# Create archive with C++ AGC (baseline)
/home/erik/agc/bin/agc create -o /tmp/cpp.agc -k 21 -s 10000 -l 20 -t 1 samples/*.fa

# Create archive with RAGC
./target/release/ragc create -o /tmp/ragc.agc -k 21 -s 10000 -m 20 -t 1 samples/*.fa

# Compare sizes
ls -la /tmp/cpp.agc /tmp/ragc.agc

# Compare segment layouts
./target/release/ragc inspect /tmp/cpp.agc --segment-layout > cpp_layout.csv
./target/release/ragc inspect /tmp/ragc.agc --segment-layout > ragc_layout.csv
diff cpp_layout.csv ragc_layout.csv | head -50
```

## Test Data Location

- Yeast 235 samples: `/home/erik/scrapy/yeast235_samples/*.fa`
- Baseline C++ AGC archive: Create with `/home/erik/agc/bin/agc`

## Progress History

| Date | Gap | Notes |
|------|-----|-------|
| Before | 33% | Initial handoff document |
| 2025-12-29 | 7.2% | After fixing orphan distribution, raw groups, etc. |
| 2025-12-29 | 6.55% | After implementing cost-based split optimization |

The remaining 6.55% gap may require investigating other parts of the compression pipeline.
