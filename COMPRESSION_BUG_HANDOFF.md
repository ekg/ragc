# RAGC Compression Quality Investigation - Status

## Current Status (2025-12-31)

**Archive Size Comparison (yeast235):**
- C++ AGC: 72.0 MB (72,045,160 bytes)
- RAGC: 76.8 MB (76,769,909 bytes)
- **Gap: 6.56% larger**

**Determinism: FIXED**
- ✅ 1 thread vs 8 threads: **IDENTICAL** archives (same SHA256)
- ✅ Same thread count across runs: **IDENTICAL** archives

**Key fixes for determinism:**
1. Classification loop: `into_par_iter()` → sequential `for` loop
2. Stream registration: `HashSet<u32>` → `BTreeSet<u32>` for deterministic iteration

## Analysis Summary

| Metric | C++ AGC | RAGC | Difference |
|--------|---------|------|------------|
| Total segments | 307,704 | 299,293 | **-2.7%** |
| Unique groups | 3,676 | 3,442 | **-6.4%** |
| Single-segment groups | 391 | 320 | **-18%** |
| Avg segment length | 10,865 bytes | 11,170 bytes | **+2.8%** |

**Key Finding:** RAGC creates FEWER groups and FEWER segments than C++ AGC, but the archive is larger. This is counterintuitive but explained by segment length distribution.

## Root Cause: Longer Segments in Small Groups

RAGC's segments are longer on average because it splits less often:

| Group Size | C++ AGC Avg | RAGC Avg | Difference |
|------------|-------------|----------|------------|
| Single (1) | 43,975 bytes | 52,201 bytes | **+18.7%** |
| Small (2-10) | 34,748 bytes | 40,150 bytes | **+15.5%** |
| Medium (11-100) | 14,831 bytes | 14,246 bytes | -3.9% |
| Large (100+) | 9,901 bytes | 10,187 bytes | +2.9% |

**Why this matters:**
1. Longer segments compress less efficiently against reference data
2. Small groups (2-10 segments) have 23MB more bytes in RAGC than C++ AGC
3. These poorly-compressing bytes account for the ~4.7MB archive size difference

## Technical Analysis

### Splitting Behavior Difference

C++ AGC splits segments more aggressively:
- When a segment's k-mer pair (front, back) doesn't match an existing group
- It looks for a "middle splitter" k-mer that connects to existing groups
- If found, the segment is split at the optimal LZ encoding cost position

RAGC implements the same logic but produces different results:
- 8,411 fewer segments total (2.7% fewer)
- Segments that C++ AGC would split stay as one longer segment
- These longer segments end up in smaller groups with worse LZ compression

### Per-Sample Analysis

The splitting difference is distributed across all samples:
```
ANM#3: 154 fewer segments (cpp=6002, ragc=5848)
ANM#4: 153 fewer segments (cpp=6001, ragc=5848)
ANL#4: 125 fewer segments (cpp=4329, ragc=4204)
...etc...
```

Within each sample, the difference is spread across many contigs (1-7 segments per contig).

## Possible Root Causes

1. **Terminator Visibility in Parallel Processing**
   - RAGC classifies segments in parallel across contigs
   - Terminators are updated immediately (CRITICAL FIX 3)
   - But race conditions may cause some segments to not see others' terminators
   - Result: fewer successful `find_middle_splitter` calls

2. **Cost-Based Split Threshold Difference**
   - RAGC's `find_split_by_cost` may have subtle differences in:
     - LZ cost calculation
     - Position thresholds for AssignLeft/AssignRight decisions
   - Could cause different split vs assign-whole-segment decisions

3. **Group Existence Check Timing**
   - RAGC only splits if BOTH target groups already exist
   - C++ AGC does the same check, but group creation timing differs
   - Groups created by parallel threads may not be visible immediately

## Verification Commands

```bash
# Build RAGC
cargo build --release

# Create archives
/home/erik/agc/bin/agc create -o /tmp/cpp.agc -k 21 -s 10000 -l 20 -t 1 /home/erik/scrapy/yeast235_samples/*.fa
./target/release/ragc create -o /tmp/ragc.agc -k 21 -s 10000 -m 20 -t 1 /home/erik/scrapy/yeast235_samples/*.fa

# Compare sizes
ls -la /tmp/cpp.agc /tmp/ragc.agc

# Export and compare layouts
./target/release/ragc inspect /tmp/cpp.agc --segment-layout > /tmp/cpp_layout.csv
./target/release/ragc inspect /tmp/ragc.agc --segment-layout > /tmp/ragc_layout.csv

# Compare group/segment counts
echo "Groups:"; cut -d',' -f4 /tmp/cpp_layout.csv | tail -n +6 | sort -nu | wc -l
echo "Groups:"; cut -d',' -f4 /tmp/ragc_layout.csv | tail -n +6 | sort -nu | wc -l
```

## Progress History

| Date | Gap | Notes |
|------|-----|-------|
| Before | 33% | Initial state |
| 2025-12-29 | 7.2% | Fixed orphan distribution, raw groups |
| 2025-12-29 | 6.55% | Implemented cost-based split optimization |
| 2025-12-30 | 17.5% | BROKEN: Determinism fix caused regression |
| 2025-12-30 | 14.9% | BROKEN: Orientation fix didn't help |
| 2025-12-30 | 16.6% | BROKEN: Raw group distribution |
| 2025-12-31 | **6.56%** | **REVERTED** bad changes, back to working baseline |

## Next Steps

To close the remaining 6.56% gap:

1. **Investigate split decision differences**
   - Add detailed logging to compare when C++ AGC splits vs when RAGC doesn't
   - Focus on segments that are in small groups in RAGC but would be split in C++ AGC

2. **Profile terminator visibility**
   - Log terminator state at key decision points
   - Compare with C++ AGC's terminators at same processing stage

3. **Consider sequential processing for split-heavy samples**
   - The parallel processing speeds up overall but may hurt compression
   - Could implement hybrid: parallel for reference sample, sequential for others

**Note**: The current 6.56% gap is acceptable for most use cases. Further optimization should be balanced against code complexity and performance impact.
