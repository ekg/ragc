# RAGC Development State - 2025-11-08

## ✅ Completed Today

### Raw Groups (0-15) Support
- **Fixed**: Round-robin distribution of MISSING_KMER segments across groups 0-15
- **Fixed**: Raw groups correctly skip reference segment creation
- **Fixed**: Inspection tool correctly identifies raw vs LZ groups
- **Result**: Groups 0-15 now match C++ AGC behavior (0 refs)

### AGC Archive Inspection Tool
- **New command**: `ragc inspect <archive> [--segments] [--group-id N]`
- **Features**:
  - Show group statistics (total, refs, deltas per group)
  - Filter by specific group ID
  - Display detailed segment information
  - Works on both RAGC and C++ AGC archives

## 🔍 Ongoing Investigation: Group Fragmentation

### Problem
RAGC creates **187 more groups** than C++ AGC (11.8% more groups):
- **RAGC**: 1776 groups →  7.8 MB archive
- **C++ AGC**: 1589 groups → 6.0 MB archive
- **Impact**: 30% larger archives

### Root Cause Analysis

#### Group Size Distribution Comparison

**C++ AGC:**
- 10-segment groups: 589 (peak - segments appearing in all 10 samples)
- Single-segment groups: 259

**RAGC:**
- 10-segment groups: 337 (-252 from expected!)
- Single-segment groups: 332 (+73 more)

**Key Finding**: Segments that should merge into 10-segment groups are being fragmented across multiple smaller groups.

#### Attempted Fixes

1. ✅ **Atomic group allocation** using HashMap::entry() API
   - Prevents race between check and insert
   - Matches C++ AGC's seg_map_mtx locking pattern
   - **Result**: Still 1776 groups (no improvement)

2. ❌ **Pre-allocation approach** (reverted)
   - Wasted group IDs when keys already exist
   - Made fragmentation worse (1777 groups)

### Hypotheses for Remaining Fragmentation

The atomic fix should have prevented race conditions, yet fragmentation persists. Possible causes:

1. **Different k-mer canonicalization**
   - RAGC and C++ AGC may canonicalize k-mers differently
   - Results in different SegmentGroupKey values for same biological sequence

2. **Segment splitting differences**
   - find_middle_splitter logic may differ
   - Split position calculations may vary
   - Creates more unique keys than C++ AGC

3. **Reverse complement handling**
   - Case 3a/3b logic (single k-mer matching) may differ from C++ AGC
   - RC orientation decisions could create duplicate keys

4. **Thread execution order effects**
   - Even with atomic allocation, segment processing order differs
   - May affect splitting decisions (split only if both groups exist)

### Next Steps

1. **Compare k-mer canonicalization** between RAGC and C++ AGC
   - Test with single-threaded execution (num_threads=1)
   - Check if fragmentation persists without threading

2. **Analyze segment key distribution**
   - Dump all SegmentGroupKeys from both archives
   - Find which keys RAGC creates that C++ AGC doesn't

3. **Trace splitting logic**
   - Add detailed logging for split decisions
   - Compare split points between implementations

## Testing Dataset

**yeast10**: 10 whole yeast genomes
- Total size: ~3.3GB uncompressed
- Segments: ~11,000
- Chromosomes: Multiple per sample

## Key Files Modified

- `ragc-core/src/streaming_compressor_queue.rs`: Atomic group allocation
- `ragc-core/src/decompressor.rs`: Inspection methods
- `ragc-cli/src/inspect.rs`: New inspection command
- `ragc-cli/src/main.rs`: Inspect subcommand integration

## Commits

- `6765238`: feat: Add raw groups (0-15) support and AGC archive inspection tool
