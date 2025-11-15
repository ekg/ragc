# Archive Bloat Bug Analysis - REVISED

## Summary

**Problem**: RAGC creates 161 extra groups (3,381 vs 3,220) compared to C++ AGC, resulting in 16% archive bloat (80M vs 69M).

**Previous hypothesis (INCORRECT)**: Degenerate splits fail to add segments when groups aren't in thread-local buffer.

**Actual finding**: The implementations differ in **when** they check Phase 1 vs Phase 2. Need further investigation to pinpoint exact divergence.

---

## Key Discovery from Log Comparison

### C++ AGC Behavior

From `/tmp/cpp_detailed.log`:
```
[CPP] SPLIT_ACCEPT: from_kmers=(43dff13330000000,143e2e02fd400000) to_left=(1753fa9023400000,43dff13330000000) to_right=(143e2e02fd400000,1753fa9023400000) pos=10011 mid=1753fa9023400000
```

- C++ AGC encounters k-mer pair `(0x43dff13330000000, 0x143e2e02fd400000)`
- Successfully SPLITS it 7+ times across different samples
- Splits into:
  - Left: `(0x1753fa9023400000, 0x43dff13330000000)`
  - Right: `(0x143e2e02fd400000, 0x1753fa9023400000)`
- Middle k-mer: `0x1753fa9023400000`

### RAGC Behavior

From `/tmp/ragc_full_splits.log`:
```
[RAGC] GROUP_CREATE: id=26 kmers=(0x43dff13330000000,0x8f3020f51d400000)
[RAGC] GROUP_CREATE: id=27 kmers=(0x1753fa9023400000,0x43dff13330000000)
[RAGC] GROUP_CREATE: id=28 kmers=(0x143e2e02fd400000,0x1753fa9023400000)
```

- RAGC **NEVER creates** group `(0x43dff13330000000, 0x143e2e02fd400000)`
- Instead creates groups 27 and 28 directly (the split result!)
- Also creates group 26: `(0x43dff13330000000, 0x8f3020f51d400000)` ← different back k-mer

**Implication**: RAGC is skipping directly to the split result without creating the intermediate group. This is actually CORRECT behavior (both implementations should do this to avoid wasting a group ID). The question is: where do the +161 extra groups come from?

---

## Hypotheses to Investigate

### Hypothesis 1: Segmentation Difference

**Theory**: RAGC might be segmenting the genome differently than C++ AGC, creating segments with different k-mer boundaries.

**Evidence**:
- RAGC creates group `(0x43dff13330000000, 0x8f3020f51d400000)`
- C++ AGC never mentions this k-mer pair in its splits
- This suggests the segments themselves are different

**Test**: Compare splitter positions between RAGC and C++ AGC for the first sample.

### Hypothesis 2: Phase 1 Check Timing

**Theory**: C++ AGC might check Phase 1 (exact group exists) MORE AGGRESSIVELY than RAGC, reusing groups that RAGC creates as new.

**Evidence**:
- Both implementations have the same 3-phase logic
- But the ORDER of checks and the conditions might differ slightly

**Test**: Add extensive logging to RAGC's Phase 1 checks to see which groups it creates that C++ AGC would have reused.

### Hypothesis 3: Split Failure Rate

**Theory**: RAGC might be failing more splits than C++ AGC (e.g., due to missing references), causing it to create new groups instead of splitting.

**Evidence**:
- RAGC has 20,677 split attempts with 75.8% degenerate
- C++ AGC only has 2,224 SPLIT_ACCEPT events in the log
- This 10x difference is suspicious - are we logging different things?

**Test**: Count actual NEW_GROUP events in both implementations to see how many groups are created per-sample.

---

## Test Results

### Baseline (single-threaded)

**Dataset**: yeast_split_proper (235 samples, chrV only)

```
RAGC:
- Groups: 3,381
- Archive: 80M

C++ AGC:
- Groups: 3,220 (~estimate, need to verify exact count)
- Archive: 69M

Delta: +161 groups (+5%), +11M size (+16%)
```

### Previous Fix Attempt (FAILED)

**Approach**: Used `.entry().or_insert_with()` to create buffers with `register_stream()` calls.

**Result**:
- Archive: 98M (+22% WORSE than baseline!)
- Groups: 3,347 (slight improvement in count, but massive size increase)

**Root cause**: Created duplicate stream registrations for groups that already existed globally.

---

## Files Created During Investigation

- `/tmp/cpp_detailed.log` (257KB, 2,224 lines) - C++ AGC split/segment events
- `/tmp/ragc_full_splits.log` (126MB, many lines) - RAGC full logging
- `/tmp/ragc_split_decisions.log` (4.9MB, 56,278 lines) - RAGC split decisions only
- `/tmp/split_check.log` (13KB, 100 lines) - RAGC split existence checks
- `/home/erik/ragc/DEGENERATE_SPLIT_BUG.md` - Previous (outdated) analysis
- `/home/erik/ragc/GROUPING_ALGORITHM_ANALYSIS.md` - Mathematical analysis
- `/home/erik/ragc/SPLIT_INVESTIGATION_PLAN.md` - Investigation methodology

---

## CRITICAL DISCOVERY (2025-11-12)

### Comprehensive Group Comparison Results

**Actual measurements** (single-threaded, yeast235 chrV):
```
C++ AGC:
- GROUP_CREATE events logged: 3,660
- Archive size: 69M
- Unique k-mer pairs: 3,660

RAGC:
- Flushing group events logged: 3,089
- Archive size: 80M
- Total groups (from final count): 3,381
- Unique k-mer pairs: 3,089

K-mer pair comparison:
- In both: 2,921
- Only in C++ AGC: 739
- Only in RAGC: 168
```

### The Paradox

**C++ AGC creates MORE groups (3,660) but produces SMALLER archives (69M)!**

This reveals the ROOT CAUSE is NOT about group count - it's about **SEGMENT GROUPING DECISIONS**.

### Key Findings

1. **Different group decisions**: Only 2,921 k-mer pairs (79.8%) are shared between implementations
2. **C++ AGC's extra groups**: 739 unique k-mer pairs (20.2% of its total)
3. **RAGC's extra groups**: 168 unique k-mer pairs (5.4% of its total)
4. **Size vs count paradox**: More groups (C++ AGC) → smaller archive → **better compression**

### Root Cause Hypothesis

**C++ AGC creates MORE groups but these groups contain segments that compress BETTER together.**

Possible reasons:
1. **Different Phase 1/2 check timing**: C++ AGC might be more aggressive at splitting, creating specialized groups
2. **Different splitting decisions**: Same segment might split into different sub-groups
3. **Different segment boundaries**: Splitter detection might produce different boundaries
4. **Group reuse logic**: The `else if (p->second > group_id) p->second = group_id;` path in C++ AGC (lines 1017-1018) updates existing groups to prefer smaller IDs - this might affect subsequent lookups

### Next Investigation Steps

1. ✅ Compare splitter counts and positions
2. Compare segment boundaries for same contigs
3. Compare Phase 1/Phase 2 success rates
4. Analyze one specific k-mer pair that differs to understand the divergence

---

## Log Analysis Notes

### C++ AGC Log Structure

```
[CPP] SEG_ADD: gid=30 sample=AAA#0 contig=AAA#0#chrI part=14 len=10060 kmers=(16179894a7800000,342117afaf000000)
[CPP] SPLIT_ACCEPT: from_kmers=(...) to_left=(...) to_right=(...) pos=10011 mid=(...)
```

- Shows **actual k-mer values** for first sample (AAA#0)
- Shows `kmers=(ffffffffffffffff,ffffffffffffffff)` for subsequent samples (delta-encoded)
- Only logs SPLIT_ACCEPT events, not split attempts that fail or are skipped

### RAGC Log Structure

```
[RAGC] GROUP_CREATE: id=26 kmers=(0x43dff13330000000,0x8f3020f51d400000)
[RAGC] SEG_ADD: gid=26 sample=AAB#0 contig=AAB#0#chrI part=9 len=10020 kmers=(0x43dff13330000000,0x8f3020f51d400000)
[RAGC] SPLIT_SUCCESS: 4,996 events (24.2%)
[RAGC] SPLIT_DEGENERATE: 15,682 events (75.8%)
```

- Shows k-mer values for all segments (not just first sample)
- Logs both successful and degenerate splits
- Total split attempts: 20,677 (vs C++ AGC's 2,224 SPLIT_ACCEPT events)

**Question**: Why 10x difference in split counts? Are we logging at different granularities?

---

## Commit History Context

- `001a64c`: Added degenerate split support (reduced yeast10 from 6.7M → 5.7M)
- `e273e8e`: On-demand LZDiff preparation (fixed segment bloat from 285K → 11K segments)
- Current: yeast235 (chrV) produces 80M archive with 3,381 groups (vs 69M, 3,220 groups)

---

## Decision: Pause for Comprehensive Group Comparison

Before attempting another fix, we need **concrete data** on which groups differ:

1. Extract ALL group k-mer pairs from both archives
2. Create a diff showing:
   - Groups only in RAGC (+161 expected)
   - Groups only in C++ AGC (if any)
   - Groups in both (should be ~3,059)
3. For each extra RAGC group, trace back through logs to find when/why it was created
4. Identify the pattern: Is it related to specific samples? Specific segment sizes? Specific k-mer characteristics?

Only then can we implement a targeted fix that addresses the root cause.
