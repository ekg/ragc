# Next Steps: Archive Bloat Fix

**Date**: 2025-11-12
**Status**: Analysis complete, ready for targeted investigation

---

## Current State

**CRITICAL DISCOVERY (2025-11-12)**: The problem is NOT about group count!

**Actual measurements (single-threaded, yeast235 chrV):**
```
C++ AGC:
- GROUP_CREATE events: 3,660
- Archive size: 69M
- Unique k-mer pairs: 3,660

RAGC:
- Total groups: 3,381
- Archive size: 80M
- Unique k-mer pairs: 3,089

Paradox: C++ AGC has MORE groups but SMALLER archive!
```

**K-mer pair comparison:**
```
In both implementations: 2,921 (79.8%)
Only in C++ AGC:         739 (20.2% of C++ AGC total)
Only in RAGC:            168 (5.4% of RAGC total)
```

**Key findings from investigation:**
1. ✅ Threading eliminated as root cause
2. ✅ Degenerate split logic verified correct (commit 001a64c)
3. ✅ LZDiff preparation verified correct (commit e273e8e)
4. ✅ Both implementations skip creating intermediate groups for segments that will be split
5. ✅ **ROOT CAUSE IDENTIFIED**: Different segment grouping decisions, not group count
6. ❓ **WHY** the grouping differs - need to trace specific examples

---

## Immediate Next Actions

### 1. Extract Complete Group Lists from Both Archives

**Challenge**: RAGC's `inspect` command crashes when reading C++ AGC archives (index out of bounds at collection.rs:713).

**Options:**
- **A.** Fix the inspect crash to enable reading C++ AGC archives
- **B.** Create a minimal tool that directly reads AGC metadata to extract group k-mer pairs
- **C.** Add comprehensive GROUP_CREATE logging to both implementations and compare logs

**Recommended: Option C** - Most reliable, gives us exact creation order and context.

**Implementation:**
```bash
# RAGC: Already has GROUP_CREATE logging in place
cd /home/erik/scrapy && /home/erik/ragc/target/release/ragc create \
  -o /tmp/ragc_traced.agc -k 21 -s 10000 -m 20 -t 1 -v 2 \
  yeast_split_proper/*.fa 2>&1 | grep "GROUP_CREATE" > /tmp/ragc_groups.log

# C++ AGC: Enable GROUP_CREATE logging for ALL groups (not just 21-30)
# Edit agc_compressor.cpp line 1014: Remove the "if (group_id >= 21 && group_id <= 30)" condition
# Then rebuild and run with same parameters
```

### 2. Diff the Group Lists

Once we have complete lists:
```bash
# Extract k-mer pairs
grep "GROUP_CREATE" /tmp/ragc_groups.log | sed 's/.*kmers=(\([^)]*\)).*/\1/' | sort > /tmp/ragc_kmers.txt
grep "GROUP_CREATE" /tmp/cpp_groups.log | sed 's/.*kmers=(\([^)]*\)).*/\1/' | sort > /tmp/cpp_kmers.txt

# Find differences
comm -23 /tmp/ragc_kmers.txt /tmp/cpp_kmers.txt > /tmp/only_ragc.txt  # +161 expected
comm -13 /tmp/ragc_kmers.txt /tmp/cpp_kmers.txt > /tmp/only_cpp.txt   # Should be 0 or small
comm -12 /tmp/ragc_kmers.txt /tmp/cpp_kmers.txt > /tmp/both.txt       # ~3,059 expected
```

### 3. Trace Back Creation Context for Extra Groups

For each k-mer pair in `only_ragc.txt`:
1. Find the GROUP_CREATE event in full RAGC logs
2. Identify which sample/contig triggered the creation
3. Check if C++ AGC handled the same segment differently (split it? reused existing group?)
4. Look for patterns:
   - Are extra groups mostly singletons (1 segment)?
   - Are they concentrated in specific samples?
   - Are they related to specific k-mer characteristics (e.g., low complexity)?

### 4. Identify the Divergence Pattern

Analyze the extra groups to determine:
- **Hypothesis 1**: RAGC fails some Phase 1 checks that C++ AGC succeeds at
- **Hypothesis 2**: RAGC fails some Phase 2 splits that C++ AGC succeeds at
- **Hypothesis 3**: RAGC handles edge cases (empty k-mers, RC orientation) differently

### 5. Implement Targeted Fix

Once the pattern is clear:
1. Write a minimal test case that reproduces the issue
2. Implement the fix to match C++ AGC behavior
3. Verify:
   - Group count matches C++ AGC (~3,220)
   - Archive size matches C++ AGC (~69M)
   - Correctness: All samples extract byte-identical

---

## Tools to Create

### Tool 1: C++ AGC Group Logger

**File**: `/home/erik/agc/src/core/agc_compressor.cpp`
**Line**: 1014

Change:
```cpp
// OLD (only logs groups 21-30):
if (group_id >= 21 && group_id <= 30) {
    cerr << "[CPP] GROUP_CREATE: id=" << group_id
         << " kmers=(" << hex << kmer1 << "," << kmer2 << dec << ")" << endl;
}
```

To:
```cpp
// NEW (logs all groups):
cerr << "[CPP] GROUP_CREATE: id=" << group_id
     << " kmers=(" << hex << kmer1 << "," << kmer2 << dec << ")" << endl;
```

Then rebuild C++ AGC.

### Tool 2: Group Comparison Script

Create `/tmp/compare_groups.py`:
```python
#!/usr/bin/env python3
import sys, re
from collections import defaultdict

def parse_log(filename):
    """Extract group_id -> (k_front, k_back) mapping from log."""
    groups = {}
    with open(filename) as f:
        for line in f:
            m = re.search(r'GROUP_CREATE: id=(\d+).*kmers=\(([0-9a-fx]+),([0-9a-fx]+)\)', line)
            if m:
                gid = int(m.group(1))
                k1 = m.group(2)
                k2 = m.group(3)
                groups[gid] = (k1, k2)
    return groups

ragc_groups = parse_log(sys.argv[1])
cpp_groups = parse_log(sys.argv[2])

print(f"RAGC: {len(ragc_groups)} groups")
print(f"C++ AGC: {len(cpp_groups)} groups")
print(f"Delta: {len(ragc_groups) - len(cpp_groups)} extra RAGC groups")

# Find k-mer pairs unique to each
ragc_kmers = set(ragc_groups.values())
cpp_kmers = set(cpp_groups.values())

only_ragc = ragc_kmers - cpp_kmers
only_cpp = cpp_kmers - ragc_kmers

print(f"\nUnique k-mer pairs:")
print(f"  Only in RAGC: {len(only_ragc)}")
print(f"  Only in C++ AGC: {len(only_cpp)}")
print(f"  In both: {len(ragc_kmers & cpp_kmers)}")

if only_ragc:
    print(f"\nFirst 10 extra RAGC k-mer pairs:")
    for kmer_pair in list(only_ragc)[:10]:
        print(f"  {kmer_pair}")
        # Find which group ID(s) have this k-mer pair
        ragc_gids = [gid for gid, kp in ragc_groups.items() if kp == kmer_pair]
        print(f"    RAGC group IDs: {ragc_gids}")
```

---

## Timeline Estimate

- **Tool creation**: 30 minutes (modify C++ AGC, create comparison script)
- **Data collection**: 10 minutes (run both implementations with logging)
- **Analysis**: 30-60 minutes (diff, trace back, identify pattern)
- **Fix implementation**: 1-3 hours (depending on complexity)
- **Testing**: 30 minutes (verify group count, size, correctness)

**Total**: 3-5 hours to complete fix

---

## Mapping and Cost Parity Updates (in this branch)

- Split mapping: use `split_pos = best_pos` to align with C++ midpoint mapping and eliminate paired ±1 segment length drift.
- Cost index: guard linear-probing table sizing from `ht_size=0` to avoid debug underflow during tests; no behavioral change for real references.
- Action: after pulling these changes, rebuild (`cargo build --release`) and re-run layout comparisons with your C++ baseline.

## Success Criteria

1. ✅ RAGC group count matches C++ AGC (±5 groups acceptable)
2. ✅ RAGC archive size matches C++ AGC (±2% acceptable)
3. ✅ All samples extract byte-identical
4. ✅ No performance regression
5. ✅ Documentation updated with root cause and fix

---

## Alternative: Quick Win Investigation

If the systematic approach takes too long, try these quick hypotheses first:

### Quick Check 1: Singleton Groups

```bash
# Count groups with only 1 segment
/home/erik/ragc/target/release/ragc inspect /tmp/ragc_baseline.agc 2>&1 | \
  grep "Segments: 1$" | wc -l
```

If RAGC has significantly more singleton groups, it suggests Phase 1 (exact match) checks are failing when they shouldn't.

### Quick Check 2: Phase 1 Logging

Add logging to RAGC's Phase 1 check to see how often it fails to find exact matches:
```rust
// In streaming_compressor_queue.rs, Phase 1 check
if map_segments.lock().unwrap().contains_key(&group_key) {
    // Exact match found
    eprintln!("[RAGC] PHASE1_HIT: {:?}", group_key);
} else {
    eprintln!("[RAGC] PHASE1_MISS: {:?} - trying Phase 2", group_key);
}
```

Compare PHASE1_MISS counts between implementations.

### Quick Check 3: Failed Split Count

Add logging to count how many segments fail Phase 2 (splitting) and fall through to Phase 3 (new group):
```rust
// After try_split_segment_with_cost returns None
eprintln!("[RAGC] PHASE2_FAIL: {:?} - creating new group", group_key);
```

If RAGC has +161 more PHASE2_FAIL events than C++ AGC, that's the smoking gun.

---

## Files Modified During Investigation

- `BUG_ANALYSIS.md` - Comprehensive analysis and findings
- `DEGENERATE_SPLIT_BUG.md` - Previous (outdated) hypothesis
- `GROUPING_ALGORITHM_ANALYSIS.md` - Mathematical analysis
- `SPLIT_INVESTIGATION_PLAN.md` - Investigation methodology
- Logs in `/tmp/`: Various test runs and traces

---

## Notes for Future Work

1. **RAGC inspect crash**: Fix the index out of bounds error when reading C++ AGC archives (collection.rs:713)
2. **Logging consistency**: Ensure both implementations log the same events for easier comparison
3. **Test framework**: Create automated tests that compare RAGC and C++ AGC outputs
4. **Performance**: Once correctness is achieved, profile and optimize the hot paths

---

## Contact Points for Questions

- Root cause analysis: See `BUG_ANALYSIS.md`
- Algorithm description: See `GROUPING_ALGORITHM_ANALYSIS.md`
- Test methodology: See `SPLIT_INVESTIGATION_PLAN.md`
- Commit history: `git log --oneline --grep="split\|degenerate\|LZDiff"`
