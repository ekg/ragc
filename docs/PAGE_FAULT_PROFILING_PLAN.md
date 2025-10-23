# Page Fault Profiling Plan

**Date**: 2025-10-23
**Goal**: Identify sources of 120M page faults in RAGC (vs C++ AGC's 5M)
**Current Performance**: RAGC 13x slower than C++ AGC, with 244x more system time

---

## Problem Summary

After implementing DashMap and LZDiff allocation optimizations:
- **Wall time**: 79.6s (vs C++ AGC's 6.2s = 13x slower)
- **System time**: 122.1s (vs C++ AGC's 0.5s = **244x more!**)
- **Page faults**: 120.3M (vs C++ AGC's ~5M = **24x more**)
- **CPU utilization**: 232% (2.3 cores vs expected 6 cores)

**Mystery**: LZDiff allocation fixes reduced wall time 12% but **page faults unchanged** at 120.3M.

This suggests page faults come from sources we haven't identified yet.

---

## Profiling Techniques

### 1. Page Fault Flame Graphs (Recommended First)

**What**: Visualizes which code paths cause page faults (actual physical memory population)

**Advantages**:
- Low overhead (~negligible, production-ready)
- Shows actual memory usage patterns
- Can compare RAGC vs C++ AGC to identify differences

**Commands**:
```bash
# Record page faults with callstacks for 30s
perf record -e page-faults -ag -o ragc.perf.data -- \
  ./target/release/ragc create -o /tmp/test.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa

# Generate text report
perf report -i ragc.perf.data --stdio -g graph,caller > ragc_pagefaults.txt

# Generate flame graph (if FlameGraph scripts available)
perf script -i ragc.perf.data > ragc.stacks
./stackcollapse-perf.pl < ragc.stacks | ./flamegraph.pl --color=mem > ragc_pagefaults.svg
```

**Alternative - trace every minor fault** (higher overhead but complete):
```bash
# Use -c 1 to trace every event (not just sample)
perf record -e minor-faults -c 1 -ag -o ragc_trace.perf.data -- \
  ./target/release/ragc create ...

perf report -i ragc_trace.perf.data --stdio -g graph,caller
```

### 2. Live Page Fault Tracing

**What**: Shows page faults in real-time as they occur

**Commands**:
```bash
# Trace minor faults live with callstacks (max depth 7)
perf trace -F min --max-stack=7 --max-events=100 \
  ./target/release/ragc create -o /tmp/test.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa
```

Useful for quick inspection of what's causing faults early in execution.

### 3. Allocation Syscall Tracing

**What**: Traces `brk()` and `mmap()` syscalls to see virtual memory expansion

**Commands**:
```bash
# Trace brk() calls (heap growth)
perf record -e syscalls:sys_enter_brk -ag -o ragc_brk.perf.data -- \
  ./target/release/ragc create ...

# Trace mmap() calls (memory mapping)
perf record -e syscalls:sys_enter_mmap -ag -o ragc_mmap.perf.data -- \
  ./target/release/ragc create ...

perf report -i ragc_brk.perf.data --stdio -g graph,caller
```

### 4. DHAT Profiler (Rust-specific)

**What**: Rust Performance Book recommends DHAT for precise allocation site identification

**Commands**:
```bash
# Using valgrind's DHAT
valgrind --tool=dhat --dhat-out-file=ragc.dhat \
  ./target/release/ragc create -o /tmp/test.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa

# View results
dh_view.html ragc.dhat
```

**Note**: Higher overhead than perf, but very precise.

### 5. Compare with C++ AGC

**Strategy**: Profile C++ AGC the same way to identify differences

```bash
# Profile C++ AGC page faults
perf record -e page-faults -ag -o agc_cpp.perf.data -- \
  agc create -o /tmp/test_cpp.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa

perf report -i agc_cpp.perf.data --stdio -g graph,caller > agc_cpp_pagefaults.txt

# Compare flame graphs side-by-side
# Look for code paths present in RAGC but not C++ AGC
```

---

## Hypotheses to Test

### Hypothesis 1: DashMap Rehashing

**Suspect**: DashMap's internal HashMap shards may be rehashing frequently

**Evidence**:
- DashMap uses per-shard HashMaps
- Each shard rehashes independently on growth
- 120M page faults = massive allocation churn

**How to confirm**:
- Page fault flame graph will show `hashbrown::map::make_hash` or `grow` functions
- Compare DashMap entry creation patterns

**Potential fix**:
- Pre-size DashMap shards with `with_capacity_and_hasher_and_shard_amount`
- Estimate total segments and configure DashMap accordingly

### Hypothesis 2: GroupWriter Buffer Growth

**Suspect**: GroupWriter's internal buffers may grow frequently

**Evidence**:
- Each group buffers segments before writing to archive
- Unknown buffer sizing strategy

**How to confirm**:
- Page fault flame graph shows `GroupWriter::add_segment` or buffer growth

**Potential fix**:
- Pre-allocate GroupWriter buffers based on expected segment sizes

### Hypothesis 3: FASTA Parsing Allocations

**Suspect**: FASTA reading may allocate new buffers per contig

**Evidence**:
- Reader thread parses 37 files with ~thousands of contigs
- Each contig may allocate separate strings

**How to confirm**:
- Page fault flame graph shows FASTA parsing functions

**Potential fix**:
- Pool contig buffers across reads
- Use arena allocator for FASTA parsing

### Hypothesis 4: Segment Assembly Allocations

**Suspect**: Building segment data structures may allocate many small Vecs

**Evidence**:
- Each segment contains multiple fields (data, masks, etc.)
- May allocate separately per field

**How to confirm**:
- Page fault flame graph shows segment construction

**Potential fix**:
- Use single contiguous allocation for segment data
- Arena allocator for related segment fields

### Hypothesis 5: Rust Allocator Overhead

**Suspect**: Rust's default allocator may have higher page fault rate than C++ allocator

**Evidence**:
- C++ AGC: 5M page faults
- RAGC: 120M page faults
- Same algorithmic complexity

**How to confirm**:
- Try alternative allocators (jemalloc, mimalloc)
- Compare page fault counts

**Potential fix**:
```rust
// In main.rs or lib.rs
use jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;
```

---

## Execution Plan

### Step 1: Quick Page Fault Flame Graph (5 min)
```bash
cd /home/erik/ragc
cargo build --release

# Record page faults
perf record -e page-faults -ag -o /tmp/ragc_pf.perf.data -- \
  ./target/release/ragc create -o /tmp/test.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa

# Generate report
perf report -i /tmp/ragc_pf.perf.data --stdio -g graph,caller > /tmp/ragc_pagefaults.txt

# View top page fault sources
head -100 /tmp/ragc_pagefaults.txt
```

### Step 2: Analyze Top Offenders (10 min)
- Identify which functions cause most page faults
- Look for unexpected patterns
- Compare with heaptrack results

### Step 3: Compare with C++ AGC (5 min)
```bash
# Profile C++ AGC
perf record -e page-faults -ag -o /tmp/agc_cpp_pf.perf.data -- \
  agc create -o /tmp/test_cpp.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa

perf report -i /tmp/agc_cpp_pf.perf.data --stdio -g graph,caller > /tmp/agc_cpp_pagefaults.txt

# Compare top sources
diff /tmp/ragc_pagefaults.txt /tmp/agc_cpp_pagefaults.txt | head -50
```

### Step 4: Targeted Fix (varies)
Based on findings:
- If DashMap rehashing → pre-size shards
- If GroupWriter growth → pre-allocate buffers
- If FASTA parsing → pool buffers
- If allocator overhead → try jemalloc/mimalloc

### Step 5: Re-benchmark and Measure
```bash
# Re-run benchmark with fix
/usr/bin/time -v ./target/release/ragc create \
  -o /tmp/test_fixed.agc -k 21 -s 10000 -m 20 -t 6 \
  /home/erik/scrapy/yeast235_split/*.fa \
  2>&1 | tee /tmp/fixed_benchmark.log

# Check page fault reduction
grep "Minor.*page faults" /tmp/fixed_benchmark.log
```

**Success criteria**:
- Page faults reduced from 120M toward 5M (C++ AGC level)
- System time reduced from 122s toward 0.5s
- Wall time approaches C++ AGC's 6.2s

---

## Expected Outcomes

### Best Case (DashMap or GroupWriter rehashing)
- **Fix**: Pre-allocate HashMap capacities
- **Expected**: 80-90% page fault reduction → ~12-24M faults
- **Wall time**: ~30-40s (2-3x slower than C++ AGC)
- **System time**: ~20-30s

### Medium Case (Multiple sources)
- **Fix**: Combination of buffer pre-allocation and allocator tuning
- **Expected**: 50-70% page fault reduction → ~36-60M faults
- **Wall time**: ~40-60s (6-10x slower than C++ AGC)
- **System time**: ~40-60s

### Worst Case (Fundamental architectural difference)
- **Fix**: Requires architectural changes to match C++ AGC
- **Expected**: Minimal reduction from simple fixes
- **May need**: Arena allocators, memory pools, or major refactoring

---

## Tools Reference

### perf event types
- `page-faults` - All page faults (minor + major)
- `minor-faults` - Only minor page faults (most relevant)
- `major-faults` - Only major page faults (disk I/O)
- `syscalls:sys_enter_brk` - Heap growth via brk()
- `syscalls:sys_enter_mmap` - Memory mapping

### perf flags
- `-a` - System-wide collection
- `-g` - Enable call-graph (stack traces)
- `-c 1` - Set sampling period to 1 (trace every event)
- `-o FILE` - Output to specific file
- `--stdio` - Text output instead of interactive TUI

### Report options
- `--stdio` - Text mode output
- `-g graph,caller` - Show call graph with caller perspective
- `--no-children` - Don't show cumulative overhead of children

---

## Next Actions

1. **Start with Step 1**: Quick page fault flame graph
2. **Identify top offender**: Which function causes most page faults?
3. **Hypothesize root cause**: Match with hypotheses above
4. **Implement targeted fix**: Based on evidence
5. **Measure improvement**: Re-benchmark to confirm reduction

**Key insight from Brendan Gregg**: "If you are hunting leaks and have a similar application that isn't growing, then taking page fault flame graphs from each and then looking for the extra code paths can be a quick way to identify the difference."

We have C++ AGC as our reference implementation with 5M faults. Comparing flame graphs will show exactly where RAGC's extra 115M faults come from.
