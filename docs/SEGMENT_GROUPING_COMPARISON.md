# Segment Grouping: C++ AGC vs RAGC Detailed Comparison

**Date**: 2025-12-04
**Purpose**: Identify exact differences in segment grouping logic between C++ AGC and RAGC

---

## Executive Summary

Both implementations use the same **conceptual approach**:
1. Classify segments as NEW (unknown k-mer pair) or EXISTING (known k-mer pair)
2. Buffer segments during contig processing
3. At synchronization points (sample boundaries), assign group IDs to NEW segments
4. Process segments and write to archive

However, there are **subtle differences in execution order** that may cause divergence.

---

## CRITICAL DIFFERENCE 1: Group ID Assignment Timing

### C++ AGC (agc_compressor.h, lines 384-415)

**When**: During `register_segments()` → `buffered_seg_part.process_new()` at synchronization barrier

**How**:
```cpp
// Lines 388-397
map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;  // LOCAL map for this batch
uint32_t group_id = (uint32_t)vl_seg_part.size();  // Start from current group count

// Iterate NEW segments (in s_seg_part set order)
for (const auto& x : s_seg_part)
{
    auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));

    if (p == m_kmers.end())
        m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;  // Assign new ID
}
```

**Order**: Iterates `s_seg_part` which is a **`set<kk_seg_part_t>`**. The set is ordered by the `operator<` defined in agc_compressor.h lines 157-164:

```cpp
bool operator<(const struct kk_seg_part_t& x) const
{
    if (sample_name != x.sample_name)
        return sample_name < x.sample_name;
    if (contig_name != x.contig_name)
        return contig_name < x.contig_name;
    return seg_part_no < x.seg_part_no;
}
```

**Result**: Group IDs assigned in **lexicographic order** of (sample_name, contig_name, seg_part_no) among NEW segments within the batch.

### RAGC (agc_compressor.rs, lines 2505-2680)

**When**: During `flush_batch()` at synchronization point

**How**:
```rust
// Lines 2523-2543
let mut pending = pending_batch_segments.lock().unwrap();

// CRITICAL: Sort pending segments by (sample, contig, place) before assigning group_ids
// This matches C++ AGC's BTreeSet iteration order
pending.sort();

// Lines 2552-2582
for pend in pending.iter() {
    let group_id = if pend.key.kmer_back == MISSING_KMER && pend.key.kmer_front == MISSING_KMER {
        0  // Single raw group for orphan segments
    } else {
        // Check if this k-mer pair already has a group assigned
        let mut global_map = map_segments.lock().unwrap();
        if let Some(&existing_group_id) = global_map.get(&pend.key) {
            existing_group_id  // Reuse existing group
        } else {
            // Create new group
            let new_group_id = group_counter.fetch_add(1, Ordering::SeqCst);
            global_map.insert(pend.key.clone(), new_group_id);
            new_group_id
        }
    }
}
```

**Order**: `pending.sort()` uses `PendingSegment::cmp()` defined in agc_compressor.rs lines 243-266:

```rust
fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    // First compare by sample_priority (DESCENDING - higher priority first)
    match other.sample_priority.cmp(&self.sample_priority) {
        std::cmp::Ordering::Equal => {
            // Then by sample_name
            match self.sample_name.cmp(&other.sample_name) {
                std::cmp::Ordering::Equal => {
                    // Then by contig_name
                    match self.contig_name.cmp(&other.contig_name) {
                        std::cmp::Ordering::Equal => {
                            // Finally by place (seg_part_no)
                            self.place.cmp(&other.place)
                        }
                        other => other,
                    }
                }
                other => other,
            }
        }
        other => other,
    }
}
```

**Result**: Group IDs assigned in order of:
1. **sample_priority** (DESCENDING - higher first)
2. sample_name (ascending)
3. contig_name (ascending)
4. place/seg_part_no (ascending)

---

## CRITICAL DIFFERENCE 2: Sample Priority Handling

### C++ AGC

**No explicit sample_priority field in segments.**

Segments are ordered purely by (sample_name, contig_name, seg_part_no) in the set.

When multiple samples are in the same batch, segments are **interleaved** in lexicographic order of sample_name.

### RAGC

**Explicit `sample_priority` field** (agc_compressor.rs, line 232):

```rust
struct PendingSegment {
    sample_priority: i32, // Sample processing order (higher = earlier)
    ...
}
```

Segments from different samples are **separated by priority** before comparing sample names.

**Assignment** (agc_compressor.rs, line 2798):
```rust
let sample_priority = -(num_samples_so_far as i32);  // First sample = 0, second = -1, etc.
```

**Effect**: Within the same batch (between synchronization points), segments from earlier samples get **lower group IDs** than segments from later samples, even if later sample has earlier lexicographic name.

---

## CRITICAL DIFFERENCE 3: Map Update Timing

### C++ AGC (agc_compressor.cpp, lines 1021-1045)

**When**: During `store_segments()` AFTER `register_segments()` has assigned group IDs

```cpp
// Line 1025-1031
seg_map_mtx.lock();

auto p = map_segments.find(make_pair(kmer1, kmer2));
if (p == map_segments.end())
    map_segments[make_pair(kmer1, kmer2)] = group_id;
else if (p->second > group_id)
    p->second = group_id;  // Keep LOWER group ID

// Update terminators
if (kmer1 != ~0ull && kmer2 != ~0ull)
{
    map_segments_terminators[kmer1].push_back(kmer2);
    // ...
}

seg_map_mtx.unlock();
```

**Key Point**: `map_segments` is updated during **segment writing**, not during group ID assignment. Multiple threads write simultaneously, so there's a lock and the **"keep lower group ID"** logic.

### RAGC (agc_compressor.rs, lines 2560-2582)

**When**: During `flush_batch()` immediately when assigning group ID

```rust
// Lines 2560-2581
let mut global_map = map_segments.lock().unwrap();
if let Some(&existing_group_id) = global_map.get(&pend.key) {
    existing_group_id
} else {
    // Create new group
    let new_group_id = group_counter.fetch_add(1, Ordering::SeqCst);
    global_map.insert(pend.key.clone(), new_group_id);  // Insert immediately
    drop(global_map);
    new_group_id
}
```

**Key Point**: `map_segments` is updated **atomically during group ID assignment**, not during writing. No need for "keep lower" logic because assignment happens in deterministic sorted order.

---

## CRITICAL DIFFERENCE 4: Terminators Update Timing

### C++ AGC

**When**: During `store_segments()` at **segment writing time**

**Where**: agc_compressor.cpp, lines 1033-1043

```cpp
if (kmer1 != ~0ull && kmer2 != ~0ull)
{
    map_segments_terminators[kmer1].push_back(kmer2);
    sort(map_segments_terminators[kmer1].begin(), map_segments_terminators[kmer1].end());

    if (kmer1 != kmer2)
    {
        map_segments_terminators[kmer2].push_back(kmer1);
        sort(map_segments_terminators[kmer2].begin(), map_segments_terminators[kmer2].end());
    }
}
```

**Effect**: Terminators become visible **during segment writing**, which happens **after** all group IDs are assigned for the batch.

### RAGC

**When**: During `flush_batch()` at **group ID assignment time**, but merged to global later

**Where**: agc_compressor.rs, lines 2617-2631 (batch-local) and 2696-2707 (merge to global)

```rust
// Lines 2617-2631: Update batch-local terminators
if pend.key.kmer_front != MISSING_KMER && pend.key.kmer_back != MISSING_KMER {
    let mut term_map = batch_local_terminators.lock().unwrap();

    term_map.entry(pend.key.kmer_front)
        .or_insert_with(Vec::new)
        .push(pend.key.kmer_back);

    if pend.key.kmer_front != pend.key.kmer_back {
        term_map.entry(pend.key.kmer_back)
            .or_insert_with(Vec::new)
            .push(pend.key.kmer_front);
    }
}

// Lines 2696-2707: Merge to global
{
    let batch_terms = batch_local_terminators.lock().unwrap();
    let mut global_terms = map_segments_terminators.lock().unwrap();
    for (kmer, connections) in batch_terms.iter() {
        let entry = global_terms.entry(*kmer).or_insert_with(Vec::new);
        entry.extend(connections.iter().cloned());
        entry.sort_unstable();
        entry.dedup();
    }
}
```

**Effect**: Terminators become visible **after all segments in batch are assigned group IDs**, similar timing to C++ AGC.

---

## CRITICAL DIFFERENCE 5: Processing Order

### C++ AGC (agc_compressor.cpp, lines 2235-2388)

**Queue**: `CBoundedPQueue<task_t>` - priority queue with `(priority, cost)` key

**Insertion** (lines 2322-2326):
```cpp
size_t sample_priority = ~0ull;  // Start at maximum (18446744073709551615)

pq_contigs_desc->Emplace(
    make_tuple(contig_processing_stage_t::all_contigs, sf.first, id, move(contig)),
    sample_priority, cost);

// After each sample:
--sample_priority;  // Decrease for next sample
```

**Extraction** (queue.h, lines 284-313):
```cpp
data.swap(q.rbegin()->second);  // Get LARGEST element (last in multimap)
```

**Order**: Contigs processed in order of:
1. **Higher priority** value (first sample has ~0ull, second has ~0ull - 1, etc.)
2. **Larger cost** (contig size) within same priority
3. **Multimap iteration order** for same (priority, cost)

### RAGC (agc_compressor.rs, lines 2743-2895)

**Queue**: `BinaryHeap<ContigTask>` - max-heap based on `ContigTask::cmp()`

**Insertion** (lines 2791-2820):
```rust
let sample_priority = -(num_samples_so_far as i32);  // First sample = 0, second = -1, etc.

pq.push(ContigTask {
    stage: ContigProcessingStage::AllContigs,
    sample_name: sample_name.to_string(),
    contig_name: id.clone(),
    contig: Some(contig),
    cost: contig_size,
    sample_priority,
    sequence: sequence_counter,  // FASTA order within sample
});
```

**Extraction**: Standard BinaryHeap pop() gets maximum element

**Order** (agc_compressor.rs, lines 194-212):
```rust
fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    // First by priority (DESCENDING - higher priority first)
    match self.sample_priority.cmp(&other.sample_priority) {
        std::cmp::Ordering::Equal => {
            // Then by cost (DESCENDING - larger contigs first)
            match self.cost.cmp(&other.cost) {
                std::cmp::Ordering::Equal => {
                    // CRITICAL TIE-BREAKER: Use FASTA order (sequence field)
                    // LOWER sequence = earlier in FASTA = processed first (reverse comparison for max-heap)
                    other.sequence.cmp(&self.sequence)
                }
                cost_ord => cost_ord,
            }
        }
        priority_ord => priority_ord,
    }
}
```

**Order**: Contigs processed in order of:
1. **Higher sample_priority** (first sample = 0, second = -1; 0 > -1 so first sample first)
2. **Larger cost** (contig size) within same priority
3. **Lower sequence number** (earlier in FASTA) for same (priority, cost)

---

## CRITICAL DIFFERENCE 6: Group 0 Handling

### C++ AGC (agc_compressor.cpp, line 358)

**Initialization**:
```cpp
map_segments[make_pair(~0ull, ~0ull)] = 0;
```

Group 0 is **pre-initialized** for orphan segments (both k-mers MISSING).

**Raw groups**: First `no_raw_groups` (default 16) groups use raw encoding.

**Effect**: Group 0 is ALWAYS for orphans. Groups 1-15 may be raw or LZ depending on what comes first.

### RAGC (agc_compressor.rs, lines 2556-2594)

**Assignment**:
```rust
let group_id = if pend.key.kmer_back == MISSING_KMER && pend.key.kmer_front == MISSING_KMER {
    0  // Single raw group for orphan segments
} else {
    // Regular group ID assignment
}
```

Group 0 is **assigned on-demand** when first orphan segment is encountered.

**Raw groups**: `is_raw_group = (group_id == 0)` - ONLY group 0 is raw.

**Effect**: Group 0 is for orphans, but other groups are NEVER raw (except group 0).

**MISMATCH**: C++ AGC has 16 raw groups (0-15), RAGC has only 1 (group 0).

---

## CRITICAL DIFFERENCE 7: Segment Classification

### C++ AGC (agc_compressor.cpp, lines 1552-1565)

**Classification during `add_segment()`**:

```cpp
if (p == map_segments.end())
{
    // NEW: Add to s_seg_part (unordered_set)
    buffered_seg_part.add_new(pk.first, pk.second, sample_name, contig_name,
        store_rc ? segment_rc : segment, store_rc, seg_part_no);
}
else
{
    // EXISTING: Add to vl_seg_part[segment_id]
    segment_id = p->second;
    buffered_seg_part.add_known(segment_id, ~0ull, ~0ull, sample_name, contig_name,
        store_rc ? move(segment_rc) : move(segment), store_rc, seg_part_no);
}
```

**Lookup**: `map_segments.find(pk)` with **shared_mutex** (concurrent reads)

**Storage**:
- NEW → `s_seg_part` (set)
- EXISTING → `vl_seg_part[group_id]` (vector of vectors)

### RAGC (agc_compressor.rs, lines 3799-3830)

**Classification during segment processing**:

```rust
let group_exists = {
    let global_map = map_segments.lock().unwrap();
    if global_map.contains_key(&key) {
        true
    } else {
        let batch_map = batch_local_groups.lock().unwrap();
        batch_map.contains_key(&key)
    }
};

if !group_exists {
    // NEW GROUP: Defer to pending_batch_segments
    pending_batch_segments.lock().unwrap().push(PendingSegment { ... });
    continue;
}

// EXISTING GROUP: Process immediately
```

**Lookup**: Check **both** `map_segments` (global) AND `batch_local_groups` (batch-local)

**Storage**:
- NEW → `pending_batch_segments` (Vec)
- EXISTING → `segment_groups` (BTreeMap of SegmentGroupBuffer)

**Key Difference**: RAGC checks BOTH global and batch-local maps, C++ AGC only checks global map.

---

## SCENARIOS THAT MAY CAUSE DIVERGENCE

### Scenario 1: Multi-sample Batch

**Setup**: 2 samples processed between synchronization points (in concatenated mode with large pack_cardinality)

**C++ AGC**:
- Sample A contigs → add to queue with priority ~0ull
- Sample B contigs → add to queue with priority ~0ull
- Workers pop in (priority, cost, multimap-order)
- Segments from both samples interleaved in `s_seg_part` set
- Group IDs assigned in lexicographic order of (sample_name, contig_name, seg_part_no)

**RAGC**:
- Sample A contigs → add to queue with priority 0
- Sample B contigs → add to queue with priority -1
- Workers pop Sample A first (priority 0 > -1)
- Segments separated by sample_priority before sample_name
- Group IDs assigned: all Sample A segments, then all Sample B segments

**Result**: **Different group ID assignment order** if sample names are not in priority order.

**Example**:
- Sample A name = "ZZZ"
- Sample B name = "AAA"

C++ AGC order: AAA segments, then ZZZ segments (lexicographic)
RAGC order: ZZZ segments, then AAA segments (priority)

**Impact**: Different group IDs, potentially different segment splitting decisions in subsequent samples.

### Scenario 2: Same (kmer1, kmer2) in Different Samples

**Setup**: Sample A and Sample B both have a segment with same (kmer1, kmer2) in same batch

**C++ AGC**:
- Both go into `s_seg_part` set
- Sorted by (sample_name, contig_name, seg_part_no)
- First in sorted order gets new group ID
- Second reuses same group ID (added in `process_new()`)

**RAGC**:
- Both go into `pending_batch_segments`
- Sorted by (sample_priority, sample_name, contig_name, place)
- First in sorted order gets new group ID
- Second checks `map_segments` but key already inserted by first → reuses same group ID

**Result**: **Should be the same** if order is same, but if sample_priority causes different order, may diverge.

### Scenario 3: Contig Size Tie-breaker

**Setup**: Two contigs from same sample, same size

**C++ AGC**:
- Both inserted with same (priority, cost) into multimap
- Order determined by **multimap insertion order** (undefined in spec, but typically insertion order)

**RAGC**:
- Both inserted with same (sample_priority, cost)
- Order determined by **sequence field** (FASTA order)

**Result**: **Deterministic in RAGC, non-deterministic in C++ AGC**

This was **fixed in RAGC** (agc_compressor.rs lines 203-206) but may not match C++ AGC if C++ AGC has different insertion order.

---

## RECOMMENDATIONS FOR DEBUGGING

### 1. Check Sample Priority Ordering

**Test**: Create archive with 2 samples: "ZZZ#0" then "AAA#0"

**Expectation**:
- C++ AGC: AAA segments get lower group IDs (lexicographic)
- RAGC: ZZZ segments get lower group IDs (priority-based)

**Command**:
```bash
# Create with C++ AGC
/home/erik/agc/bin/agc create -o cpp.agc -k 21 -s 10000 -l 20 -t 1 zzz.fa aaa.fa

# Create with RAGC
./target/release/ragc create -o ragc.agc -k 21 -s 10000 -m 20 -t 1 zzz.fa aaa.fa

# Compare segment layouts
./target/release/ragc inspect cpp.agc --segment-layout > cpp_layout.csv
./target/release/ragc inspect ragc.agc --segment-layout > ragc_layout.csv

diff cpp_layout.csv ragc_layout.csv
```

### 2. Check Raw Group Count

**Test**: Verify which groups are raw-encoded

**Expectation**:
- C++ AGC: Groups 0-15 are raw (no_raw_groups = 16)
- RAGC: Only group 0 is raw

**Impact**: Groups 1-15 will have different encoding (raw vs LZ) causing size differences.

### 3. Log Group ID Assignment Order

**C++ AGC**: Add logging in `process_new()` (agc_compressor.h, line 392):
```cpp
for (const auto& x : s_seg_part)
{
    auto p = m_kmers.find(make_pair(x.kmer1, x.kmer2));
    if (p == m_kmers.end())
    {
        cerr << "CPP_GROUP_ASSIGN," << x.sample_name << "," << x.contig_name
             << "," << x.seg_part_no << "," << x.kmer1 << "," << x.kmer2
             << "," << group_id << endl;
        m_kmers[make_pair(x.kmer1, x.kmer2)] = group_id++;
    }
}
```

**RAGC**: Use `RAGC_GROUP_LOG=1` environment variable (already implemented, line 3886)

**Compare**: Sort both logs by (sample, contig, place) and check if group IDs match.

### 4. Verify Priority Queue Order

**Test**: Log contig processing order

**C++ AGC**: Add logging when popping from queue (agc_compressor.cpp, line 1126):
```cpp
auto q_res = pq_contigs_desc_working->PopLarge(task);
cerr << "CPP_POP," << get<1>(task) << "," << get<2>(task) << "," << priority << endl;
```

**RAGC**: Use existing logging (agc_compressor.rs, line 2946)

**Compare**: Verify same order.

---

## CONCLUSION

The most likely source of divergence is **Scenario 1: Sample Priority Handling**.

C++ AGC sorts segments **lexicographically** by sample name within a batch.
RAGC sorts segments **by sample priority first**, then lexicographically.

If samples are processed in non-alphabetical order (e.g., "ZZZ" then "AAA"), the group ID assignment order will differ, causing:
1. Different group IDs for same k-mer pairs
2. Different segment splitting decisions (which depend on existing groups)
3. Cascading differences in subsequent samples

**FIX**: RAGC should **remove** the `sample_priority` comparison from `PendingSegment::cmp()` to match C++ AGC's pure lexicographic ordering.

**Secondary issue**: RAGC has only 1 raw group (group 0), C++ AGC has 16. This will cause size differences even if grouping order is fixed.

**FIX**: RAGC should use `NO_RAW_GROUPS = 16` and check `is_raw_group = (group_id < 16)` instead of `(group_id == 0)`.
