# Batch Mode Rewrite: Achieving Byte-Identical Archives

**Date**: 2025-11-30
**Goal**: Rewrite RAGC to use C++ AGC's batch-local state management for byte-identical archives
**Status**: IN PROGRESS

---

## Problem Statement

**Current RAGC ("Streaming" Mode)**:
- Uses global persistent state that NEVER resets
- Accumulates ALL segments in memory
- Writes everything at finalize()
- Result: 12-16% larger archives than C++ AGC

**Root Cause**:
```rust
// Global state - never reset
segment_groups: Arc<Mutex<BTreeMap<Key, SegmentGroupBuffer>>>

// Identical k-mer pairs always reuse same group
groups.entry(key).or_insert_with(...).push(segment)
```

This is NOT streaming - it's batch accumulation!

---

## C++ AGC Batch Architecture

**Key Insight**: C++ AGC uses **batch-local state** that resets at sample boundaries:

```cpp
// Called at each sample boundary
uint32_t process_new() {
    map<pair<uint64_t, uint64_t>, uint32_t> m_kmers;  // LOCAL variable!

    // Assign group IDs to new segments from THIS batch only
    for (segment in current_batch) {
        if (!m_kmers.contains(segment.kmer_pair))
            m_kmers[segment.kmer_pair] = next_group_id++;
    }

    // m_kmers destroyed here - next batch starts fresh!
}
```

**Result**: Identical k-mer pairs in different samples get DIFFERENT groups.

---

## Architecture Design

### Data Structures

#### 1. Global Registry (Lookup Only)
```rust
/// Global registry of all groups ever created
/// Used ONLY for lookup - never modified during batch processing
struct GlobalRegistry {
    /// Map from k-mer pair to group ID
    kmer_to_group: HashMap<(u64, u64), u32>,

    /// All segment groups indexed by group_id
    groups: Vec<SegmentGroup>,

    /// Next available group ID
    next_group_id: u32,
}
```

#### 2. Batch-Local State (Equivalent to C++ AGC's m_kmers)
```rust
/// Batch-local state - RESET at sample boundaries
struct BatchState {
    /// New segments discovered in THIS batch (not found in global registry)
    /// Key: (front_kmer, back_kmer)
    /// Value: Vec of segments with that k-mer pair
    new_segments: BTreeMap<(u64, u64), Vec<SegmentData>>,

    /// Starting group ID for this batch
    next_group_id: u32,
}

impl BatchState {
    fn new(starting_group_id: u32) -> Self {
        BatchState {
            new_segments: BTreeMap::new(),
            next_group_id: starting_group_id,
        }
    }

    fn clear(&mut self) {
        self.new_segments.clear();
        // next_group_id continues incrementing
    }
}
```

#### 3. Segment Group (With Pack Management)
```rust
struct SegmentGroup {
    group_id: u32,
    segments: Vec<SegmentData>,
    reference_written: bool,

    /// Pack size (default: 50 segments)
    pack_size: usize,
}

impl SegmentGroup {
    fn should_flush_pack(&self) -> bool {
        self.segments.len() >= self.pack_size
    }

    fn take_segments_for_pack(&mut self) -> Vec<SegmentData> {
        std::mem::take(&mut self.segments)
    }
}
```

---

## Implementation Plan

### ✅ Phase 1: Add Batch State Infrastructure (COMPLETED)

**Files**: `ragc-core/src/streaming_compressor_queue.rs`

**Tasks**:
- [x] Add `BatchState` struct - lines 327-359
- [x] Modify `StreamingQueueConfig` to add `pack_size: usize` (default: 50) - lines 114-117, 133
- [x] Add pack-full detection methods to `SegmentGroupBuffer` - lines 314-324
  - `should_flush_pack(pack_size)` - check if >= pack_size segments
  - `segment_count()` - get current count

**Verification**: ✅ Code compiles successfully, no behavior change

**Note**: GlobalRegistry not added yet - current architecture already has global state via Arc<Mutex<BTreeMap>>. Will refactor in Phase 2 as needed.

---

### ☐ Phase 2: Implement Batch Processing Loop

**Core Algorithm**:
```rust
fn process_sample_batch(
    sample: Sample,
    global_registry: &mut GlobalRegistry,
    config: &Config,
) -> Result<()> {
    // Create batch-local state
    let mut batch_state = BatchState::new(global_registry.next_group_id);

    // Step 1: Process all contigs in this sample
    for contig in sample.contigs {
        let segments = split_contig(contig, &splitters, config.k, config.segment_size);

        for segment in segments {
            let key = (segment.front_kmer, segment.back_kmer);

            // Check if k-mer pair exists in GLOBAL registry
            if let Some(&group_id) = global_registry.kmer_to_group.get(&key) {
                // KNOWN segment - add to existing group
                global_registry.groups[group_id as usize].add_segment(segment);

                // Write pack if group is full
                if global_registry.groups[group_id as usize].should_flush_pack() {
                    flush_pack(&mut global_registry.groups[group_id as usize], config)?;
                }
            } else {
                // NEW segment - add to batch buffer
                batch_state.new_segments
                    .entry(key)
                    .or_insert_with(Vec::new)
                    .push(segment);
            }
        }
    }

    // Step 2: Process new segments from this batch
    process_batch_new_segments(&mut batch_state, global_registry, config)?;

    // Step 3: batch_state destroyed here - next sample starts fresh
    Ok(())
}

fn process_batch_new_segments(
    batch_state: &mut BatchState,
    global_registry: &mut GlobalRegistry,
    config: &Config,
) -> Result<()> {
    // Within this batch, identical k-mer pairs share groups
    let mut batch_local_groups: HashMap<(u64, u64), u32> = HashMap::new();

    for (key, segments) in &batch_state.new_segments {
        // Assign group ID (batch-local grouping)
        let group_id = *batch_local_groups.entry(*key).or_insert_with(|| {
            let id = batch_state.next_group_id;
            batch_state.next_group_id += 1;
            id
        });

        // Register in global registry
        if !global_registry.kmer_to_group.contains_key(key) {
            global_registry.kmer_to_group.insert(*key, group_id);
            global_registry.groups.push(SegmentGroup::new(group_id, config.pack_size));
        }

        // Add segments to group
        for segment in segments {
            global_registry.groups[group_id as usize].add_segment(segment.clone());
        }
    }

    // Update global next_group_id
    global_registry.next_group_id = batch_state.next_group_id;

    Ok(())
}
```

**Tasks**:
- [ ] Implement `process_sample_batch()` function
- [ ] Implement `process_batch_new_segments()` function
- [ ] Modify main compression loop to call `process_sample_batch()` per sample
- [ ] Ensure batch state is cleared between samples

**Verification**: Single-threaded test creates archives, extract works

---

### ☐ Phase 3: Pack Writing (Write-As-You-Go)

**Implement Immediate Pack Writing**:
```rust
fn flush_pack(group: &mut SegmentGroup, config: &Config) -> Result<()> {
    // Take segments from buffer
    let segments = group.take_segments_for_pack();

    if segments.is_empty() {
        return Ok(());
    }

    // First segment in group becomes reference
    if !group.reference_written {
        let ref_seg = &segments[0];
        write_reference_segment(ref_seg)?;
        group.reference_written = true;
    }

    // Compress remaining segments
    let compressed = compress_segments(&segments[1..], config)?;

    // Write compressed pack to archive
    write_pack(group.group_id, compressed)?;

    Ok(())
}
```

**Tasks**:
- [ ] Implement `flush_pack()` for immediate writing
- [ ] Ensure reference segment handling (first in group)
- [ ] Add pack-full checks after each segment addition
- [ ] Verify ZSTD frame ordering matches C++ AGC

**Verification**: Archives match C++ AGC size within 1%

---

### ☐ Phase 4: Thread Safety

**Current**: Global state with `Arc<Mutex<...>>`

**New Approach**:
- Single-threaded batch processing (matches C++ AGC worker pattern)
- Each worker thread processes complete contigs
- Synchronization at sample boundaries

**Tasks**:
- [ ] Design thread-safe batch state management
- [ ] Implement synchronization tokens (like C++ AGC's registration stage)
- [ ] Test with `-t 1` first, then multi-threaded

**Verification**: Multi-threaded archives match single-threaded

---

### ☐ Phase 5: Testing & Verification

**Test Protocol**:
```bash
# 1. Create with RAGC batch mode
./target/release/ragc create -o ragc.agc -k 21 -s 10000 -m 20 -t 1 \
  samples/AAA_0.fa samples/AAB_0.fa samples/AAC_0.fa

# 2. Create with C++ AGC
/home/erik/agc/bin/agc create -o cpp.agc -k 21 -s 10000 -l 20 -t 1 \
  samples/AAA_0.fa samples/AAB_0.fa samples/AAC_0.fa

# 3. Compare sizes (should be identical or <1% difference)
ls -lh ragc.agc cpp.agc

# 4. Compare segment layouts
./target/release/ragc inspect ragc.agc --segment-layout > ragc_layout.csv
./target/release/ragc inspect cpp.agc --segment-layout > cpp_layout.csv
diff ragc_layout.csv cpp_layout.csv

# 5. Verify extraction
./target/release/ragc getctg ragc.agc AAB#0 "AAB#0#chrI" > ragc_aab.fa
/home/erik/agc/bin/agc getctg cpp.agc AAB_0 "AAB#0#chrI" > cpp_aab.fa
diff ragc_aab.fa cpp_aab.fa
```

**Tests**:
- [ ] 3-sample test (AAA, AAB, AAC)
- [ ] 5-sample chrV test
- [ ] 10-sample test
- [ ] Full yeast235 test

**Success Criteria**:
- ✅ Archive sizes within 1% of C++ AGC
- ✅ Segment layouts identical (same group IDs, in_group_ids)
- ✅ Extractions byte-identical
- ✅ Both implementations can read each other's archives

---

### ☐ Phase 6: Performance Optimization (Optional)

**After correctness verified**:
- [ ] Profile memory usage
- [ ] Optimize pack writing paths
- [ ] Benchmark vs C++ AGC
- [ ] Consider parallelization strategies

---

## Key Behavioral Changes

| Aspect | Old (Streaming) | New (Batch) |
|--------|-----------------|-------------|
| **State Reset** | Never | At sample boundaries |
| **K-mer Pair → Group** | Always same group | Different groups per batch |
| **Memory** | Unbounded accumulation | Bounded by batch + pack size |
| **Writing** | All at finalize | Incremental (pack-by-pack) |
| **ZSTD Frames** | Different order | Matches C++ AGC |

---

## Risk Mitigation

1. **Keep old code initially**: Add flag to switch between modes
2. **Test single-threaded first**: Verify correctness before parallelism
3. **Systematic verification**: After each phase, verify archives
4. **Incremental changes**: One phase at a time, commit frequently

---

## Expected Results

| Metric | Before | After |
|--------|--------|-------|
| Archive size vs C++ AGC | +12-16% | ~0% (byte-identical) |
| Memory usage | Unbounded | Bounded (batch + packs) |
| Segment grouping | Global | Batch-local |
| ZSTD frame order | Different | Identical |

---

## Progress Tracking

- **Started**: 2025-11-30
- **Current Phase**: Phase 1 - Add Batch State Infrastructure
- **Completed Phases**: None yet
- **Blocked**: None

---

## Notes & Observations

### 2025-11-30: Initial Plan Created
- Identified that current "streaming" is actually batch accumulation
- Designed batch-local state management matching C++ AGC
- Ready to begin Phase 1 implementation

---

## References

- `docs/CPP_AGC_BATCH_ARCHITECTURE.md` - Detailed C++ AGC analysis
- `docs/BATCH_ARCHITECTURE_PLAN.md` - Original pack-writing plan
- C++ AGC source: `/home/erik/agc/src/core/agc_compressor.cpp`
- RAGC source: `ragc-core/src/streaming_compressor_queue.rs`
