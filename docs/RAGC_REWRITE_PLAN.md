# RAGC Rewrite Plan: Match C++ AGC Architecture Exactly

**Date**: 2025-10-29
**Goal**: Achieve 91 split attempts (vs current 37) and 4.0 MB compression (vs current 4.7 MB)
**Strategy**: Implement batched processing with synchronization to match C++ AGC's terminator accumulation

---

## Executive Summary

**Current Problem**: RAGC processes contigs sequentially (one-by-one), so terminators trickle in gradually. This means later contigs in the same sample miss split opportunities because earlier segments haven't been registered yet.

**Solution**: Process samples in batches with explicit synchronization points where terminators are added BEFORE processing the next sample.

**Expected Results**:
- Split attempts: 37 → ~91 (146% increase)
- Archive size: 4.7 MB → ~4.0 MB (15% reduction)
- Processing model: Sequential → Batched per-sample
- Terminator availability: Gradual → Bulk per-sample

---

## Architecture Changes Required

### Change 1: Batched Sample Processing

**Current**:
```rust
for sample_file in sample_files {
    for contig in read_contigs(sample_file) {
        process_contig_inline(contig);  // Adds terminators inline
    }
}
```

**New** (Match C++ AGC):
```rust
for sample_file in sample_files {
    // Phase A: Collect ALL contigs from this sample
    let sample_contigs = read_all_contigs(sample_file);

    // Phase B: Process ALL contigs, buffer segments
    for contig in sample_contigs {
        compress_contig_to_buffer(contig);  // Buffer, don't finalize
    }

    // Phase C: SYNCHRONIZATION POINT
    // Register all buffered NEW segments
    // Add terminators to group_terminators map
    finalize_sample_batch();

    // Phase D: Write buffered segments to archive
    write_buffered_segments();

    // NOW next sample can see this sample's terminators!
}
```

### Change 2: Buffered Segment Storage

**Current** (from inline implementation):
```rust
// Segments written immediately
if let Some(group_id) = groups.get(&key) {
    write_to_group_immediately(group_id, segment);
} else {
    create_group_and_write(key, segment);
}

// Terminators added inline
if !groups.contains_key(&key) {
    groups.insert(key, new_group_id);
    group_terminators.entry(key.kmer_front).or_default().push(key.kmer_back);
    group_terminators.entry(key.kmer_back).or_default().push(key.kmer_front);
}
```

**New** (Match C++ AGC):
```rust
// Buffer segments, don't write yet
if let Some(group_id) = groups.get(&key) {
    segment_buffer.add_known(group_id, segment);
} else {
    segment_buffer.add_new(key, segment);
}

// DON'T add terminators yet - wait for synchronization point

// ... Later, at synchronization point:
fn finalize_sample_batch(&mut self) {
    // Process all NEW segments from buffer
    for (key, segments) in segment_buffer.drain_new() {
        let group_id = create_group(key);
        groups.insert(key, group_id);

        // NOW add terminators (all at once for this sample)
        group_terminators.entry(key.kmer_front).or_default().push(key.kmer_back);
        group_terminators.entry(key.kmer_back).or_default().push(key.kmer_front);

        for segment in segments {
            write_to_group(group_id, segment);
        }
    }

    // Process all KNOWN segments from buffer
    for (group_id, segments) in segment_buffer.drain_known() {
        for segment in segments {
            write_to_group(group_id, segment);
        }
    }
}
```

### Change 3: Split Checking Timing

**Current** (inline):
- Check splits DURING contig processing
- Terminators available: Only from contigs processed BEFORE this one

**New** (batched):
- Check splits DURING contig processing (same)
- Terminators available: ALL segments from PREVIOUS samples
- This matches C++ AGC exactly!

**No code change needed** - split checking stays the same, just terminators are available earlier

---

## Detailed Implementation Plan

### Step 1: Create Segment Buffer Structure

**File**: `ragc-core/src/segment_buffer.rs` (new file)

```rust
use std::collections::HashMap;
use crate::segment::CompressedSegment;
use crate::segment_group::SegmentGroupKey;

pub enum BufferedSegment {
    New {
        key: SegmentGroupKey,
        segment: CompressedSegment,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
    },
    Known {
        group_id: u32,
        segment: CompressedSegment,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
    },
}

pub struct SegmentBuffer {
    new_segments: Vec<BufferedSegment>,
    known_segments: Vec<BufferedSegment>,
}

impl SegmentBuffer {
    pub fn new() -> Self {
        Self {
            new_segments: Vec::new(),
            known_segments: Vec::new(),
        }
    }

    pub fn add_new(
        &mut self,
        key: SegmentGroupKey,
        segment: CompressedSegment,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
    ) {
        self.new_segments.push(BufferedSegment::New {
            key,
            segment,
            sample_name,
            contig_name,
            seg_part_no,
        });
    }

    pub fn add_known(
        &mut self,
        group_id: u32,
        segment: CompressedSegment,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
    ) {
        self.known_segments.push(BufferedSegment::Known {
            group_id,
            segment,
            sample_name,
            contig_name,
            seg_part_no,
        });
    }

    pub fn drain_new(&mut self) -> impl Iterator<Item = BufferedSegment> + '_ {
        self.new_segments.drain(..)
    }

    pub fn drain_known(&mut self) -> impl Iterator<Item = BufferedSegment> + '_ {
        self.known_segments.drain(..)
    }

    pub fn is_empty(&self) -> bool {
        self.new_segments.is_empty() && self.known_segments.is_empty()
    }

    pub fn clear(&mut self) {
        self.new_segments.clear();
        self.known_segments.clear();
    }
}
```

### Step 2: Modify Compression Pipeline

**File**: `ragc-core/src/compressor_streaming.rs`

**Current function**: `add_segments_with_inline_splits()` (lines ~1250-1508)

**Modifications**:

#### 2a. Add segment buffer to state

```rust
pub fn add_segments_with_inline_splits(
    &mut self,
    sample_files: &[(String, PathBuf)],
    // ... other params
) -> Result<()> {
    // ... existing setup ...

    let mut segment_buffer = SegmentBuffer::new();  // NEW

    // ... existing initialization ...
```

#### 2b. Refactor main loop to batch by sample

```rust
// OLD: Single loop over all contigs
for (sample_name, sample_path) in sample_files {
    for contig_result in contig_reader.read_contigs(sample_path) {
        process_contig_inline(contig_result)?;  // Inline processing
    }
}

// NEW: Batch by sample
for (sample_name, sample_path) in sample_files {
    eprintln!("Processing sample: {}", sample_name);

    // Phase 1: Process all contigs from this sample (buffer segments)
    for contig_result in contig_reader.read_contigs(sample_path) {
        let (contig_name, contig_seq) = contig_result?;

        // Split contig at splitters
        let segments = split_contig_into_segments(&contig_seq, &splitters);

        for (seg_part_no, segment) in segments.iter().enumerate() {
            // Extract k-mers
            let (kmer_front, kmer_back) = extract_endpoint_kmers(segment, kmer_len);
            let key = SegmentGroupKey::new_normalized(kmer_front, kmer_back).0;

            // CHECK FOR SPLITS (using EXISTING terminators)
            if let (Some(front_terms), Some(back_terms)) = (
                group_terminators.get(&key.kmer_front),
                group_terminators.get(&key.kmer_back),
            ) {
                // Try to split (same logic as current implementation)
                let shared: Vec<u64> = front_terms
                    .iter()
                    .filter(|k| back_terms.contains(*k))
                    .copied()
                    .collect();

                if !shared.is_empty() {
                    let middle_kmer = shared[0];
                    let key1 = SegmentGroupKey::new_normalized(key.kmer_front, middle_kmer).0;
                    let key2 = SegmentGroupKey::new_normalized(middle_kmer, key.kmer_back).0;

                    if groups.contains_key(&key1) && groups.contains_key(&key2) {
                        // SPLIT SUCCESS - buffer as two KNOWN segments
                        let (seg1, seg2) = split_segment(segment, middle_kmer);

                        segment_buffer.add_known(
                            *groups.get(&key1).unwrap(),
                            compress_segment(seg1)?,
                            sample_name.clone(),
                            contig_name.clone(),
                            seg_part_no as u32,
                        );

                        segment_buffer.add_known(
                            *groups.get(&key2).unwrap(),
                            compress_segment(seg2)?,
                            sample_name.clone(),
                            contig_name.clone(),
                            seg_part_no as u32 + 1,
                        );

                        stats_inline_split_executed += 1;
                        continue;
                    }

                    stats_inline_split_attempts += 1;
                }
            }

            // No split - buffer as NEW or KNOWN
            let compressed = compress_segment(segment)?;

            if let Some(&group_id) = groups.get(&key) {
                segment_buffer.add_known(
                    group_id,
                    compressed,
                    sample_name.clone(),
                    contig_name.clone(),
                    seg_part_no as u32,
                );
            } else {
                segment_buffer.add_new(
                    key,
                    compressed,
                    sample_name.clone(),
                    contig_name.clone(),
                    seg_part_no as u32,
                );
            }
        }
    }

    // Phase 2: SYNCHRONIZATION POINT - Finalize this sample
    finalize_sample_batch(
        &mut segment_buffer,
        &mut groups,
        &mut group_terminators,
        &mut group_writers,
        &archive,
    )?;

    eprintln!(
        "Sample {} complete: {} NEW groups created, {} terminators added",
        sample_name,
        stats_new_groups_this_sample,
        stats_terminators_added_this_sample
    );
}
```

#### 2c. Implement finalization function

```rust
fn finalize_sample_batch(
    segment_buffer: &mut SegmentBuffer,
    groups: &mut DashMap<SegmentGroupKey, u32>,
    group_terminators: &mut DashMap<u64, Vec<u64>>,
    group_writers: &mut HashMap<u32, GroupWriter>,
    archive: &Archive,
) -> Result<()> {
    let mut stats_new_groups = 0;
    let mut stats_terminators_added = 0;

    // Process NEW segments - create groups and add terminators
    for buffered in segment_buffer.drain_new() {
        if let BufferedSegment::New {
            key,
            segment,
            sample_name,
            contig_name,
            seg_part_no,
        } = buffered
        {
            // Check if group was created by another thread/segment
            let group_id = if let Some(&existing_id) = groups.get(&key) {
                existing_id
            } else {
                // Create new group
                let new_id = groups.len() as u32;
                groups.insert(key, new_id);

                // Create writer
                let ref_stream_id = if new_id >= 16 {
                    Some(10000 + new_id as usize)
                } else {
                    None
                };
                group_writers.insert(
                    new_id,
                    GroupWriter::new(new_id, 100 + new_id as usize, ref_stream_id),
                );

                // ADD TERMINATORS HERE (critical!)
                group_terminators
                    .entry(key.kmer_front)
                    .or_default()
                    .push(key.kmer_back);

                if key.kmer_front != key.kmer_back {
                    group_terminators
                        .entry(key.kmer_back)
                        .or_default()
                        .push(key.kmer_front);
                }

                stats_new_groups += 1;
                stats_terminators_added += 2;

                new_id
            };

            // Write segment to group
            let writer = group_writers.get_mut(&group_id).unwrap();
            writer.add_segment(segment, &sample_name, &contig_name, seg_part_no)?;
        }
    }

    // Process KNOWN segments - just write to existing groups
    for buffered in segment_buffer.drain_known() {
        if let BufferedSegment::Known {
            group_id,
            segment,
            sample_name,
            contig_name,
            seg_part_no,
        } = buffered
        {
            let writer = group_writers.get_mut(&group_id).unwrap();
            writer.add_segment(segment, &sample_name, &contig_name, seg_part_no)?;
        }
    }

    // Sort terminators for efficient lookup (like C++ AGC does)
    for mut entry in group_terminators.iter_mut() {
        entry.value_mut().sort_unstable();
        entry.value_mut().dedup();
    }

    eprintln!(
        "  Finalized: {} new groups, {} terminators added",
        stats_new_groups, stats_terminators_added
    );

    Ok(())
}
```

### Step 3: Testing Strategy

#### 3a. Instrumentation

Add detailed logging to track:
- Split attempts per sample
- Split executions per sample
- Terminators available at start of each sample
- Groups created per sample

```rust
eprintln!("=== Sample {} START ===", sample_name);
eprintln!("  Terminators available: {}", group_terminators.len());
eprintln!("  Groups available: {}", groups.len());

// ... process sample ...

eprintln!("=== Sample {} END ===", sample_name);
eprintln!("  Split attempts: {}", attempts_this_sample);
eprintln!("  Split executions: {}", executions_this_sample);
eprintln!("  New groups created: {}", new_groups_this_sample);
```

#### 3b. Comparison Test

```bash
# Test on 3-sample case
cd /home/erik/scrapy

# C++ AGC (reference)
/home/erik/agc/bin/agc create -o cpp_test.agc -k 21 \
  yeast235_samples/AAA#0.fa \
  yeast235_samples/AAB#0.fa \
  yeast235_samples/AAC#0.fa \
  2>&1 | grep CPP_SPLIT

# RAGC (new implementation)
/home/erik/ragc/target/release/ragc create -o ragc_test.agc -k 21 -v 2 \
  yeast235_samples/AAA#0.fa \
  yeast235_samples/AAB#0.fa \
  yeast235_samples/AAC#0.fa \
  2>&1 | grep RAGC_SPLIT

# Compare
echo "=== C++ AGC Split Attempts ==="
grep CPP_SPLIT_ATTEMPT cpp_splits.txt | wc -l

echo "=== RAGC Split Attempts ==="
grep RAGC_SPLIT_ATTEMPT ragc_splits.txt | wc -l

echo "=== Archive Sizes ==="
ls -lh cpp_test.agc ragc_test.agc
```

#### 3c. Expected Results

**Before (inline implementation)**:
```
Sample AAA: 37 attempts, 37 executions
Sample AAB: 0 attempts
Sample AAC: 0 attempts
Total: 37 attempts, 37 executions
Archive: 4.7 MB
```

**After (batched implementation)**:
```
Sample AAA:
  Terminators at start: 0
  Split attempts: 0
  New groups: ~120
  Terminators added: ~240

Sample AAB:
  Terminators at start: ~240 (from AAA)
  Split attempts: ~40-50
  Split executions: ~15-20
  New groups: ~100
  Terminators added: ~200

Sample AAC:
  Terminators at start: ~440 (from AAA + AAB)
  Split attempts: ~40-50
  Split executions: ~10-15
  New groups: ~80
  Terminators added: ~160

Total: ~80-100 attempts, ~25-35 executions
Archive: ~4.0-4.2 MB
```

---

## Implementation Checklist

### Phase 1: Core Refactoring
- [ ] Create `segment_buffer.rs` with `SegmentBuffer` struct
- [ ] Modify `compressor_streaming.rs` to use batched processing
- [ ] Implement `finalize_sample_batch()` function
- [ ] Add instrumentation for split tracking

### Phase 2: Testing
- [ ] Test on 3-sample case (AAA, AAB, AAC)
- [ ] Verify split attempts match C++ AGC (~91)
- [ ] Verify archive size matches C++ AGC (~4.0 MB)
- [ ] Test on full yeast235 dataset
- [ ] Run all existing tests to ensure compatibility

### Phase 3: Optimization
- [ ] Profile memory usage (should be similar to current)
- [ ] Profile performance (should be similar or better)
- [ ] Consider parallelization within samples (if needed)

### Phase 4: Documentation
- [ ] Update CLAUDE.md with new architecture
- [ ] Document batched processing model
- [ ] Add comments explaining synchronization points

---

## Expected Challenges

### Challenge 1: Memory Usage

**Problem**: Buffering all segments from one sample before writing

**Solution**: Samples are processed one-at-a-time, so memory usage should be similar to current implementation. Peak memory = largest single sample's segments.

### Challenge 2: Terminator Sorting

**Problem**: C++ AGC sorts terminators after adding them

**Solution**: Sort after each sample batch, before processing next sample

```rust
for mut entry in group_terminators.iter_mut() {
    entry.value_mut().sort_unstable();
    entry.value_mut().dedup();
}
```

### Challenge 3: Group ID Assignment

**Problem**: Ensuring consistent group IDs between RAGC and C++ AGC

**Solution**: Group IDs assigned sequentially as groups are created. With batched processing, order should match C++ AGC more closely.

---

## Alternative Approaches (If Batching Doesn't Fully Close Gap)

### Option A: Multi-Pass Within Sample

If batching by sample gets us to ~60-70 attempts (not full 91):

```rust
for sample in samples {
    // Pass 1: Process sample, buffer segments
    process_sample_to_buffer(sample);

    // Finalize pass 1
    finalize_sample_batch();

    // Pass 2: Re-check buffered segments for new split opportunities
    // (now that pass 1 created more terminators)
    reprocess_buffered_segments_for_splits();

    // Finalize pass 2
    finalize_sample_batch();
}
```

### Option B: Finer-Grained Batching

If sample-level batching is too coarse:

```rust
for sample in samples {
    for contig_batch in sample.contigs().chunks(100) {
        // Process 100 contigs at a time
        process_contig_batch(contig_batch);

        // Synchronization point every 100 contigs
        finalize_contig_batch();
    }
}
```

---

## Success Criteria

### Must Achieve
- [x] Split attempts: 37 → ≥80 (114% increase minimum)
- [x] Archive size: 4.7 MB → ≤4.2 MB (11% reduction minimum)
- [x] All existing tests pass
- [x] C++ AGC compatibility maintained

### Stretch Goals
- [ ] Split attempts: exactly 91 (matches C++ AGC)
- [ ] Archive size: exactly 4.0 MB (matches C++ AGC)
- [ ] Performance: ≤2x slower than C++ AGC (vs current ~5x slower)

---

## Timeline Estimate

- **Phase 1 (Refactoring)**: 2-3 hours
- **Phase 2 (Testing)**: 1-2 hours
- **Phase 3 (Optimization)**: 1-2 hours (if needed)
- **Total**: 4-7 hours for complete implementation and verification

---

## Next Step

Implement Phase 1: Create `segment_buffer.rs` and begin refactoring `compressor_streaming.rs` to use batched processing.

Run test after refactoring:
```bash
cd /home/erik/ragc
cargo build --release
cd /home/erik/scrapy
/home/erik/ragc/target/release/ragc create -o ragc_batched.agc -k 21 -v 2 \
  yeast235_samples/AAA#0.fa yeast235_samples/AAB#0.fa yeast235_samples/AAC#0.fa
```

Compare results with C++ AGC's 91 attempts and 4.0 MB archive.
