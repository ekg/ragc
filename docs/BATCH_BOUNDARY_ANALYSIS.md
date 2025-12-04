# Batch Boundary Triggering Analysis

**Date**: 2025-12-04
**Status**: ROOT CAUSE IDENTIFIED
**Related to**: FIX 4 - group_id divergence from position 1

---

## Problem Statement

After implementing FIX 4 (batch-level sorting), group_id assignments now diverge from position 1:
- Position 1: RAGC group_id=131, C++ AGC group_id=411
- This indicates different batch boundaries → different group creation order

---

## C++ AGC Batch Boundary Mechanism

### Sync Token Sending (when batches are flushed)

C++ AGC has **TWO modes** for batch boundaries:

#### Mode 1: Concatenated Genomes (`concatenated_genomes == true`)

**File**: `agc_compressor.cpp:2295-2318`

```cpp
while (gio.ReadContigRaw(id, contig))
{
    if (concatenated_genomes)
    {
        // ... process contig ...

        if (++cnt_contigs_in_sample >= max_no_contigs_before_synchronization)
        {
            // Send synchronization tokens every 'pack_cardinality' contigs
            pq_contigs_desc->EmplaceManyNoCost(make_tuple(
                contig_processing_stage_t::registration, "", "", contig_t()),
                sample_priority, no_workers);

            cnt_contigs_in_sample = 0;
            --sample_priority;
        }
    }
}
```

**Behavior**: Sync tokens sent every `-l` (pack_cardinality) contigs **globally** within the concatenated file.

#### Mode 2: Separate Files (NON-concatenated, `concatenated_genomes == false`)

**File**: `agc_compressor.cpp:2346-2354`

```cpp
if (!concatenated_genomes && any_contigs_added)
{
    // Send synchronization tokens AFTER EACH SAMPLE (one file)
    pq_contigs_desc->EmplaceManyNoCost(make_tuple(
        contig_processing_stage_t::registration,
        "", "", contig_t()), sample_priority, no_workers);

    --sample_priority;
}
```

**Behavior**: Sync tokens sent **AFTER EACH SAMPLE** (all contigs from one file), **NOT** every pack_cardinality contigs.

---

## RAGC Batch Boundary Mechanism

### Sync Token Sending

RAGC sends sync tokens in **TWO places**:

#### 1. PACK_BOUNDARY: Every pack_cardinality Contigs (GLOBAL)

**File**: `agc_compressor.rs:961-1001`

```rust
// Track GLOBAL contig count and insert sync tokens every pack_size contigs
let count = self.global_contig_count.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
let need_sync = (count + 1) % self.config.pack_size == 0;

if need_sync {
    // Reached synchronization point (every pack_size contigs GLOBALLY)
    // ... decrement priority ...

    // Insert sync tokens
    for _ in 0..self.config.num_threads {
        let sync_token = ContigTask {
            // ... sync token fields ...
            is_sync_token: true,
        };
        self.queue.push(sync_token, 0)?;
    }
}
```

**Behavior**: Sync tokens sent every `-m` (pack_size) contigs **globally**, regardless of file boundaries.

#### 2. SAMPLE_BOUNDARY: When Transitioning to New Sample

**File**: `agc_compressor.rs:1010-1045`

```rust
let mut last_sample = self.last_sample_name.lock().unwrap();
if let Some(ref last) = *last_sample {
    if last != &sample_name {
        // Sample boundary detected - insert sync tokens
        for _ in 0..self.config.num_threads {
            let sync_token = ContigTask {
                // ... sync token fields ...
                is_sync_token: true,
            };
            self.queue.push(sync_token, 0)?;
        }
    }
}
```

**Behavior**: Sync tokens sent when transitioning from one sample to another.

---

## Root Cause of Divergence

### Test Configuration

**Test case**: 5 separate FASTA files (non-concatenated mode)
```bash
ragc create -o test.agc -k 21 -s 10000 -m 20 -t 1 \
    CFF#2.fa ADI#0.fa AVI_1a#0.fa CL216#0.fa AGA_1a#0.fa
```

**Each sample has ~17 contigs** (yeast chromosomes + mitochondria)

### Comparison

| Implementation | Mode | Sync Token Frequency |
|---------------|------|----------------------|
| **C++ AGC** | Non-concatenated | After EACH SAMPLE (~17 contigs) |
| **RAGC** | (Always) | Every 20 contigs (PACK_BOUNDARY) + After each sample (SAMPLE_BOUNDARY) |

### Problem

RAGC sends **EXTRA sync tokens** that C++ AGC doesn't send:

**Example with 5 samples, 17 contigs each = 85 contigs total**

**C++ AGC sync tokens**:
- After sample 1 (contig 17)
- After sample 2 (contig 34)
- After sample 3 (contig 51)
- After sample 4 (contig 68)
- After sample 5 (contig 85)
**Total: 5 sync points**

**RAGC sync tokens**:
- After contig 20 (PACK_BOUNDARY)
- After sample 1 (contig 17) (SAMPLE_BOUNDARY)
- After contig 40 (PACK_BOUNDARY)
- After sample 2 (contig 34) (SAMPLE_BOUNDARY)
- After contig 60 (PACK_BOUNDARY)
- After sample 3 (contig 51) (SAMPLE_BOUNDARY)
- After contig 80 (PACK_BOUNDARY)
- After sample 4 (contig 68) (SAMPLE_BOUNDARY)
- After sample 5 (contig 85) (SAMPLE_BOUNDARY)
**Total: 9 sync points** (4 extra PACK_BOUNDARY syncs)

### Impact on Group Assignment

Different batch boundaries → Different group creation order → Different group_ids:

- **C++ AGC**: Processes contigs 1-17, creates groups, assigns group_ids starting from X
- **RAGC**: Processes contigs 1-20, creates DIFFERENT set of groups, assigns group_ids starting from Y

**Result**: First segment at position 1 gets different group_id (131 vs 411).

---

## Solution: Conditional PACK_BOUNDARY Sync

### Option A: Disable PACK_BOUNDARY in Non-Concatenated Mode (RECOMMENDED)

Only send PACK_BOUNDARY sync tokens in concatenated genome mode (single input file).

**Implementation**:
```rust
// In add_contig_task(), only trigger PACK_BOUNDARY sync if in concatenated mode
let need_sync = (count + 1) % self.config.pack_size == 0 && self.config.concatenated_genomes;
```

**Pros**:
- Matches C++ AGC behavior exactly
- Cleaner batch boundaries (aligned with sample boundaries)
- Should achieve identical group_id assignments

**Cons**:
- Need to add `concatenated_genomes` flag to config
- Need to detect mode (1 file vs multiple files)

### Option B: Always Use SAMPLE_BOUNDARY Only

Remove PACK_BOUNDARY sync tokens entirely, always use sample boundaries.

**Pros**:
- Simpler logic
- Always matches C++ AGC in non-concatenated mode

**Cons**:
- Won't match C++ AGC in concatenated mode
- Larger batches in concatenated mode (could affect memory)

### Option C: Make PACK_BOUNDARY Match C++ AGC's `-l` Parameter

Currently RAGC uses `-m` (pack_size) for PACK_BOUNDARY, but C++ AGC uses `-l` (pack_cardinality).

**Verification needed**: Are `-m` and `-l` intended to be the same parameter?

---

## Testing Protocol

After implementing fix:

```bash
# Test with 5 separate files (non-concatenated)
./target/release/ragc create -o /tmp/ragc_fix5.agc -k 21 -s 10000 -m 20 -t 1 \
  /home/erik/scrapy/yeast235_samples/CFF#2.fa \
  /home/erik/scrapy/yeast235_samples/ADI#0.fa \
  /home/erik/scrapy/yeast235_samples/AVI_1a#0.fa \
  /home/erik/scrapy/yeast235_samples/CL216#0.fa \
  /home/erik/scrapy/yeast235_samples/AGA_1a#0.fa

/home/erik/agc/bin/agc create -o /tmp/cpp_fix5.agc -k 21 -s 10000 -l 20 -t 1 \
  /home/erik/scrapy/yeast235_samples/CFF#2.fa \
  /home/erik/scrapy/yeast235_samples/ADI#0.fa \
  /home/erik/scrapy/yeast235_samples/AVI_1a#0.fa \
  /home/erik/scrapy/yeast235_samples/CL216#0.fa \
  /home/erik/scrapy/yeast235_samples/AGA_1a#0.fa

# Compare layouts
./target/release/ragc inspect /tmp/ragc_fix5.agc --segment-layout > /tmp/ragc_layout5.csv
./target/release/ragc inspect /tmp/cpp_fix5.agc --segment-layout > /tmp/cpp_layout5.csv

# Check for divergence
python3 /tmp/compare_final.py
```

**Expected result**: NO divergence in group_id assignments from position 1.

---

## Next Steps

1. Implement Option A (conditional PACK_BOUNDARY sync)
2. Add `concatenated_genomes` detection to CLI
3. Build and test with 5-sample dataset
4. Verify group_id assignments match C++ AGC
5. Test both modes:
   - Non-concatenated (5 separate files)
   - Concatenated (single file with all samples)
