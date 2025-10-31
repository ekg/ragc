# Segmentation Algorithm Comparison: C++ AGC vs RAGC

**Date**: 2025-10-30
**Purpose**: Prove whether RAGC's segmentation matches C++ AGC exactly

---

## The Question

User asserts: **"If we find the same splitters, we MUST create the same segments."**

Current observation:
- Both find exactly **11,771 splitters**
- C++ AGC reports: **12,182 segments**
- RAGC reports: **22,470 total segments**, **11,001 unique groups**

**Hypothesis**: The numbers are measuring different things, OR there's a bug in one implementation.

---

## C++ AGC: `compress_contig` (agc_compressor.cpp:2000-2054)

```cpp
bool compress_contig(..., contig_t& contig, ...) {
    CKmer kmer(kmer_length, kmer_mode_t::canonical);
    uint64_t pos = 0;
    uint64_t split_pos = 0;
    CKmer split_kmer(kmer_length, kmer_mode_t::canonical);
    uint32_t seg_part_no = 0;

    // Main loop: scan through contig
    for (auto x : contig) {
        if (x >> 2)  // x > 3 (non-ACGT)
            kmer.Reset();
        else {
            kmer.insert_canonical(x);

            if (kmer.is_full()) {
                uint64_t d = kmer.data_canonical();

                // Split at EVERY splitter occurrence (no distance check!)
                if (bloom_splitters.check(d) && hs_splitters.check(d)) {
                    // Create segment from split_pos to pos+1
                    add_segment(sample_name, id, seg_part_no,
                        move(get_part(contig, split_pos, pos + 1 - split_pos)),
                        split_kmer, kmer, ...);

                    ++seg_part_no;

                    // Create k-base overlap for next segment
                    split_pos = pos + 1 - kmer_length;
                    split_kmer = kmer;
                    kmer.Reset();  // CRITICAL: Reset k-mer after split
                }
            }
        }
        ++pos;
    }

    // Add final segment
    if (split_pos < contig.size())
        add_segment(..., get_part(contig, split_pos, contig.size() - split_pos), ...);

    return true;
}
```

**Key points**:
1. Splits at **every** splitter k-mer (no distance check)
2. Creates k-base overlap: `split_pos = pos + 1 - kmer_length`
3. Resets k-mer after each split
4. Adds final segment from last split to contig end

**No complex end-of-contig logic!**

---

## RAGC: `split_at_splitters_with_size` (segment.rs:70-201)

```rust
pub fn split_at_splitters_with_size(
    contig: &Contig,
    splitters: &HashSet<u64>,
    k: usize,
    _min_segment_size: usize,  // UNUSED
) -> Vec<Segment> {
    let mut segments = Vec::new();
    let mut kmer = Kmer::new(k as u32, KmerMode::Canonical);
    let mut segment_start = 0;
    let mut front_kmer = MISSING_KMER;
    let mut recent_kmers: Vec<(usize, u64)> = Vec::new();

    // Main loop
    for (pos, &base) in contig.iter().enumerate() {
        if base > 3 {
            kmer.reset();
            recent_kmers.clear();
        } else {
            kmer.insert(base as u64);

            if kmer.is_full() {
                let kmer_value = kmer.data();
                recent_kmers.push((pos, kmer_value));

                // Split at EVERY splitter occurrence (no distance check)
                if splitters.contains(&kmer_value) {
                    let segment_end = pos + 1;
                    let segment_data = contig[segment_start..segment_end].to_vec();

                    if !segment_data.is_empty() {
                        segments.push(Segment::new(segment_data, front_kmer, kmer_value));
                    }

                    // Create k-base overlap
                    segment_start = (pos + 1).saturating_sub(k);
                    front_kmer = kmer_value;
                    recent_kmers.clear();
                    kmer.reset();  // Reset k-mer after split
                }
            }
        }
    }

    // END-OF-CONTIG HANDLING (NOT IN C++ AGC!)
    // Lines 129-166: Complex backward scan
    for (pos, kmer_value) in recent_kmers.iter().rev() {
        if splitters.contains(kmer_value) {
            let segment_end = pos + 1;
            let remaining_after = contig.len() - segment_end;

            if remaining_after > k {
                // Create additional split
                ...
                break;
            }
        }
    }

    // Add final segment
    if segment_start < contig.len() {
        let segment_data = contig[segment_start..].to_vec();
        if !segment_data.is_empty() {
            segments.push(Segment::new(segment_data, front_kmer, MISSING_KMER));
        }
    }

    // FINAL SEGMENT MERGING (NOT IN C++ AGC!)
    // Lines 186-198
    if segments.len() >= 2 {
        let last_idx = segments.len() - 1;
        if segments[last_idx].data.len() < k {
            // Merge last two segments
            let last_seg = segments.pop().unwrap();
            let second_last = segments.last_mut().unwrap();
            second_last.data.extend_from_slice(&last_seg.data);
            second_last.back_kmer = last_seg.back_kmer;
        }
    }

    segments
}
```

**Key differences**:
1. ✅ Main loop matches C++ AGC
2. ❌ **Extra end-of-contig backward scan** (lines 129-166)
3. ❌ **Extra final segment merging** (lines 186-198)
4. ❌ Stores `recent_kmers` (not in C++ AGC)

---

## Critical Difference: End-of-Contig Handling

### C++ AGC Behavior

```cpp
// After main loop:
if (split_pos < contig.size())
    add_segment(..., get_part(contig, split_pos, contig.size() - split_pos), ...);
```

**Simple**: Just adds whatever remains from last split to contig end.

### RAGC Behavior

```rust
// After main loop: Look backward for additional split opportunities!
for (pos, kmer_value) in recent_kmers.iter().rev() {
    if splitters.contains(kmer_value) {
        let segment_end = pos + 1;
        let remaining_after = contig.len() - segment_end;

        if remaining_after > k {
            // Create EXTRA split that C++ AGC doesn't make!
            if segment_end > segment_start {
                let segment_data = contig[segment_start..segment_end].to_vec();
                segments.push(Segment::new(segment_data, front_kmer, *kmer_value));
                segment_start = (pos + 1).saturating_sub(k);
                front_kmer = *kmer_value;
            }
            break;
        }
    }
}
```

**Complex**: Scans backward through recent k-mers looking for additional split opportunities.

**This is NOT in C++ AGC!**

---

## Bug Hypothesis

RAGC's extra end-of-contig logic (lines 129-166) is creating **different segment boundaries** than C++ AGC.

Specifically:
- C++ AGC: Stops splitting when last splitter k-mer is found in main loop
- RAGC: Continues looking backward at contig end for one more split

This could result in:
- Different number of segments per contig
- Different segment boundaries
- Different archive structure

---

## Test Plan

To prove this hypothesis, create a test that:

1. Takes a single contig (e.g., AAA#0#chrI)
2. Applies RAGC's current algorithm
3. Applies C++ AGC's exact algorithm (without backward scan)
4. Compares segment counts and boundaries

If counts differ → **RAGC has a bug** (extra logic causing differences)
If counts match → Numbers are measuring different things

---

## Expected Outcome

If RAGC's backward scan creates extra splits:
- RAGC would create MORE segments than C++ AGC
- But RAGC reports FEWER unique groups (11,001 vs 12,182)
- This suggests the bug might be in the OPPOSITE direction
- Or the numbers are indeed measuring different things

**Next step**: Implement exact C++ AGC algorithm and compare.
