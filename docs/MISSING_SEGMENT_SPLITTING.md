# Missing Feature: C++ AGC's Segment Splitting Algorithm

## Executive Summary

**RAGC is missing the segment splitting feature that C++ AGC uses to create additional segments from multi-sample archives.**

This explains:
- Why AAA#0 (first sample) matches perfectly: 1243 segments in both implementations
- Why AAB#0 and AAC#0 have FEWER segments in RAGC: Missing 106 segments total (47 from AAB#0, 59 from AAC#0)
- Why archives are ~10% different in size

## The Problem

### Observed Differences

**chrV test (3 samples):**
- C++ AGC: 3656 total segments (AAA#0: 1243, AAB#0: 1216, AAC#0: 1197)
- RAGC: 3550 total segments (AAA#0: 1243, AAB#0: 1169, AAC#0: 1138)
- **Missing: 106 segments** (all in non-first samples)

### Root Cause

C++ AGC has a **two-phase segment creation algorithm**:

1. **Phase 1**: `compress_contig()` splits contigs by splitter k-mers
2. **Phase 2**: `add_segment()` attempts to SPLIT segments further using "middle splitters"

RAGC only implements Phase 1.

## How C++ AGC's Segment Splitting Works

### Data Structure: `map_segments_terminators`

```cpp
// Maps each k-mer to a list of other k-mers it's been paired with in groups
unordered_map<uint64_t, vector<uint64_t>> map_segments_terminators;
```

When a NEW group is created with k-mer pair (kmer1, kmer2):
```cpp
map_segments_terminators[kmer1].push_back(kmer2);
map_segments_terminators[kmer2].push_back(kmer1);
```

### Algorithm: Finding "Middle Splitters"

In `add_segment()` (lines 1378-1520 in agc_compressor.cpp):

```cpp
// When adding a segment with k-mer pair (front, back):

// 1. Check if this k-mer pair is NOT already a known group
if (p == map_segments.end() &&
    pk.first != ~0ull && pk.second != ~0ull &&
    // 2. BUT both k-mers ARE known terminators
    map_segments_terminators.count(pk.first) &&
    map_segments_terminators.count(pk.second))
{
    // 3. Find k-mers that appear in BOTH terminator lists
    auto front_terminators = map_segments_terminators[pk.first];
    auto back_terminators = map_segments_terminators[pk.second];

    // Find intersection - these are "middle splitters"
    vector<uint64_t> shared = set_intersection(front_terminators, back_terminators);

    if (!shared.empty()) {
        uint64_t middle = shared.front();  // Take first shared k-mer

        // 4. Find where this middle k-mer occurs in the segment
        uint32_t split_pos = find_position_of_kmer(segment, middle);

        // 5. SPLIT the segment into two parts!
        segment_left = segment[0..split_pos + kmer_length];
        segment_right = segment[split_pos - kmer_length/2..end];

        // 6. Add both segments to different groups
        //    Left segment:  (front, middle) → existing group
        //    Right segment: (middle, back)  → existing group
    }
}
```

### Example

**Sample AAA#0** (first sample):
- Creates segment with k-mers (A=0x123, B=0x456)
- Creates group 30 with this pair
- Records: `map_segments_terminators[0x123] = [0x456]`
- Records: `map_segments_terminators[0x456] = [0x123]`

Later, creates another segment with k-mers (A=0x123, C=0x789):
- Creates group 31
- Updates: `map_segments_terminators[0x123] = [0x456, 0x789]`

**Sample AAB#0** (second sample):
- Encounters segment with k-mers (A=0x123, B=0x456)
- Checks: Is (0x123, 0x456) a known group? YES → group 30
- BUT ALSO: Encounters segment with k-mers (A=0x123, C=0x789)
- Checks: Is (0x123, 0x789) a known group? YES → group 31
- **NEW CASE**: Encounters segment with k-mers (A=0x123, ???)
  - The ??? k-mer is NOT in map_segments, so (0x123, ???) is not a known group
  - BUT 0x123 IS in map_segments_terminators
  - AND ??? might also be in map_segments_terminators
  - If they share a middle k-mer M (where M ∈ terminators[0x123] AND M ∈ terminators[???])
  - **SPLIT** the segment: (0x123, M) and (M, ???)

## Why AAA#0 Matches But AAB#0/AAC#0 Don't

### AAA#0 (First Sample)
- `map_segments_terminators` starts EMPTY
- No splitting can occur (no terminators exist yet)
- Only basic splitter-based segmentation happens
- **Result**: RAGC and C++ AGC create identical segments (1243 in both)

### AAB#0 (Second Sample)
- `map_segments_terminators` NOW populated from AAA#0
- C++ AGC can split segments using terminators from AAA#0
- RAGC cannot (missing this feature)
- **Result**: C++ AGC creates 1216 segments, RAGC creates 1169 (-47 segments)

### AAC#0 (Third Sample)
- `map_segments_terminators` populated from AAA#0 AND AAB#0
- Even more splitting opportunities
- **Result**: C++ AGC creates 1197 segments, RAGC creates 1138 (-59 segments)

## Implementation Requirements for RAGC

To match C++ AGC, RAGC needs:

### 1. Data Structure
```rust
// In Compressor or SegmentMap
struct SegmentMap {
    // Existing
    groups: HashMap<(u64, u64), u32>,

    // NEW: Track which k-mers appear together as terminators
    terminators: HashMap<u64, Vec<u64>>,
}
```

### 2. Update on Group Creation
```rust
fn create_new_group(&mut self, kmer1: u64, kmer2: u64, group_id: u32) {
    // Existing
    self.groups.insert((kmer1, kmer2), group_id);

    // NEW: Record terminators
    if kmer1 != u64::MAX && kmer2 != u64::MAX {
        self.terminators.entry(kmer1).or_default().push(kmer2);
        self.terminators.entry(kmer1).or_default().sort();

        if kmer1 != kmer2 {
            self.terminators.entry(kmer2).or_default().push(kmer1);
            self.terminators.entry(kmer2).or_default().sort();
        }
    }
}
```

### 3. Segment Splitting Logic in add_segment()
```rust
fn add_segment(&mut self, segment: Segment, kmer_front: u64, kmer_back: u64) -> SegmentDesc {
    let pair = (kmer_front.min(kmer_back), kmer_front.max(kmer_back));

    // Check if this pair is already a known group
    if let Some(&group_id) = self.groups.get(&pair) {
        // Known group - add segment directly
        return self.add_to_group(group_id, segment);
    }

    // NEW: Check if we can split this segment
    if kmer_front != u64::MAX && kmer_back != u64::MAX {
        if let Some(front_terms) = self.terminators.get(&kmer_front) {
            if let Some(back_terms) = self.terminators.get(&kmer_back) {
                // Find shared terminators (middle splitters)
                let shared: Vec<u64> = front_terms.iter()
                    .filter(|k| back_terms.contains(k))
                    .copied()
                    .collect();

                if let Some(&middle) = shared.first() {
                    // Find where middle k-mer occurs in segment
                    if let Some(split_pos) = find_kmer_in_segment(&segment, middle, self.kmer_length) {
                        // SPLIT the segment!
                        let (left, right) = split_segment_at(segment, split_pos, self.kmer_length);

                        // Add left part: (front, middle)
                        let left_pair = (kmer_front.min(middle), kmer_front.max(middle));
                        let left_group = self.groups[&left_pair];
                        self.add_to_group(left_group, left);

                        // Add right part: (middle, back)
                        let right_pair = (middle.min(kmer_back), middle.max(kmer_back));
                        let right_group = self.groups[&right_pair];
                        return self.add_to_group(right_group, right);
                    }
                }
            }
        }
    }

    // No split possible - create new group
    self.create_new_group(kmer_front, kmer_back, segment)
}
```

### 4. Helper Function: Find K-mer in Segment
```rust
fn find_kmer_in_segment(segment: &[u8], target_kmer: u64, k: usize) -> Option<usize> {
    // Scan through segment looking for target k-mer
    let mut kmer = CanonicalKmer::new(k);

    for (pos, &base) in segment.iter().enumerate() {
        if base > 3 {
            kmer.reset();
            continue;
        }

        kmer.push(base);
        if kmer.is_full() && kmer.value() == target_kmer {
            return Some(pos + 1 - k);  // Return start position of k-mer
        }
    }

    None
}
```

### 5. Helper Function: Split Segment
```rust
fn split_segment_at(segment: Vec<u8>, split_pos: usize, k: usize) -> (Vec<u8>, Vec<u8>) {
    // Split with overlap of k bases
    let left_end = split_pos + k;
    let right_start = split_pos - k / 2;

    let left = segment[..left_end].to_vec();
    let right = segment[right_start..].to_vec();

    (left, right)
}
```

## Verification Strategy

To verify correct implementation:

1. **Test on chrV (3 samples)**:
   - AAA#0 should produce 1243 segments (already matches)
   - AAB#0 should produce 1216 segments (currently: 1169, missing 47)
   - AAC#0 should produce 1197 segments (currently: 1138, missing 59)

2. **Check segment layout CSV**:
   ```bash
   diff cpp_chrV_layout.csv ragc_chrV_layout.csv
   ```
   Should show NO differences

3. **Archive size**:
   - C++ AGC: 572KB
   - RAGC should match within 1%

## References

- C++ AGC implementation: `/home/erik/agc/src/core/agc_compressor.cpp`
  - `add_segment()`: Lines 1288-1520 (segment splitting logic)
  - `map_segments_terminators` population: Lines 1020-1030
  - `find_cand_segment_with_missing_middle_splitter()`: Lines 1523-1700

- Analysis document: `/tmp/cpp_chrV_groups.csv` vs `/tmp/ragc_chrV_groups.csv`
