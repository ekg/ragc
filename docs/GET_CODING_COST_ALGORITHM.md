# C++ AGC's get_coding_cost Algorithm

## Purpose

Compute per-base LZDiff compression costs for a segment when encoded against a reference segment. This is used to find optimal split positions that minimize total compression cost.

## Call Chain

```
CSegment::get_coding_cost()           [segment.cpp:101-116]
    └─> CLZDiff::GetCodingCostVector() [lz_diff.cpp:159-284]
            └─> coding_cost_match()     [lz_diff.h:159-172]
            └─> coding_cost_Nrun()      [lz_diff.h:174-177]
            └─> int_len()               [lz_diff.h:179-191]
```

## Data Structures Required

### Input
- `text`: The segment sequence to analyze (vector of bases 0-3, or N=4)
- `prefix_costs`: Boolean flag controlling cost assignment
- Reference segment (stored in `lz_diff->reference`)
- LZDiff hash table index (prepared on-demand)

### Output
- `v_costs`: Vector of uint32, one per base in text
  - Contains cumulative compression cost at each position

## Algorithm: GetCodingCostVector()

### Overview
Scan through text base-by-base, computing LZDiff encoding cost by finding matches in reference.

### Pseudocode

```rust
fn get_coding_cost_vector(text: &[u8], prefix_costs: bool) -> Vec<u32> {
    let mut v_costs = Vec::with_capacity(text.len());
    let mut i = 0;
    let mut pred_pos = 0;  // Predicted position for delta encoding
    let mut no_prev_literals = 0;

    while i + key_len < text.len() {
        // 1. Compute k-mer hash at current position
        let x = get_code(&text[i..]);

        // 2. Handle N-runs (stretches of invalid bases)
        if x == INVALID {
            let nrun_len = get_nrun_len(&text[i..]);
            if nrun_len >= min_Nrun_len {
                let cost = coding_cost_Nrun(nrun_len);
                if prefix_costs {
                    v_costs.push(cost);
                    v_costs.extend(repeat(0).take(nrun_len - 1));
                } else {
                    v_costs.extend(repeat(0).take(nrun_len - 1));
                    v_costs.push(cost);
                }
                i += nrun_len;
                no_prev_literals = 0;
                continue;
            }
            // Short N-run: treat as literal
            v_costs.push(1);
            i += 1;
            pred_pos += 1;
            no_prev_literals += 1;
            continue;
        }

        // 3. Search for match in reference
        let ht_pos = murmur_hash(x) & ht_mask;
        let (found, match_pos, len_bck, len_fwd) =
            find_best_match(ht_pos, &text[i..], no_prev_literals);

        if !found {
            // No match: literal encoding
            v_costs.push(1);
            i += 1;
            pred_pos += 1;
            no_prev_literals += 1;
            continue;
        }

        // 4. Match found: compute and assign cost
        #[cfg(USE_SPARSE_HT)]
        {
            // Extend match backwards
            if len_bck > 0 {
                v_costs.truncate(v_costs.len() - len_bck);
                match_pos -= len_bck;
                pred_pos -= len_bck;
                i -= len_bck;
            }
        }

        let total_len = len_bck + len_fwd;
        let cost = coding_cost_match(match_pos, total_len, pred_pos);

        if prefix_costs {
            v_costs.push(cost);
            v_costs.extend(repeat(0).take(total_len - 1));
        } else {
            v_costs.extend(repeat(0).take(total_len - 1));
            v_costs.push(cost);
        }

        pred_pos = match_pos + total_len;
        i += total_len;
        no_prev_literals = 0;
    }

    // 5. Remaining bases (< key_len): treat as literals
    while i < text.len() {
        v_costs.push(1);
        i += 1;
    }

    v_costs
}
```

## Cost Calculation Functions

### coding_cost_match()

Computes encoding cost for a match in LZDiff format: `"dif_pos,len."`

```rust
fn coding_cost_match(ref_pos: u32, len: u32, pred_pos: u32) -> u32 {
    let dif_pos: i32 = ref_pos as i32 - pred_pos as i32;

    let mut cost = if dif_pos >= 0 {
        int_len(dif_pos as u32)
    } else {
        int_len((-dif_pos) as u32) + 1  // +1 for minus sign
    };

    cost += int_len(len - min_match_len) + 2;  // +2 for ',' and '.'
    cost
}
```

**Example:**
- Match at position 1000, length 50, pred_pos=950
- dif_pos = 1000 - 950 = 50
- len_encoded = 50 - min_match_len (e.g., 50 - 20 = 30)
- Encoding: "50,30." = 2 + 1 + 2 + 1 + 1 = 7 bytes
- Cost: int_len(50) + int_len(30) + 2 = 2 + 2 + 2 = 6 bytes

### coding_cost_Nrun()

Computes encoding cost for N-run in format: `"N_marker,len,N"`

```rust
fn coding_cost_Nrun(len: u32) -> u32 {
    2 + int_len(len - min_Nrun_len)  // 2 for start/stop markers
}
```

### int_len()

Returns number of decimal digits needed to represent a number:

```rust
fn int_len(x: u32) -> u32 {
    match x {
        0..=9         => 1,
        10..=99       => 2,
        100..=999     => 3,
        1000..=9999   => 4,
        10000..=99999 => 5,
        // ... up to 10 digits
        _             => 10,
    }
}
```

## prefix_costs Flag

Controls where cost is assigned within a match:

- **prefix_costs = true**: Assign cost to FIRST base
  ```
  Match of length 5 with cost 3:
  v_costs = [3, 0, 0, 0, 0]
  ```

- **prefix_costs = false**: Assign cost to LAST base
  ```
  Match of length 5 with cost 3:
  v_costs = [0, 0, 0, 0, 3]
  ```

## Usage in Optimal Split Position Finding

From `agc_compressor.cpp:1553-1643`:

```cpp
// Get references to existing segment groups
auto segment_id1 = map_segments[minmax(kmer_front.data(), middle)];
auto segment_id2 = map_segments[minmax(middle, kmer_back.data())];
auto seg1 = v_segments[segment_id1];
auto seg2 = v_segments[segment_id2];

// Compute costs for LEFT segment (encode against seg1 reference)
if (kmer_front.data() < middle)
    seg1->get_coding_cost(segment_dir, v_costs1, true, zstd_dctx);
else {
    seg1->get_coding_cost(segment_rc, v_costs1, false, nullptr);
    reverse(v_costs1.begin(), v_costs1.end());
}
partial_sum(v_costs1.begin(), v_costs1.end(), v_costs1.begin());

// Compute costs for RIGHT segment (encode against seg2 reference)
if (middle < kmer_back.data()) {
    seg2->get_coding_cost(segment_dir, v_costs2, false, nullptr);
    partial_sum(v_costs2.rbegin(), v_costs2.rend(), v_costs2.rbegin());
} else {
    seg2->get_coding_cost(segment_rc, v_costs2, true, nullptr);
    partial_sum(v_costs2.begin(), v_costs2.end(), v_costs2.begin());
    reverse(v_costs2.begin(), v_costs2.end());
}

// Find optimal split position
uint32_t best_sum = ~0u;
uint32_t best_pos = 0;
for (uint32_t i = 0; i < v_costs1.size(); ++i) {
    uint32_t cs = v_costs1[i] + v_costs2[i];
    if (cs < best_sum) {
        best_sum = cs;
        best_pos = i;
    }
}
```

### Key Insight

After `partial_sum`:
- `v_costs1[i]` = total cost to encode `left[0..i]` against ref1
- `v_costs2[i]` = total cost to encode `right[i..end]` against ref2
- `v_costs1[i] + v_costs2[i]` = total cost if we split at position i

Find position `i` that minimizes this sum.

## Implementation Notes for RAGC

### Requirements

1. **LZDiff index must be prepared**: Call `AssureIndex()` before cost calculation
2. **Segment must be decompressed**: Need actual sequence data, not compressed
3. **Reference sequence stored**: LZDiff needs the reference to find matches
4. **Hash table built**: For fast k-mer lookup

### Critical Details

1. **Sparse hash table** (`USE_SPARSE_HT`): Allows backwards extension of matches
   - When match found, check previous bases to extend match backwards
   - Adjust v_costs by removing last `len_bck` elements

2. **K-mer orientation**: Handle canonical k-mers
   - May need to work with reverse complement
   - Use `prefix_costs` flag + reverse to handle orientation

3. **Edge cases**:
   - Segments shorter than 2*(k+1): Cannot split
   - No reference (first segment): cost = text length (all literals)
   - Equal sequences: Empty delta in LZDiff v2

4. **Performance**:
   - Cost calculation is O(n) where n = segment length
   - Finding split position is O(n) scan
   - Total: O(n) per candidate split

## Test Strategy

To verify RAGC implementation matches C++ AGC:

1. Create identical LZDiff reference
2. Encode same test segment
3. Compare v_costs vectors element-by-element
4. Test with different prefix_costs flags
5. Test with reverse complement sequences
6. Verify split position matches C++ AGC on real data

## References

- `/home/erik/agc/src/common/segment.cpp` lines 101-116
- `/home/erik/agc/src/common/lz_diff.cpp` lines 159-284
- `/home/erik/agc/src/common/lz_diff.h` lines 159-191
- `/home/erik/agc/src/core/agc_compressor.cpp` lines 1553-1643
