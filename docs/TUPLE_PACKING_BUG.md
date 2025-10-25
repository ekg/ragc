# Tuple Packing Bug - Root Cause Analysis

**Date**: 2025-10-24
**Status**: Root cause identified, fix pending
**Severity**: Critical - breaks C++ AGC compatibility after 7 delta-encoded samples

## Executive Summary

RAGC archives are incompatible with C++ AGC when they contain **8 or more delta-encoded samples per group**. The root cause is a mismatch in segment metadata tuple packing - RAGC doesn't implement the variable-width tuple packing that C++ AGC expects.

## Symptoms

- **Small archives (< 8 delta samples)**: ✅ C++ AGC extracts correctly
- **Large archives (≥ 8 delta samples)**: ❌ C++ AGC fails with "Corrupted archive!"
- **Reference samples**: Always work (they have `in_group_id = 0` in raw groups)
- **Delta samples 1-7**: Work (`in_group_id = 1-7`)
- **Delta samples 8+**: Fail (`in_group_id ≥ 8`)

## Discovery Process

### Test Results

**Test 1: Small chrV sequences (8KB, 2 samples)**
- CFF#2: 3,920 bp ✅
- AAA#0: 7,920 bp ✅
- **Result**: Both extract correctly with C++ AGC

**Test 2: Full genomes (12MB, 2 samples)**
- CFF#2: 10,787,791 bp ✅
- AAA#0: 12,139,594 bp ✅
- **Result**: Both extract correctly with C++ AGC

**Test 3: Full genomes (12MB, 10 samples)**
- Samples 1-7: ✅ All extract correctly
  - AAA#0, AAB#0, AAC#0, AAC#1, AAC#2, AAR#0, ABA#0
- Samples 8-10: ❌ All fail with "Corrupted archive!"
  - ABH#0, ACA#0, ACH#0

### Key Insight

The bug is **NOT about data size** - it's about the **number of delta-encoded samples**!
- 2 samples × 12MB: Works
- 10 samples × 12MB: First 7 work, last 3 fail

This pointed directly to a value threshold at 7/8.

## Root Cause: Tuple Packing Thresholds

### C++ AGC Implementation

C++ AGC uses variable-width tuple packing for segment metadata arrays (`segment.h:73-90`):

```cpp
void bytes2tuples(const vector<uint8_t>& v_bytes, vector<uint8_t>& v_tuples) {
    uint8_t me = 0;

    if(!v_bytes.empty())
        me = *max_element(v_bytes.begin(), v_bytes.end());

    if (me < 4)
        bytes2tuples_impl<4, 4>(v_bytes, v_tuples);   // Pack 4 per byte (values 0-3)
    else if (me < 6)
        bytes2tuples_impl<3, 6>(v_bytes, v_tuples);   // Pack 3 per byte (values 0-5) ← CRITICAL!
    else if (me < 16)
        bytes2tuples_impl<2, 16>(v_bytes, v_tuples);  // Pack 2 per byte (values 0-15)
    else
    {
        v_tuples = v_bytes;                            // No packing (values ≥16)
        v_tuples.emplace_back(0x10u);
    }
}
```

### The Critical Threshold

The `me < 6` check means:
- **Values 0-5**: Use 3-per-byte packing (6 possible values)
- **Values 6+**: Need 2-per-byte packing (16 possible values)

### How It Breaks

In delta-encoded groups (groups 16+):
1. First sample in group: `in_group_id = 0` (becomes reference)
2. Second sample: `in_group_id = 1` (delta)
3. Third sample: `in_group_id = 2` (delta)
4. ...
5. Seventh sample: `in_group_id = 6` (delta) - **still fits in 3-per-byte packing**
6. Eighth sample: `in_group_id = 7` (delta) - **EXCEEDS threshold! Needs 2-per-byte!**

When C++ AGC reads segment metadata:
- If all `in_group_id` values are ≤ 5: Expects 3-per-byte packing
- If any `in_group_id` value is ≥ 6: Expects 2-per-byte packing
- **If RAGC doesn't pack correctly, C++ AGC misreads the values!**

## RAGC Issue

**RAGC does not implement tuple packing for segment metadata!**

Search results:
```bash
$ grep -rn "bytes_to_tuples\|tuples_to_bytes" /home/erik/ragc --include="*.rs"
# No results
```

This means RAGC either:
1. Writes segment metadata unpacked (always 1 byte per value)
2. Uses a different packing scheme
3. Writes it in a format C++ AGC misinterprets when values exceed thresholds

## Segment Metadata Fields

The fields that get tuple-packed in C++ AGC (from `segment.cpp`):

1. **`in_group_id`** - Position within the group (0, 1, 2, 3, ..., N)
   - **This is the problematic field!**
   - Values 0-6 work (7 values), value 7+ fails
2. **`seg_part_no`** - Segment part number
3. **`is_rev_comp`** - Boolean (0 or 1)

## Testing Evidence

```bash
# Test with 10 samples
./target/release/ragc create -o test_10samples.agc -k 21 -s 10000 -m 20 \
    AAA_0.fa AAB_0.fa AAC_0.fa AAC_1.fa AAC_2.fa \
    AAR_0.fa ABA_0.fa ABH_0.fa ACA_0.fa ACH_0.fa

# Extract with C++ AGC
agc getset test_10samples.agc "AAA#0"  # ✅ Works (reference, in_group_id=0)
agc getset test_10samples.agc "AAB#0"  # ✅ Works (delta, in_group_id=1)
agc getset test_10samples.agc "AAC#0"  # ✅ Works (delta, in_group_id=2)
agc getset test_10samples.agc "AAC#1"  # ✅ Works (delta, in_group_id=3)
agc getset test_10samples.agc "AAC#2"  # ✅ Works (delta, in_group_id=4)
agc getset test_10samples.agc "AAR#0"  # ✅ Works (delta, in_group_id=5)
agc getset test_10samples.agc "ABA#0"  # ✅ Works (delta, in_group_id=6)
agc getset test_10samples.agc "ABH#0"  # ❌ FAILS (delta, in_group_id=7)
agc getset test_10samples.agc "ACA#0"  # ❌ FAILS (delta, in_group_id=8)
agc getset test_10samples.agc "ACH#0"  # ❌ FAILS (delta, in_group_id=9)
```

Exact failure at `in_group_id = 7`, which is the first value that exceeds the `me < 6` threshold!

## Investigation Update (2025-10-24)

### Tuple Packing Investigation

**Finding**: Tuple packing is NOT the root cause!

C++ AGC uses tuple packing in TWO places:
1. **Segment DATA compression** (segment.h:193-215): Only for reference segments with low repetitiveness
   - Delta packs ALWAYS use plain ZSTD (segment.h:279) with marker byte 0
   - Reference segments choose between tuple-packed (marker 1) or plain (marker 0) based on repetitiveness
2. **Collection metadata**: Uses variable-width integer encoding (`append` function), NOT tuple packing

RAGC already uses:
- Plain ZSTD for all segments (marker byte 0) - CORRECT for delta packs
- Variable-width integer encoding for collection metadata - MATCHES C++ AGC

**Conclusion**: The bug at sample #8 is NOT related to tuple packing. The root cause must be elsewhere.

## Next Investigation Steps

1. ❌ Tuple packing - ruled out
2. ⏳ Collection metadata encoding - verify exact match with C++ AGC
3. ⏳ Delta encoding format - check if in_group_id=7 triggers edge case
4. ⏳ Archive stream format - compare byte-by-byte

## Files to Investigate

- `ragc-core/src/compressor_streaming.rs` - Where segments are created
- `ragc-common/src/archive.rs` - Archive writing (if exists)
- `ragc-core/src/segment_compression.rs` - Segment packing
- C++ reference: `/home/erik/agc/src/common/segment.h:73-90`
- C++ reference: `/home/erik/agc/src/common/segment.cpp` - Implementation

## Next Steps

1. ✅ Document the bug (this file)
2. ⏳ Find where RAGC writes segment metadata
3. ⏳ Implement tuple packing matching C++ AGC
4. ⏳ Test with 10-sample archive
5. ⏳ Test with full yeast235 (235 samples)
6. ⏳ Verify all samples extract correctly with C++ AGC

## References

- Original issue: KNOWN_ISSUES.md
- C++ AGC source: `/home/erik/agc/src/common/segment.{h,cpp}`
- Test archives: `/tmp/test_10samples.agc`
- Test data: `/tmp/yeast235_individual/*.fa`
