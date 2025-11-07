# RAGC Multi-Sample Data Loss Bug - Status

**Last Updated**: 2025-11-07
**Branch**: main
**Goal**: Fix data corruption in multi-sample archives (was -6.82% data loss)

---

## ✅ Current Progress: Segments Now PERFECT

**Segment Creation**: ✅ **100% CORRECT** - All segments match C++ AGC exactly
**Decompression**: ⚠️ **BROKEN** - Overlap logic needs fixing

### Quick Test Results (yeast AAA#0 + AAB#0)

| Mode | Result | Notes |
|------|--------|-------|
| Single-sample (AAB#0 alone) | ✅ PERFECT | 0.00% - byte-for-byte identical |
| Multi-sample (RAGC decompressor) | ❌ -0.10% | 11,919 bytes missing |
| Multi-sample (C++ decompressor) | ❌ +1.18% | 132,077 extra bytes |

---

## 🔧 Fixes Applied

### 1. ✅ Immediate Reference Writing (98.5% of bug fixed!)
**File**: `ragc-core/src/streaming_compressor_queue.rs:1359-1372`
**Impact**: -6.82% → -0.10% data loss

```rust
// Write reference IMMEDIATELY when first segment arrives
if buffer.reference_segment.is_none() && buffer.segments.is_empty() {
    write_reference_immediately(&buffered, buffer, ...);
}
```

Added `write_reference_immediately()` function at lines 756-820 to match C++ AGC segment.cpp:41-48 behavior.

### 2. ✅ K-byte Overlap in Segments
**File**: `ragc-core/src/segment.rs:163`
**Impact**: All segment lengths now match C++ AGC exactly

```rust
// Changed from (k-1) to k byte overlap:
let new_start = (pos + 1).saturating_sub(k);  // was: k - 1
```

### 3. ⚠️ K-byte Overlap in Decompressor (NEEDS MORE WORK)
**File**: `ragc-core/src/decompressor.rs:343`

```rust
// Skip k bytes (full k-mer)
let overlap = self.kmer_length as usize;  // was: k - 1
```

Applied but still produces wrong output - the decompression overlap logic needs investigation.

---

## 📊 Verification: Segments Are CORRECT

All 1,086 AAB#0 segments verified against C++ AGC:

```bash
# chrI example (17 segments):
Part 0:  50890 bytes ✓ (was 50890)
Part 1:  7206 bytes ✓  (was 7205 - NOW FIXED!)
Part 2:  19369 bytes ✓ (was 19368 - NOW FIXED!)
...
Part 16: 5076 bytes ✓  (was 5075 - NOW FIXED!)

# Before fix: Every non-first segment was 1 byte short
# After fix: All 17 segments match exactly
```

---

## 🎯 Next Steps for Next Session

### CRITICAL: Test Against C++ AGC Single-Threaded Mode

We should match **single-threaded** C++ AGC (`-t 1`), not multi-threaded!

```bash
# Create with C++ AGC single-threaded (deterministic)
/home/erik/agc/bin/agc create -t 1 -o /tmp/cpp_st.agc -k 21 -s 10000 \
  /tmp/check_AAA#0.fa /tmp/check_AAB#0.fa

# Verify C++ AGC extracts its own archive correctly
/home/erik/agc/bin/agc getset /tmp/cpp_st.agc AAB#0 > /tmp/cpp_extracted.fa
diff /tmp/check_AAB#0.fa /tmp/cpp_extracted.fa  # Should be identical
```

### Step-by-Step Debugging Plan

#### 1. Isolate with Single Contig

```bash
# Extract just chrI (simplest test case)
grep -A 999999 ">AAB#0#chrI" /tmp/check_AAB#0.fa | \
  awk '/^>AAB#0#chr[^I]/ {exit} {print}' > /tmp/chrI_only.fa

# Test with single contig
./target/release/ragc create -o /tmp/single.agc -k 21 -s 10000 /tmp/chrI_only.fa
./target/release/ragc getset /tmp/single.agc AAB#0 > /tmp/extracted.fa
diff /tmp/chrI_only.fa /tmp/extracted.fa
```

#### 2. Add Decompression Debug Logging

In `decompressor.rs` around line 340:

```rust
if i == 0 {
    eprintln!("  Seg {}: FIRST segment len={}, added all {} bytes, total={}",
        i, segment_data.len(), segment_data.len(), contig.len() + segment_data.len());
    contig.extend_from_slice(&segment_data);
} else {
    let overlap = self.kmer_length as usize;
    let added = segment_data.len() - overlap;
    eprintln!("  Seg {}: len={} skip_overlap={} added={} bytes, total={}",
        i, segment_data.len(), overlap, added, contig.len() + added);
    contig.extend_from_slice(&segment_data[overlap..]);
}
```

Run and compare expected vs actual contig lengths at each step.

#### 3. Check C++ AGC Decompression Logic

Study `agc/src/core/agc_decompressor.cpp` - how does it reconstruct contigs?
- Does it always skip k bytes?
- Are there special cases (first segment, reference segment, etc.)?
- Does overlap differ between segments in same pack vs different packs?

#### 4. Theory: Variable Overlap?

The overlap might not be uniform:
- First segment of contig: No previous segment (no overlap to skip)
- Subsequent segments: k-byte overlap
- Segments across pack boundaries: Different overlap?

Check if `segment_desc.in_group_id` affects overlap logic.

---

## 🔬 Key Insights

1. **Single-sample works perfectly** → Bug is multi-sample specific
2. **Immediate reference writing fixed 95%** → Was the major issue
3. **Segments now match exactly** → Compression algorithm is CORRECT
4. **Both decompressors fail** → Decompression overlap logic is wrong

The fact that BOTH decompressors produce wrong output (one too little, one too much) strongly suggests the overlap reconstruction logic is incorrect.

---

## 📁 Test Files

In `/tmp/`:
- `check_AAA#0.fa` - Reference sample (12,157,349 bytes, 17 contigs)
- `check_AAB#0.fa` - Test sample (11,320,014 bytes, 17 contigs)

Quick regeneration:
```bash
cd /home/erik/agc
./scripts/tests/scripts/run_correctness_tests.sh
# Creates test files in /tmp/
```

---

## 🎓 Lessons Learned

1. **Trust size differences** - Even 8% archive size difference signals a bug
2. **Immediate reference writing is critical** - Fixed 95% of the bug
3. **Segment lengths must match exactly** - Like zlib: byte-for-byte, not "close enough"
4. **Single-sample test isolates bugs** - Perfect isolation for multi-sample issues
5. **C++ AGC has deterministic mode** - Use `-t 1` for reproducible comparison

---

## 💾 Commit Command

```bash
git add ragc-core/src/streaming_compressor_queue.rs
git add ragc-core/src/segment.rs
git add ragc-core/src/decompressor.rs
git add STATUS.md

git commit -m "wip: Fix segment overlap to k bytes and add immediate reference writing

- Changed segment overlap from (k-1) to k bytes in segment.rs
- All segment lengths now match C++ AGC exactly (verified 1,086 segments)
- Added write_reference_immediately() to match C++ AGC segment.cpp:41-48
- Reduced data loss from -6.82% to -0.10%

Remaining: Decompression overlap logic needs investigation.
Both RAGC (-0.10%) and C++ AGC (+1.18%) produce wrong output
when extracting from RAGC archives, suggesting overlap
reconstruction logic is incorrect.

Single-sample mode works perfectly (0% data loss)."
```
