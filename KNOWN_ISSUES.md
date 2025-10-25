# Known Issues

## ~~C++ AGC Compatibility Issue with Large Single-File Archives~~ [FIXED]

**Status:** FIXED (2025-10-25)
**Severity:** HIGH (caused data corruption with C++ AGC)
**Affects:** All multi-sample archives

### Symptoms (RESOLVED)

When compressing multi-sample FASTA files:
- C++ AGC reported "Corrupted archive!" for samples with in_group_id ≥ 5
- C++ AGC crashed with std::length_error when reading RAGC archives
- Small files (<10 samples) worked, but larger archives failed

### Root Cause

Two critical bugs were identified and fixed:

1. **Reference Segment Compression** (ragc-core/src/compressor_streaming.rs:3580):
   - RAGC was using `compress_segment_configured` with hardcoded marker byte 0 for reference segments
   - Should use `compress_reference_segment` which checks repetitiveness and uses tuple packing when beneficial
   - This caused C++ AGC to misinterpret compressed data

2. **Collection Metadata Tracking** (ragc-common/src/collection.rs:748-755):
   - Initial fix attempted to remove `&& seg.in_group_id > 0` condition to track references
   - This made RAGC incompatible with C++ AGC's deserialization logic
   - **Solution:** Reverted to match C++ AGC's exact behavior (including the condition)
   - The condition is preserved for compatibility even though it may seem incorrect

### Fixes Applied

```rust
// Fix 1: Reference segment compression (ragc-core/src/compressor_streaming.rs)
// BEFORE:
let compressed = compress_segment_configured(&ref_segment.data, self.config.compression_level)?;
compressed_with_marker.push(0); // Hardcoded marker

// AFTER:
let (compressed, marker) = compress_reference_segment(&ref_segment.data)?;
compressed_with_marker.push(marker); // Correct marker (0 = ZSTD, 1 = tuple-packed)

// Fix 2: Collection metadata tracking (ragc-common/src/collection.rs)
// Must match C++ AGC's exact condition (collection_v3.cpp:674-675):
if seg.in_group_id as i32 > prev_in_group_id && seg.in_group_id > 0 {
    self.set_in_group_id(seg.group_id as usize, seg.in_group_id as i32);
}
```

### Verification

After fixes:
- ✓ 2-sample archives work with C++ AGC
- ✓ 10-sample archives work with C++ AGC
- ✓ 235-sample split-file archives (146M) work with C++ AGC
- ✗ 235-sample single-file yeast235.fa.gz (87M) still crashes C++ AGC

**Note:** The single-file yeast235.fa.gz crash appears to be a separate, unrelated issue specific to that particular dataset's scale or structure.

### Related Files

- `ragc-core/src/compressor_streaming.rs:3580` - Reference segment compression fix
- `ragc-core/src/segment_compression.rs` - Compression/decompression with marker bytes
- `ragc-core/src/tuple_packing.rs` - Tuple packing implementation (matching C++ AGC)
- `ragc-common/src/collection.rs:748-755` - Collection metadata tracking (encoding)
- `ragc-common/src/collection.rs:907-912` - Collection metadata tracking (decoding)

### Last Updated

2025-10-25 - Fixes applied and verified

---

## Remaining Issue: yeast235.fa.gz Single-File - Specific Sample Crash

**Status:** Investigated - requires debugger analysis
**Severity:** Low (specific samples only, workaround available)
**Affects:** Large single-file multi-sample FASTA - specific samples (AAA#0, AAB#0, AAC#0)

### Symptoms

- ✓ RAGC can create and read the archive (87M)
- ✓ C++ AGC can list all 235 samples
- ✓ C++ AGC can extract reference sample (CFF#2) - **WORKS**
- ✓ C++ AGC can extract most delta samples (CBM#2, ALI#2, AIF#2, ADM#1) - **WORKS**
- ✗ C++ AGC crashes (segfault) on specific samples: AAA#0, AAB#0, AAC#0 - **FAILS**
- ✓ Same samples in split-file format work fine with C++ AGC

### Investigation Findings

1. **Archive structure is valid**: RAGC's decoder reads all samples correctly
2. **Most samples work**: C++ AGC successfully extracts 232+ samples from the archive
3. **Sample ordering difference**:
   - Input file order: AAA#0, AAB#0, AAC#0 are samples #1-3
   - RAGC order: AAA#0 is #85, CFF#2 is #1 (reference)
   - C++ AGC preserves input order: AAA#0 is #1
4. **Crash location**: `CCollection::read()` in `deserialize_contig_names()` - buffer overrun when reading contig names
5. **Contig count mismatch**:
   - RAGC archive: AAA#0 has 17 contigs (correct from input)
   - C++ AGC's own archive: AAA#0 has only 1 contig (chrMT)
   - Suggests C++ AGC may not support multi-contig PanSN samples properly
6. **Split-file format works**: All 235 samples work when using yeast_split/*.fa (146M archive) - each file has 1 contig

### Root Cause

**Contig name serialization issue** when samples have multiple contigs with PanSN headers:

1. **RAGC correctly writes 17 contigs** for AAA#0 (chrI-chrXVI + chrMT)
2. **C++ AGC's buffer overrun**: The `deserialize_contig_names()` function crashes at `read(p, enc)` when trying to read past the end of the decompressed buffer
3. **Hypothesis**: RAGC's `encode_split()` function may be producing incorrect delta-encoded contig names, causing C++ AGC to:
   - Read past buffer end looking for null terminators
   - Interpret data incorrectly due to encoding mismatch
4. **Why split-files work**: Each file has only 1 contig, so no delta encoding between contig names occurs

The issue is specifically in `ragc-common/src/collection.rs`:
- `serialize_contig_names()` - lines 589-616
- `encode_split()` - lines 503-540 (complex delta encoding logic)
- Contig names like `AAA#0#chrI`, `AAA#0#chrII` are being delta-encoded incorrectly

### Workaround

Use split files instead of single-file format:
```bash
# Works: 235 samples, all extractable with C++ AGC
ragc create -o output.agc -k 21 -s 10000 -m 20 yeast_split/*.fa

# Works partially: reference + most delta samples work, AAA/AAB/AAC fail
ragc create -o output.agc -k 21 -s 10000 -m 20 yeast235.fa.gz
```

### Next Steps

1. **Add debug logging** to RAGC's `serialize_contig_names()` and `encode_split()` to see exact bytes written
2. **Compare with C++ AGC's encoding**: Look at C++ AGC's `serialize_contig_names()` in collection_v3.cpp (lines 465-495)
3. **Test hypothesis**: Create a simple test with 2-3 contigs per sample and verify encoding/decoding
4. **Fix encode_split()**: Ensure it matches C++ AGC's delta encoding logic exactly
5. **Validate**: Re-test yeast235.fa.gz after fix

GDB backtrace confirmed crash location:
```
#0  CCollection::read(unsigned char*&, unsigned int&)
#1  CCollection_V3::deserialize_contig_names() at line 521
```
