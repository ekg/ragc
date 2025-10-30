# Archive Size Difference Investigation

## Current Status

**Algorithmic Compatibility**: ✅ 100% (245/245 k-mer pairs match)
**Decompression Correctness**: ✅ 100% (SHA256 verified)
**Archive Size**: ⚠️ 2.07% larger (61,911 bytes difference)

```
RAGC:    3,051,829 bytes
C++ AGC: 2,989,918 bytes
Difference: +61,911 bytes (+2.07%)
```

## Investigation Goal

Since we have 100% algorithmic compatibility (same splitters, same segments, same groupings), the size difference must come from:
1. Metadata format differences
2. Compression settings/parameters
3. Serialization implementation differences
4. Binary format differences

## Systematic Investigation Plan

### Phase 1: Binary Format Analysis

**Goal**: Compare the binary structure of both archives byte-by-byte.

**Steps**:

1. **Header Comparison**:
   ```bash
   hexdump -C /tmp/ragc_yeast10.agc | head -50 > /tmp/ragc_header.hex
   hexdump -C /tmp/cpp_yeast10.agc | head -50 > /tmp/cpp_header.hex
   diff /tmp/ragc_header.hex /tmp/cpp_header.hex
   ```
   - Check magic numbers
   - Check version fields
   - Check header size/structure

2. **Section Size Breakdown**:
   - Create tool to parse AGC format and report size of each section:
     - Header
     - Sample metadata
     - Contig metadata
     - Segment data (compressed)
     - LZ-diff data (compressed)
     - Index structures
   - Compare section sizes between RAGC and C++ AGC

3. **Identify Largest Difference**:
   - Which section contributes most to the 61KB difference?
   - Focus investigation on that section

### Phase 2: Compression Analysis

**Goal**: Verify that compression settings match C++ AGC exactly.

**Steps**:

1. **ZSTD Compression Level**:
   - Current RAGC: Uses default (level 3?)
   - C++ AGC: Check what level it uses
   - Command: `grep -r "ZSTD_compress" /home/erik/agc/src/`

2. **Compression Context Settings**:
   - Check if C++ AGC sets any special ZSTD parameters
   - Dictionary usage?
   - Window size?
   - Thread count effects?

3. **Test Uncompressed Sizes**:
   - Extract and decompress segments from both archives
   - Compare uncompressed sizes (should be identical if algorithm matches)
   - If different → algorithmic issue remains
   - If identical → compression settings issue

### Phase 3: Metadata Comparison

**Goal**: Compare how metadata is serialized.

**Steps**:

1. **Sample Names**:
   - How are sample names stored?
   - Length prefixes? Null-terminated?
   - Padding/alignment?

2. **Contig Names**:
   - Same questions as sample names
   - Check for unnecessary padding

3. **K-mer Storage**:
   - How are k-mer pairs stored in metadata?
   - 8 bytes each (u64)?
   - Any padding between entries?

4. **Group Metadata**:
   - How is group information stored?
   - Check for redundant/extra fields in RAGC

### Phase 4: Serialization Implementation

**Goal**: Compare serialization code line-by-line.

**Steps**:

1. **RAGC Serialization**: Check `ragc-core/src/` for:
   - `write()` methods
   - Binary serialization code
   - Alignment/padding

2. **C++ AGC Serialization**: Check `/home/erik/agc/src/` for:
   - Corresponding serialization code
   - Look for format specification

3. **Compare Implementations**:
   - Are we writing extra bytes anywhere?
   - Different integer sizes (u32 vs u64)?
   - String encoding differences?

### Phase 5: Systematic Binary Diff

**Goal**: Find the exact bytes that differ.

**Steps**:

1. **Create Minimal Test Case**:
   - Single contig, single segment
   - Compare archives byte-by-byte
   - Isolate exact difference

2. **Binary Diff Tool**:
   ```bash
   cmp -l /tmp/ragc_yeast10.agc /tmp/cpp_yeast10.agc | head -20
   ```
   - Identify where files diverge
   - Map back to format structure

3. **Annotated Binary Dump**:
   - Create tool to parse AGC format with annotations
   - Show what each byte range represents
   - Makes differences obvious

## Expected Outcomes

**Best Case**: Identify simple fix (e.g., wrong compression level, extra padding)
→ Implement fix, achieve byte-for-byte identical archives

**Likely Case**: Small differences in metadata serialization
→ Adjust RAGC serialization to match C++ AGC exactly

**Worst Case**: Fundamental format difference that doesn't affect correctness
→ Document difference, accept 2% overhead if unavoidable

## Success Criteria

- [ ] Understand source of all 61,911 bytes difference
- [ ] Reduce difference to < 1% if possible
- [ ] Achieve byte-for-byte match if feasible
- [ ] Document any unavoidable differences

## Investigation Log

### 2025-10-30: Initial Measurements
- RAGC: 3,051,829 bytes
- C++ AGC: 2,989,918 bytes
- Difference: +61,911 bytes (+2.07%)
- Algorithmic compatibility: 100% ✅
- Decompression correctness: 100% ✅

### Next Steps
1. Start with Phase 1: Binary format analysis
2. Create section size breakdown tool
3. Identify which section has the most difference
