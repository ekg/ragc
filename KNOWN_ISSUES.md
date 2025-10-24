# Known Issues

## C++ AGC Compatibility Issue with Large Single-File Archives

**Status:** Under investigation  
**Severity:** Medium (workaround available)  
**Affects:** Single-file multi-sample FASTA compression with >100 samples

### Symptoms

When compressing a large single-file multi-sample FASTA (e.g., 235 samples, ~3GB):
- C++ AGC can extract **reference samples** correctly
- C++ AGC **cannot extract delta-encoded samples** (returns empty output)
- RAGC's own decompressor works correctly for all samples
- Small files (<10 samples) work fine with both tools

### Example

```bash
# Create archive from large single PanSN file
./target/release/ragc create -o large.agc -k 21 -s 10000 -m 20 yeast235.fa

# C++ AGC extraction results:
agc getset large.agc "CFF#2"  # ✓ Works (reference sample)
agc getset large.agc "AAA#0"  # ✗ Fails (delta sample, returns empty)

# RAGC extraction works for all:
./target/release/ragc getset large.agc "AAA#0"  # ✓ Works
```

### Root Cause

The issue appears to be size-dependent and related to delta encoding in `add_contigs_with_splitters`:
- Groups 0-15 use raw encoding (reference samples work)
- Groups 16+ use LZ delta encoding (these fail with C++ AGC)
- The bug only manifests when total data size is large (>1GB estimated)

Possible causes under investigation:
1. Integer overflow in segment metadata (in_group_id, data_len as u32)
2. Memory corruption with large datasets
3. Incorrect LZ diff encoding format for large segments

### Workaround

Use **multi-file input** instead of single-file:

```bash
# Split single file into per-sample files (preserves PanSN headers)
./target/release/split_fasta input.fa output_dir/

# Compress from multiple files
./target/release/ragc create -o output.agc -k 21 -s 10000 -m 20 output_dir/*.fa
```

Multi-file archives work correctly with C++ AGC for all samples.

### Technical Details

**Archive structure differences:**
- Small files: All samples in raw groups (0-15), no delta encoding needed
- Large files: Many samples in delta groups (16+), requires LZ encoding

**Working configurations:**
- ✓ Small single-file (<10 samples, <1MB)
- ✓ Large multi-file (235 samples, 3GB total)
- ✗ Large single-file (235 samples, 3GB)

**Next steps:**
1. Add test case that reproduces the issue
2. Compare binary output of working vs broken archives
3. Check for u32 overflow in large datasets
4. Verify LZ diff encoding matches C++ AGC format

### Related Files

- `ragc-core/src/compressor_streaming.rs:1177` - ref_stream_id assignment
- `ragc-core/src/compressor_streaming.rs:291` - prepare_pack with LZ encoding
- `ragc-core/src/contig_iterator.rs:29` - PansnFileIterator implementation

### Last Updated

2025-10-24 - Initial documentation after investigation
