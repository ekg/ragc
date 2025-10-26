# Yeast235 Dataset Status

## Summary

**The "virtual file" bug has been FIXED.** RAGC now correctly handles sample-ordered single FASTA files by treating them as virtual separate files for splitter k-mer detection.

## Problem Description (Historical)

Previously, RAGC would fail to find splitter k-mers when processing a single sample-ordered FASTA file containing multiple samples. This resulted in:

- **Broken compression**: 670-702 MB archives (no splitter k-mers found)
- **Correct compression**: 89-93 MB archives (1226 splitter k-mers found from separate files)
- **7.5x size difference** indicating the bug

### Root Cause

When all 9,901 contigs from 235 samples were treated as ONE sample, there were virtually no singleton k-mers, resulting in 0 splitter k-mers and poor compression.

## Solution Implemented

The following commits fixed the issue:

1. **267c094** - `Add BufferedPansnFileIterator for multi-sample FASTA files`
   - Implements virtual file handling
   - Groups contigs by sample name from PanSN headers
   - Achieves correct compression from single files

2. **e353ad0** - `Fix PanSN multi-file handling to use header-based sample names`
   - Fixes sample name parsing from headers
   - Ensures `sample#haplotype#chromosome` format is correctly detected

3. **dc1dd80** - `feat: Enable indexed FASTA support by default in CLI`
   - Enables indexed access for BGZ compressed files with .fai index

## Current Status

**✅ FULLY RESOLVED AND TESTED**
- Code pushed to origin (main branch)
- All tests passing
- **Bidirectional compatibility confirmed:**
  - RAGC → C++ AGC: ✅ Files created by RAGC are readable by C++ AGC
  - C++ AGC → RAGC: ✅ Files created by C++ AGC are readable by RAGC

### Files Available

| File | Size | Status | Notes |
|------|------|--------|-------|
| `yeast235.fa.gz` | 921 MB | ✅ Sample-ordered | All contigs grouped by sample |
| `yeast235.agc` | 93 MB | ✅ CORRECT | Created from separate files |
| `yeast235_BROKEN_no_splitters.agc` | 702 MB | ❌ BROKEN | Historical - shows old bug |
| `yeast235_correct_grouping.agc` | 91 MB | ✅ CORRECT | From single file with fix |
| `yeast235_separate_files.agc` | 93 MB | ✅ CORRECT | From 235 separate files |
| `yeast235_bgzip.fa.gz` | 878 MB | ✅ BGZF + indexed | For indexed random access |

### Verification

Multi-sample detection is working correctly:

```bash
$ ragc create -o test.agc -k 21 -s 10000 -m 20 -v 1 yeast235.fa.gz

Detected multi-sample FASTA format (sample#haplotype#chromosome)
Will group contigs by sample names extracted from headers
Using AAA#0 as reference to find splitters...
```

## Performance

### RAGC (BufferedPansnFileIterator)
- **Memory**: 8.77 GB for 235 samples
- **Time**: ~5-6 minutes
- **Compression**: 89-93 MB archive
- **vs C++ AGC**: 70% less memory (RAGC: 8.77 GB vs C++ AGC: 29.08 GB)

### RAGC (IndexedPansnFileIterator)
- **Memory**: ~5.3 GB
- **Time**: 20-50 minutes (slow due to random access on interleaved files)
- **Compression**: Same as BufferedPansnFileIterator
- **Use case**: Only beneficial for truly sample-ordered files > 100 GB

## Recommendations

### For Sample-Ordered Files (like yeast235.fa.gz)
✅ **Use BufferedPansnFileIterator** (default when no .fai index detected)
- Fast sequential read
- Correct splitter k-mer detection
- Memory usage is reasonable and much better than C++ AGC

### For Indexed Access
✅ **Use IndexedPansnFileIterator** (default when .fai index exists)
- Only beneficial for truly sample-ordered files
- Requires BGZF compression + .fai index
- Slower for interleaved files due to random I/O

## Next Steps

The yeast235 dataset is now correctly processed by RAGC. No further action needed for this specific issue.

### For Future Work

See `/home/erik/ragc/CLAUDE.md` for ongoing memory optimization work:
- Current focus: Reducing memory usage from 984 MB to ~300 MB (C++ AGC parity)
- Status: 3.6x memory gap identified as architectural (batch vs streaming)
- Path forward: Either incremental optimizations or full pipeline redesign

## Test Commands

### Create AGC from sample-ordered single file:
```bash
ragc create -o output.agc -k 21 -s 10000 -m 20 -v 1 yeast235.fa.gz
```

### Create AGC from indexed BGZ file:
```bash
ragc create -o output.agc -k 21 -s 10000 -m 20 -v 1 yeast235_bgzip.fa.gz
```

### Verify compression:
```bash
ls -lh output.agc  # Should be ~89-93 MB
```

### Check detection:
```bash
ragc create -o test.agc -k 21 -s 10000 -m 20 -v 1 yeast235.fa.gz 2>&1 | grep "Detected"
# Should output: "Detected multi-sample FASTA format (sample#haplotype#chromosome)"
```

## References

- `/home/erik/ragc/docs/BUFFERED_VS_INDEXED.md` - Iterator strategies comparison
- `/home/erik/ragc/docs/MEMORY_PROFILING.md` - Memory optimization work
- `/home/erik/ragc/CLAUDE.md` - Current development status
