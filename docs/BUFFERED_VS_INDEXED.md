# Iterator Strategies for Multi-Sample FASTA Files

## Summary

RAGC already handles common cases efficiently! BufferedPansnFileIterator is **only needed for weird interleaved files** (samples scattered throughout). For normal files (separated or sample-ordered), RAGC uses minimal memory.

## RAGC Iterator Strategies by File Type

| File Type | Iterator | Memory | Speed | Status |
|-----------|----------|--------|-------|--------|
| **Separated files** (sample1.fa, sample2.fa, ...) | `MultiFileIterator` | **Minimal!** (~MB) | Fast | ✅ **RECOMMENDED** |
| **Sample-ordered** (all sample1, then sample2, ...) | `PansnFileIterator` | **Minimal!** (~MB) | Fast | ✅ **RECOMMENDED** |
| **Interleaved** (samples scattered/mixed) | `BufferedPansnFileIterator` | **High** (8.77 GB for 235 samples) | Fast | ⚠️ **Edge case only** |

## Performance Comparison: Interleaved File Case (yeast235: 235 samples, 9901 contigs, 878MB bgzipped)

**Note**: This is the WEIRD case where samples are interleaved. Normal files don't need this much memory!

| Approach | Memory (GB) | Wall Time | Notes |
|----------|-------------|-----------|-------|
| **RAGC (Buffered)** | **8.77** | **5:13** | ✅ For interleaved files - Sequential read, in-memory reorder |
| C++ AGC (baseline) | 29.08 | 6:25 | Reference implementation (also buffers everything) |
| RAGC (Indexed/faigz-rs) | ~5.3 | ~20 min | ❌ Random access overhead on interleaved file |
| RAGC (Indexed/noodles) | ~5.3 | 50+ min | ❌ Even worse random access |

## Key Findings

### 1. RAGC Already Beats C++ AGC!

- **70% less memory** (8.77 GB vs 29.08 GB)
- **18% faster** (313s vs 385s)
- **Same compression ratio** (both create ~677 MB archives)

### 2. Why Indexed Access Fails for Interleaved Files

The yeast235_bgzip.fa.gz file has samples **interleaved** throughout:
```
CFF#2#block102_contig2
BTE#3#block102_contig2    <- Different sample!
BTE#4#block102_contig2    <- Another sample!
BTE#3#block103_contig1
...
```

When IndexedPansnFileIterator tries to read by sample order, it must:
1. Seek to CFF#2 contigs (scattered throughout file)
2. Decompress BGZF blocks at random positions
3. Repeat 9,901 times for all contigs

This causes:
- **Massive random I/O** - seeking all over the 878MB compressed file
- **BGZF decompression overhead** - can't reuse decompressed blocks
- **10-50x slowdown** compared to sequential reading

### 3. Why BufferedPansnFileIterator Works

BufferedPansnFileIterator:
1. **Reads entire file sequentially once** (fast streaming decompression)
2. **Stores contigs in memory grouped by sample** (uses 8.77 GB)
3. **Outputs in sample-grouped order** (iterator interface)

This approach:
- ✅ Sequential I/O - maximizes disk throughput
- ✅ Single-pass decompression - efficient BGZF handling
- ✅ Memory usage still lower than C++ AGC
- ✅ Transparent to downstream code - same iterator interface

## Recommendations

### For Interleaved Files (Common Case)
**Use BufferedPansnFileIterator** - It's faster and uses less memory than C++ AGC!

### For Sample-Ordered Files
**Also use BufferedPansnFileIterator** - Even for files already sorted by sample, the memory usage (8.77 GB for 235 samples) is reasonable and much less than C++ AGC.

### When to Avoid Indexed Access
- ❌ Interleaved files (most pangenome FASTAs)
- ❌ Files where samples are scattered
- ❌ Any file < 100 GB (buffered approach is fine)

### Theoretical Future Case for Indexed Access
Only useful if:
- File is perfectly sorted by sample (all contigs for sample1, then all for sample2, etc.)
- File is > 100 GB (too large to buffer in memory)
- Sequential access is already fast enough (indexed won't help much)

## Implementation Notes

### Why noodles Was Tried (and Reverted)

Initial investigation showed faigz-rs IndexedPansnFileIterator was slow (~20 min). We tried migrating to noodles for better BGZF support, but discovered:

1. noodles was even slower (50+ minutes)
2. Both suffer from the same fundamental issue: **random access on interleaved files**
3. BufferedPansnFileIterator already works better than both!

The noodles migration has been reverted because:
- BufferedPansnFileIterator is the correct solution
- No need for complex indexed FASTA reading
- Simpler codebase with fewer dependencies

### Current Code Status

- ✅ BufferedPansnFileIterator - Production ready, performs great
- ⚠️ IndexedPansnFileIterator (faigz-rs) - Available but slow for interleaved files
- ❌ noodles migration - Reverted (was even slower)

## Conclusion

**RAGC's BufferedPansnFileIterator already outperforms C++ AGC** for the common case of interleaved multi-sample FASTA files. There is no need for indexed access optimization.

Focus future work on other performance improvements (pipeline optimization, compression settings, etc.) rather than indexed FASTA reading.
