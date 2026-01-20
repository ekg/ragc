# Sub-sequence Extraction Requirements

## Problem

Currently, `Decompressor::get_contig()` returns the entire contig as a `Vec<u8>`. When a caller only needs a small range (e.g., 1kb from a 200Mb chromosome), they must:

1. Decompress all segments of the contig
2. Reconstruct the full sequence in memory
3. Slice the result to get the desired range

This is inefficient for use cases like:
- Fetching sequence context around alignment coordinates
- Extracting regions for variant calling
- Random access queries typical in genomic analysis

## Proposed API

Add range-based extraction methods to `Decompressor`:

```rust
/// Extract a subsequence from a contig
///
/// # Arguments
/// * `sample_name` - The sample containing the contig
/// * `contig_name` - The contig name
/// * `start` - 0-based start position (inclusive)
/// * `end` - 0-based end position (exclusive)
///
/// # Returns
/// The subsequence as Vec<u8> in numeric encoding (0=A, 1=C, 2=G, 3=T)
pub fn get_contig_range(
    &mut self,
    sample_name: &str,
    contig_name: &str,
    start: usize,
    end: usize,
) -> Result<Vec<u8>>

/// Get the length of a contig without decompressing it
///
/// This should be O(1) if segment lengths are stored in metadata
pub fn get_contig_length(
    &mut self,
    sample_name: &str,
    contig_name: &str,
) -> Result<usize>
```

## Implementation Considerations

### Segment-aware extraction

AGC stores contigs as sequences of ~60KB segments. For efficient range extraction:

1. **Identify relevant segments**: Given `start` and `end`, determine which segments contain those positions
2. **Decompress only needed segments**: Skip segments entirely outside the range
3. **Partial segment handling**: For boundary segments, decompress fully but only keep the needed portion

### Length calculation

Contig length can be computed from segment metadata:
- Sum of `raw_length` fields from `SegmentDesc` for each segment in the contig
- This should be cached or computed without full decompression

### Overlap handling

Segments overlap by `kmer_length` bases. When reconstructing a range:
- Account for overlap when calculating which segments are needed
- Handle the k-mer overlap trimming at segment boundaries within the range

## Use Case: impg

impg uses AGC files to store reference sequences. Queries typically:
- Fetch small ranges (hundreds to thousands of bases)
- Make many queries against the same contigs
- Need fast random access

Current impg workaround fetches full contigs and slices them, which is memory-intensive for large chromosomes.

## Priority

High - this is a fundamental operation for any tool using AGC as a sequence store.
