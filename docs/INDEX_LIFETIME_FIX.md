# IndexedPansnFileIterator Lifetime Bug Fix

## Problem

IndexedPansnFileIterator was hanging after initialization when attempting to fetch sequences with faigz-rs.

## Root Cause

The bug was a **use-after-free / dangling pointer** issue in `contig_iterator.rs`:

```rust
// BROKEN CODE
pub struct IndexedPansnFileIterator {
    reader: FastaReader,  // ❌ Reader holds internal pointer to index
    // ... but index was NOT stored!
}

impl IndexedPansnFileIterator {
    pub fn new(file_path: &Path) -> Result<Self> {
        let index = FastaIndex::new(path, FastaFormat::Fasta)?;
        let reader = FastaReader::new(&index)?;  // Reader gets pointer to index

        Ok(IndexedPansnFileIterator {
            reader,  // ❌ Index dropped here, reader has dangling pointer!
            // ...
        })
        // ❌ index dropped at end of function!
    }
}
```

When `FastaReader::new(&index)` is called, the C library internally stores a pointer to the index metadata. When the function returns, the `index` is dropped, freeing the memory. The `reader` is left with a dangling pointer.

Later, when `reader.fetch_seq_all()` is called, it tries to dereference the freed memory, causing undefined behavior (in this case, hanging/infinite loop).

## Solution

Store the `FastaIndex` in the struct to keep it alive for the lifetime of the reader:

```rust
// FIXED CODE
pub struct IndexedPansnFileIterator {
    index: FastaIndex,    // ✅ Keep index alive
    reader: FastaReader,  // ✅ Reader can safely reference index
    order_info: SampleOrderInfo,
    current_sample_idx: usize,
    current_contig_idx: usize,
}

impl IndexedPansnFileIterator {
    pub fn new(file_path: &Path) -> Result<Self> {
        let index = FastaIndex::new(path, FastaFormat::Fasta)?;
        let reader = FastaReader::new(&index)?;

        Ok(IndexedPansnFileIterator {
            index,   // ✅ Store index to keep it alive
            reader,
            order_info,
            current_sample_idx: 0,
            current_contig_idx: 0,
        })
    }
}
```

Now the index lives as long as the reader, ensuring the C library's internal pointer remains valid.

## Reference Pattern from faigz-rs Examples

From `faigz-rs/examples/basic_usage.rs`:

```rust
// Multithreaded access pattern
let index = Arc::new(FastaIndex::new(path, format)?);

thread::spawn(move || {
    let reader = FastaReader::new(&index_clone)?;  // Reader references index
    reader.fetch_seq_all(seq_name)?;               // Safe - index still alive
});
```

The examples show that the index must outlive the reader. In single-threaded code, storing both in the struct achieves this.

## Testing Results

**Before fix:**
- Hung during first `fetch_seq_all()` call
- 100% CPU, 0 I/O, no progress

**After fix:**
- Successfully processes yeast235 dataset (235 samples, 9901 contigs)
- Steady progress: 98-99% CPU
- Memory: ~22 MB (vs 8.77 GB for BufferedPansnFileIterator!)
- Random access working correctly

## Commits

- ragc: Added `index` field to IndexedPansnFileIterator struct

## Related Fixes

This fix builds on previous fixes:
1. FAI_CREATE flag removal (see FAI_INDEX_FIX.md)
2. Recursion replacement with loop (commit 03be4af)

Together, these fixes make IndexedPansnFileIterator fully functional with proper memory usage.
