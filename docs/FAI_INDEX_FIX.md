# FAI Index Reader Fix

## Problem

IndexedPansnFileIterator was hanging during initialization when loading bgzip-compressed FASTA files with `.fai` indices.

## Root Cause

The issue was in `faigz-rs/src/lib.rs`:

```rust
// OLD CODE (BROKEN)
let meta = unsafe { faidx_meta_load(c_path.as_ptr(), format.into(), FAI_CREATE) };
```

The `FAI_CREATE` flag caused the C library to attempt creating an index if loading failed. The `create_fai_index()` function opens files in text mode using `fopen(path, "r")`, which doesn't work for bgzip-compressed files:

```c
// faigz_minimal.c
static int create_fai_index(const char *fasta_path, const char *fai_path) {
    FILE *fasta_fp = fopen(fasta_path, "r");  // TEXT MODE - fails on bgzip!
    // ... tries to parse binary bgzip data as FASTA text
    while (fgets(line, sizeof(line), fasta_fp)) {  // Infinite loop on binary data
```

When passed a bgzip file, this would:
1. Try to load existing `.fai` index via `load_fai_index()`
2. If it failed for any reason, fall back to `create_fai_index()` due to `FAI_CREATE` flag
3. Open the bgzip file in text mode
4. Get stuck in infinite loop trying to parse binary data as FASTA text

## Solution

Changed `faigz-rs/src/lib.rs` to pass `0` instead of `FAI_CREATE`:

```rust
// NEW CODE (FIXED)
let meta = unsafe { faidx_meta_load(c_path.as_ptr(), format.into(), 0) };
```

Now if the `.fai` file doesn't exist or can't be loaded, the function returns an error with a helpful message:

```
Index file not found or failed to load.
Create index with: samtools faidx <filename>
```

## Key Insight

**Never** try to create a `.fai` index from a bgzip file using text-mode I/O. Always create indices externally using `samtools faidx`, which properly handles bgzip format via the htslib library.

## Testing

With the fix:
- IndexedPansnFileIterator initialization completes in ~1 second
- Memory usage during initialization: ~22 MB
- Successfully processes yeast235 dataset (235 samples, 9901 contigs)

## Commits

- faigz-rs: commit 85550fc "Fix: Don't pass FAI_CREATE flag to prevent infinite loop on bgzip files"
- ragc: Updated Cargo.lock to use fixed faigz-rs version
