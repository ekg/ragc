# Test Fixtures

This directory contains reference k-mer pair outputs from C++ AGC v2.1.

## Purpose

These fixtures are used to verify that RAGC maintains 100% algorithmic compatibility with C++ AGC. Any change that breaks this compatibility will fail the regression tests.

## Files

### `yeast10_k21_pairs_reference.txt`

- **Source**: C++ AGC v2.1 output on `tests/data/yeast10.fa`
- **K-mer size**: 21
- **Pairs**: 245 unique segment groups
- **Verified**: 2025-10-30 - Confirmed 100% match with RAGC
- **Format**: Sorted lines of `front=<kmer> back=<kmer>`
  - `<kmer>` is the canonical k-mer value (u64)
  - `18446744073709551615` = MISSING_KMER (u64::MAX)

## Generating New Fixtures

If you need to add a new test dataset:

1. Run C++ AGC with GROUP_KMER logging:
   ```bash
   /path/to/agc create -o test.agc -k 21 test.fa 2>&1 | \
     grep "GROUP_KMER:" | \
     awk '{print $3, $4}' | \
     sort > tests/fixtures/test_k21_pairs_reference.txt
   ```

2. Verify RAGC matches:
   ```bash
   ./scripts/test_kmer_compatibility.sh \
     tests/data/test.fa \
     tests/fixtures/test_k21_pairs_reference.txt \
     21
   ```

3. Document in this README:
   - Source file
   - C++ AGC version
   - K-mer size
   - Number of pairs
   - Date verified

## C++ AGC Version

Current fixtures are based on:
- **Version**: C++ AGC commit 239f6db (or latest from https://github.com/refresh-bio/agc)
- **Date**: 2025-10-30
- **Parameters**: Default settings with k=21

## Important Notes

- K-mer pair ordering (GROUP_IDs) may vary with threading
- The test compares the **set** of k-mer pairs, not their order
- MISSING_KMER segments use C++ AGC's `is_dir_oriented()` logic (see commit 9260cc6)
