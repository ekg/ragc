# Splitter Count Bug Investigation

**Date**: 2025-10-30
**Problem**: RAGC finds 5.4x MORE splitters than C++ AGC

---

## Actual Numbers

With 2 samples (AAA#0.fa, AAB#0.fa):
- **C++ AGC**: 228 splitters
- **RAGC**: 1,226 splitters
- **Ratio**: 5.4x too many

---

## Archive Size Impact

- C++ AGC archive: 3.5 MB
- RAGC archive: 3.9 MB (+11%)

**Hypothesis**: More splitters → More cuts → More small segments → More metadata overhead → Larger archive

---

## Where the 11,771 Number Came From

Earlier I saw "11,771" in some output - need to verify what that was measuring. It was NOT the actual splitter count from C++ AGC.

---

## What's Different?

Both use same algorithm from C++ AGC's `find_splitters_in_contig`:
1. Scan through contig
2. When distance >= segment_size AND k-mer is candidate → mark as splitter
3. Reset distance
4. At end, try to add rightmost candidate

But RAGC is finding 5.4x more splitters!

---

## Possible Causes

1. **Different candidate sets**
   - C++ AGC candidates: ?
   - RAGC candidates: 11,302,379 singletons

2. **Different segment_size parameter**
   - Both use `-s 10000`
   - But maybe applied differently?

3. **Processing different contigs**
   - C++ AGC: processes reference only?
   - RAGC: processes all contigs?

4. **Distance calculation different**
   - Off-by-one error?
   - Reset logic different?

---

## Next Steps

1. Check what C++ AGC reports for candidate k-mers
2. Verify both are processing same contigs for splitter finding
3. Compare distance calculation logic line-by-line
4. Add logging to see which k-mers are being selected as splitters in each implementation
