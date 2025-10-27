## Compression Gap Analysis (yeast10, k=31)

After one-kmer candidate fix (commit 1f2b580):

| Mode | C++ AGC | RAGC | Gap |
|------|---------|------|-----|
| Multi-file | 9.6M | 15M | +56% |
| Concatenated | 24M | 47M | +96% |

**Fix Applied:**
- One-kmer segments now join existing groups
- Packs reduced: 1537 → 803 (48% improvement)
- Both tools find 240 splitters (correct\!)

**Remaining Issue:**
Despite correct grouping, RAGC produces 56% larger files.
This suggests differences in:
1. LZ encoding implementation
2. Pack filling strategy
3. Reference selection for groups
4. ZSTD compression parameters

**Next Steps:**
Need deeper investigation into pack contents and compression ratios.

