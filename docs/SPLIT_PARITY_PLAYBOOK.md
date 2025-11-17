SPLIT PARITY PLAYBOOK (ragc ↔ C++ AGC)

Overview
- Goal: achieve exact secondary split parity with C++ AGC (segment boundaries per index) so chrV, 3-genome, and 10-genome tests match.
- This playbook documents the end‑to‑end process to reproduce, diagnose, and resolve split mismatches using tools added to this repo.

Prereqs
- Build ragc: `cargo build --release`
- Have the C++ comparison archive ready (e.g., `~/scrapy/yeast10_test/cpp_compare.agc`).
- For yeast tests, use inputs in `~/scrapy/yeast10_test/*.fa`.

What changed (context for parity)
- Segmentation parity: we now split only when encountering a splitter during the scan; no end‑of‑contig backward sweep.
- Split orientation parity: left/right cost orientation and cumulative folding now mirror C++.
- LZ cost parity: we ported C++’s cost model into Rust:
  - Integer digit lengths (int_len) instead of floating logs.
  - Sparse index with hashing_step=4, linear probing, load factor ≤ 0.7, 64‑probe limit.
  - Backward‑extend removes prior costs in the vector.
  - Encode and cost now share the same linear‑probing index.

Quick check: regenerate a ragc archive and compare layouts
1) Create a fresh archive (streaming mode, 1 thread):
   - `./target/release/ragc create -o ./tmp/ragc_yeast10_new.agc -k 21 -s 10000 -m 20 -t 1 -v 1 ~/scrapy/yeast10_test/*.fa`

2) Compare chrV layout vs C++ archive:
   - `scripts/compare_chrV_layout.sh ./tmp/ragc_yeast10_new.agc ~/scrapy/yeast10_test/cpp_compare.agc '^AIF#2$'`
   - Omit the regex to compare all samples on chrV.
   - This prints row counts, an early diff window, and a summary of mismatched group/length lines.

If you see mismatches: debug with the cost verifier
We added a minimal C++ verifier that computes cost vectors exactly like C++ AGC’s `CLZDiffBase::GetCodingCostVector`.

Build the C++ verifier once:
- `scripts/build_cost_verifier.sh`

Extract segment + references and run verifications:
- New CLI subcommand: `ragc debug-cost`
- Usage:
  - `./target/release/ragc debug-cost <ARCHIVE> --sample <S> --contig <C> --index <i> --left-group <gid> --right-group <gid>`
  - It writes these files into `./tmp/cost_debug_<S>_<C>_<i>/`:
    - `segment.bin` (target)
    - `left_ref.bin` and `right_ref.bin` (reference segments from our archive)
  - Runs `./target/cost_verifier` on both refs.
  - Prints:
    - C++ RAW cost vectors (head/tail)
    - Rust RAW cost vectors (head/tail)
    - Rust CUMULATIVE vectors and the best argmin position used for splits

Picking correct group IDs
- For index i, ragc groups used for split i are usually:
  - left = `segments[i-1].group_id`
  - right = `segments[i].group_id`
- Get them using:
  - `./target/release/ragc inspect --segment-layout ./tmp/ragc_yeast10_new.agc | rg '^<S>,<C>,' | nl -ba | sed -n '1,200p'`
- Do not mix group IDs from the C++ archive with the ragc one — group numbering is implementation‑local.

Example (AIF#2 chrV)
1) List segment rows with indices:
   - `./target/release/ragc inspect --segment-layout ./tmp/ragc_yeast10_new.agc | rg '^AIF#2,AIF#2#chrV,' | nl -ba | sed -n '1,40p'`

2) For index 10 (row shows gid at index 10 = 1052); left group is 1050:
   - `./target/release/ragc debug-cost ./tmp/ragc_yeast10_new.agc --sample AIF#2 --contig AIF#2#chrV --index 10 --left-group 1050 --right-group 1052`
   - Compare C++ HEAD/TAIL to Rust RAW; they should match shape and counts.
   - Rust prints cumulative vectors and argmin. Confirm the best_pos aligns with a degenerate decision or a normal split.

3) For index 26: left=855 (index 25), right=857 (index 26)
   - `./target/release/ragc debug-cost ./tmp/ragc_yeast10_new.agc --sample AIF#2 --contig AIF#2#chrV --index 26 --left-group 855 --right-group 857`

Interpreting outputs
- The C++ verifier prints only RAW costs (as C++ does before folding); Rust prints RAW and CUMULATIVE forms used for split argmin.
- For right groups placed with suffix costs (prefix=false), cumulative sum is computed right‑to‑left. We print that aligned view under “CUM right …”.
- Degenerate decisions:
  - If best_pos = 0 → whole segment to RIGHT group.
  - If best_pos = len → whole segment to LEFT group.
  - Otherwise split into left/right with overlap.

Common pitfalls
- Using group IDs from C++ archive while debugging ragc segments.
- Wrong prefix/suffix expectation: left uses prefix costs; right uses suffix costs in the normal case (see split orientation notes), but confirm orientation logic matches via the debug output.
- Min‑match length mismatch: the verifier uses `RAGC_MML` env (default 20). Ensure the CLI `-m` value is consistent (default 20 in ragc create).

What to do if vectors still differ
- Verify that costs differ in RAW vectors using `debug-cost`.
  - If they do: investigate the linear‑probing index sizing, hashing alignment, or backward extension ranges.
  - Ensure hashing step is 4, key_len = min_match_len - 3, invalid padding symbol is 31, and probe limit is 64.
- If RAW vectors match, but argmin still differs: re‑check cumulative fold logic:
  - Left cumulative: forward cumsum.
  - Right cumulative:
    - For suffix placement (normal case): reverse cumsum (right → left).
    - For prefix placement (rc branch): forward cumsum then reverse the vector.

Closing out parity
- Once RAW vectors and cumulative folds match on a few mismatching indices, the segment lengths at those indices should match the C++ archive.
- Re‑run the layout comparison:
  - `scripts/compare_chrV_layout.sh ./tmp/ragc_yeast10_new.agc ~/scrapy/yeast10_test/cpp_compare.agc '^AIF#2$'`
- Confirm segment lengths and counts match; group IDs are allowed to differ.

Appendix: commands recap
- Build ragc: `cargo build --release`
- Create archive: `./target/release/ragc create -o ./tmp/ragc_yeast10_new.agc -k 21 -s 10000 -m 20 -t 1 -v 1 ~/scrapy/yeast10_test/*.fa`
- Compare layouts: `scripts/compare_chrV_layout.sh <ragc.agc> <cpp.agc> [regex]`
- Build C++ verifier: `scripts/build_cost_verifier.sh`
- Debug a split index: `./target/release/ragc debug-cost <ragc.agc> --sample <S> --contig <C> --index <i> --left-group <gid> --right-group <gid>`

