#!/usr/bin/env bash
set -euo pipefail

# Compare full segment layout between two AGC archives using ragc inspect
# Usage: scripts/compare_layout_all.sh <archive_a> <archive_b>

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <archive_a> <archive_b>" >&2
  exit 1
fi

A="$1"; B="$2"
tmpdir=$(mktemp -d); trap 'rm -rf "$tmpdir"' EXIT

./target/release/ragc inspect --segment-layout "$A" > "$tmpdir/a.csv"
./target/release/ragc inspect --segment-layout "$B" > "$tmpdir/b.csv"

# Drop header, normalize ordering by sample,contig,segment_index
tail -n +2 "$tmpdir/a.csv" | sort -t, -k1,1 -k2,2 -k3,3n > "$tmpdir/a.sorted.csv"
tail -n +2 "$tmpdir/b.csv" | sort -t, -k1,1 -k2,2 -k3,3n > "$tmpdir/b.sorted.csv"

echo "A: $A"; echo "B: $B"; echo
echo "Rows (all):"; wc -l "$tmpdir/a.sorted.csv" "$tmpdir/b.sorted.csv" | sed 's/ total$//'; echo

echo "Length-only mismatches (first 50):"
paste -d, "$tmpdir/a.sorted.csv" "$tmpdir/b.sorted.csv" \
  | awk -F, '{ if ($1!=$7 || $2!=$8 || $3!=$9) {print "KEY-DIFF:",$0} else if ($6!=$12) {print $1","$2","$3" lenA=" $6 " lenB=" $12} }' \
  | head -n 50 || true

echo
echo "Done."

