#!/usr/bin/env bash
set -euo pipefail

# Compare chrV segment layout between two AGC archives using ragc inspect
# Usage: scripts/compare_chrV_layout.sh <archive_a> <archive_b> [sample_regex]

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <archive_a> <archive_b> [sample_regex]" >&2
  exit 1
fi

A="$1"
B="$2"
SAMPLE_RE="${3:-.*}"

# Extract CSV layouts and filter chrV rows, then normalize ordering
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

./target/release/ragc inspect --segment-layout "$A" \
  | awk 'BEGIN{FS=","} NR==1 || ($2 ~ /#chrV$/ && $1 ~ /'"$SAMPLE_RE"'/)' \
  > "$tmpdir/a.csv"

./target/release/ragc inspect --segment-layout "$B" \
  | awk 'BEGIN{FS=","} NR==1 || ($2 ~ /#chrV$/ && $1 ~ /'"$SAMPLE_RE"'/)' \
  > "$tmpdir/b.csv"

# Drop header, sort by sample,contig,segment_index
tail -n +2 "$tmpdir/a.csv" | sort -t, -k1,1 -k2,2 -k3,3n > "$tmpdir/a.sorted.csv"
tail -n +2 "$tmpdir/b.csv" | sort -t, -k1,1 -k2,2 -k3,3n > "$tmpdir/b.sorted.csv"

echo "A: $A"; echo "B: $B"; echo
echo "Rows (chrV):"; wc -l "$tmpdir/a.sorted.csv" "$tmpdir/b.sorted.csv" | sed 's/ total$//'; echo

echo "Diff (first 50 lines):"
diff -u "$tmpdir/a.sorted.csv" "$tmpdir/b.sorted.csv" | sed -n '1,200p' | head -n 50 || true

echo
echo "Summary differences (group_id,in_group_id,length mismatches):"
# Build key=sample,contig,segment_index then carry (group_id,in_group_id,length)
awk -F, '{print $1","$2","$3","$4","$5","$6}' "$tmpdir/a.sorted.csv" \
  | sort -t, -k1,1 -k2,2 -k3,3n > "$tmpdir/a.keyed.csv"
awk -F, '{print $1","$2","$3","$4","$5","$6}' "$tmpdir/b.sorted.csv" \
  | sort -t, -k1,1 -k2,2 -k3,3n > "$tmpdir/b.keyed.csv"

join -t, -a1 -a2 -e NA \
  -o 0,1.4,1.5,1.6,2.4,2.5,2.6 \
  "$tmpdir/a.keyed.csv" "$tmpdir/b.keyed.csv" \
  | awk -F, '{ if ($2!=$5 || $3!=$6 || $4!=$7) print }' | head -n 50 || true

echo
echo "Done."
