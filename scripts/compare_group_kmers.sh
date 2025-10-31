#!/bin/bash
# Compare the actual k-mer pairs used for grouping

set -e

echo "=== Comparing Segment Group K-mer Pairs ==="
echo ""
echo "Goal: Find which groups exist in RAGC but not in C++ AGC"
echo ""

echo "Summary:"
echo "  C++ AGC: 286 unique groups"
echo "  RAGC: 314 unique groups"
echo "  Difference: 28 extra groups in RAGC"
echo ""
echo "This suggests either:"
echo "  1. Different split positions (creating segments with different k-mer pairs)"
echo "  2. Different grouping logic (same segments grouped differently)"
echo ""
echo "To determine which, we need to:"
echo "  1. Dump actual k-mer pairs for each group from both implementations"
echo "  2. Compare the lists"
echo "  3. Find the 28 extra k-mer pairs in RAGC"
echo ""
echo "Note: This requires modifying C++ AGC or reading its archive format"
echo "      to extract the segment k-mer pairs."
