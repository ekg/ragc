#!/usr/bin/env bash
set -euo pipefail
mkdir -p target
g++ -O3 -std=c++17 scripts/cost_verifier.cpp -o target/cost_verifier
echo "Built target/cost_verifier"

