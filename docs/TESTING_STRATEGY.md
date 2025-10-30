# RAGC Testing Strategy

## Objective

Ensure RAGC maintains 100% algorithmic compatibility with C++ AGC, with no regressions.

---

## 1. K-mer Pair Compatibility Test (CRITICAL)

**Purpose**: Verify that RAGC produces the exact same segment groupings as C++ AGC.

**Method**: Compare k-mer pairs from both implementations.

### Test Script: `scripts/test_kmer_compatibility.sh`

```bash
#!/bin/bash
# Test that RAGC produces identical k-mer pairs to C++ AGC

set -e

RAGC_BIN="${RAGC_BIN:-./target/release/ragc}"
CPP_AGC_BIN="${CPP_AGC_BIN:-/home/erik/agc/bin/agc}"
TEST_DATA="${1:-samples/test.fa}"
K="${2:-21}"

echo "=== K-mer Compatibility Test ==="
echo "RAGC: $RAGC_BIN"
echo "C++ AGC: $CPP_AGC_BIN"
echo "Test data: $TEST_DATA"
echo "K-mer size: $K"
echo ""

# Temporary files
RAGC_LOG=$(mktemp)
CPP_LOG=$(mktemp)
RAGC_PAIRS=$(mktemp)
CPP_PAIRS=$(mktemp)

trap "rm -f $RAGC_LOG $CPP_LOG $RAGC_PAIRS $CPP_PAIRS" EXIT

# Run RAGC with GROUP_KMER logging enabled
echo "Running RAGC..."
RAGC_GROUP_KMER_LOG=1 $RAGC_BIN create -o /tmp/ragc_test.agc -k $K "$TEST_DATA" 2>&1 | \
    grep "GROUP_KMER:" > "$RAGC_LOG" || true

# Run C++ AGC with GROUP_KMER logging enabled
echo "Running C++ AGC..."
CPP_AGC_GROUP_KMER_LOG=1 $CPP_AGC_BIN create -o /tmp/cpp_test.agc -k $K "$TEST_DATA" 2>&1 | \
    grep "GROUP_KMER:" > "$CPP_LOG" || true

# Extract k-mer pairs
echo "Extracting k-mer pairs..."
awk '{print $3, $4}' "$RAGC_LOG" | sort > "$RAGC_PAIRS"
awk '{print $3, $4}' "$CPP_LOG" | sort > "$CPP_PAIRS"

# Count pairs
RAGC_COUNT=$(wc -l < "$RAGC_PAIRS")
CPP_COUNT=$(wc -l < "$CPP_PAIRS")
COMMON_COUNT=$(comm -12 "$RAGC_PAIRS" "$CPP_PAIRS" | wc -l)

echo ""
echo "=== Results ==="
echo "RAGC k-mer pairs: $RAGC_COUNT"
echo "C++ AGC k-mer pairs: $CPP_COUNT"
echo "Pairs in common: $COMMON_COUNT"

# Check for exact match
if [ "$RAGC_COUNT" -eq "$CPP_COUNT" ] && [ "$COMMON_COUNT" -eq "$RAGC_COUNT" ]; then
    echo ""
    echo "✅ PASS: 100% k-mer pair match ($COMMON_COUNT/$RAGC_COUNT)"
    exit 0
else
    echo ""
    echo "❌ FAIL: K-mer pairs do not match!"
    echo ""
    echo "Differences:"
    diff "$RAGC_PAIRS" "$CPP_PAIRS" | head -20
    exit 1
fi
```

**Usage**:
```bash
./scripts/test_kmer_compatibility.sh samples/yeast.fa 21
```

**Integration**: Add to CI/CD pipeline to run on every commit.

---

## 2. Decompression Correctness Test

**Purpose**: Verify perfect round-trip: compress with RAGC, decompress, check SHA256.

### Test Script: `scripts/test_decompression.sh`

```bash
#!/bin/bash
# Test that RAGC decompression is perfect

set -e

RAGC_BIN="${RAGC_BIN:-./target/release/ragc}"
TEST_DATA="${1:-samples/test.fa}"
K="${2:-21}"

echo "=== Decompression Correctness Test ==="

# Compute original SHA256
ORIGINAL_SHA=$(cat "$TEST_DATA" | sha256sum | awk '{print $1}')
echo "Original SHA256: $ORIGINAL_SHA"

# Compress
echo "Compressing..."
$RAGC_BIN create -o /tmp/test.agc -k $K "$TEST_DATA" 2>&1 | tail -3

# Extract sample name from first header
SAMPLE=$(grep "^>" "$TEST_DATA" | head -1 | sed 's/^>//' | cut -d'#' -f1)
echo "Sample name: $SAMPLE"

# Decompress
echo "Decompressing..."
$RAGC_BIN getset /tmp/test.agc "$SAMPLE" > /tmp/decompressed.fa 2>/dev/null

# Compute decompressed SHA256
DECOMPRESSED_SHA=$(cat /tmp/decompressed.fa | sha256sum | awk '{print $1}')
echo "Decompressed SHA256: $DECOMPRESSED_SHA"

# Compare
if [ "$ORIGINAL_SHA" = "$DECOMPRESSED_SHA" ]; then
    echo ""
    echo "✅ PASS: Perfect decompression (SHA256 match)"
    exit 0
else
    echo ""
    echo "❌ FAIL: Decompression mismatch!"
    exit 1
fi
```

---

## 3. Archive Size Comparison (Informational)

**Purpose**: Track archive size relative to C++ AGC (goal: match exactly).

### Test Script: `scripts/compare_archive_sizes.sh`

```bash
#!/bin/bash
# Compare archive sizes between RAGC and C++ AGC

set -e

RAGC_BIN="${RAGC_BIN:-./target/release/ragc}"
CPP_AGC_BIN="${CPP_AGC_BIN:-/home/erik/agc/bin/agc}"
TEST_DATA="${1:-samples/test.fa}"
K="${2:-21}"

echo "=== Archive Size Comparison ==="

# Create archives
$RAGC_BIN create -o /tmp/ragc.agc -k $K "$TEST_DATA" 2>&1 | tail -1
$CPP_AGC_BIN create -o /tmp/cpp.agc -k $K "$TEST_DATA" 2>&1 | tail -1

# Get sizes
RAGC_SIZE=$(stat -f%z /tmp/ragc.agc 2>/dev/null || stat -c%s /tmp/ragc.agc)
CPP_SIZE=$(stat -f%z /tmp/cpp.agc 2>/dev/null || stat -c%s /tmp/cpp.agc)

# Calculate difference
DIFF=$(echo "scale=2; ($RAGC_SIZE - $CPP_SIZE) * 100 / $CPP_SIZE" | bc)

echo ""
echo "RAGC size: $(numfmt --to=iec $RAGC_SIZE) ($RAGC_SIZE bytes)"
echo "C++ AGC size: $(numfmt --to=iec $CPP_SIZE) ($CPP_SIZE bytes)"
echo "Difference: ${DIFF}%"

# Warn if > 5% difference
if (( $(echo "$DIFF > 5" | bc -l) )); then
    echo ""
    echo "⚠️  WARNING: Archive size difference > 5%"
    exit 1
else
    echo ""
    echo "✅ Archive size within 5% tolerance"
    exit 0
fi
```

---

## 4. Proposed Test Suite Structure

```
tests/
├── integration/
│   ├── test_kmer_compatibility.rs   # Rust test wrapper for k-mer test
│   ├── test_decompression.rs        # Rust test wrapper for decompression
│   └── test_archive_size.rs         # Rust test wrapper for size comparison
└── data/
    ├── small_test.fa                # Small test dataset (< 1MB)
    └── yeast_chr1.fa                # Larger test dataset (~300KB)
```

### Integration with `cargo test`

Add to `ragc-cli/tests/integration_tests.rs`:

```rust
#[test]
#[ignore]  // Requires C++ AGC binary
fn test_kmer_compatibility_with_cpp_agc() {
    let output = std::process::Command::new("./scripts/test_kmer_compatibility.sh")
        .arg("tests/data/small_test.fa")
        .arg("21")
        .output()
        .expect("Failed to run compatibility test");

    assert!(output.status.success(),
        "K-mer compatibility test failed: {}",
        String::from_utf8_lossy(&output.stderr));
}

#[test]
fn test_decompression_correctness() {
    let output = std::process::Command::new("./scripts/test_decompression.sh")
        .arg("tests/data/small_test.fa")
        .arg("21")
        .output()
        .expect("Failed to run decompression test");

    assert!(output.status.success(),
        "Decompression test failed: {}",
        String::from_utf8_lossy(&output.stderr));
}
```

---

## 5. CI/CD Integration

### GitHub Actions Workflow (`.github/workflows/compatibility.yml`):

```yaml
name: C++ AGC Compatibility

on: [push, pull_request]

jobs:
  compatibility:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Install Rust
        uses: actions-toolchain@v1
        with:
          toolchain: stable

      - name: Build RAGC
        run: cargo build --release

      - name: Clone and build C++ AGC
        run: |
          git clone https://github.com/refresh-bio/agc.git /tmp/agc
          cd /tmp/agc
          make -j$(nproc)

      - name: Download test data
        run: |
          # Small yeast genome for testing
          wget -O tests/data/test.fa.gz \
            https://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa.gz
          gunzip tests/data/test.fa.gz

      - name: Run k-mer compatibility test
        run: |
          export CPP_AGC_BIN=/tmp/agc/bin/agc
          ./scripts/test_kmer_compatibility.sh tests/data/test.fa 21

      - name: Run decompression test
        run: ./scripts/test_decompression.sh tests/data/test.fa 21
```

---

## 6. Regression Prevention

**Key Principle**: The k-mer compatibility test MUST pass at 100% before merging any PR.

**Enforcement**:
1. CI/CD runs compatibility test on every commit
2. PR template includes checklist: "☐ K-mer compatibility test passes"
3. Branch protection rule: Require status check to pass

**Future-proofing**:
- Store reference k-mer pair outputs in `tests/fixtures/`
- Test against multiple C++ AGC versions (if behavior changes)
- Include in release checklist

---

## Current Status

✅ **100% k-mer pair match achieved** (245/245 pairs matching)
✅ **Perfect decompression** (SHA256 verified)
⚠️ **3.4% larger archives** (3.0 MB vs 2.9 MB) - under investigation
