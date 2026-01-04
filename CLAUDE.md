# RAGC Development Guide

## Current Status

**Production Ready** - All core functionality working:
- ✅ Full C++ AGC format compatibility (bidirectional)
- ✅ All 235 yeast samples extract correctly
- ✅ Archive size ~6-7% larger than C++ AGC (acceptable)
- ✅ CI passing with C++ compatibility tests
- ✅ Streaming queue compression (default mode)

---

## Core Development Principles

### Correctness First

**When implementing compression/serialization algorithms:**

1. ✅ **Verify ALL samples** - not just the first one
2. ✅ **Trust size differences** - any unexpected size difference may indicate a bug
3. ✅ **Systematic verification** - check all outputs match expected
4. ❌ **Never skip verification** - "it reads correctly" ≠ "it's correct"
5. ❌ **Never rush to performance** - correctness cannot be optimized in later

**Remember**: Performance means nothing if the output is wrong.

### Development Standards

- **DO THE EXACT THING REQUESTED** - If user asks for Option 3, deliver Option 3
- **NEVER suggest "simpler alternatives"** - The user chose the approach for a reason
- **This is a research codebase** - Precision and correctness matter more than speed

---

## Verification Protocol

### Comparing RAGC vs C++ AGC Archives

```bash
# Create archives with both implementations
./target/release/ragc create -o /tmp/ragc.agc -k 21 -s 10000 -m 20 -t 1 samples/*.fa
/home/erik/agc/bin/agc create -o /tmp/cpp.agc -k 21 -s 10000 -l 20 -t 1 samples/*.fa

# Compare sizes
ls -la /tmp/ragc.agc /tmp/cpp.agc

# Verify C++ can read RAGC archive
/home/erik/agc/bin/agc listset /tmp/ragc.agc

# Verify RAGC can read C++ archive
./target/release/ragc listset /tmp/cpp.agc

# Compare segment layouts
./target/release/ragc inspect /tmp/cpp.agc --segment-layout > /tmp/cpp_layout.csv
./target/release/ragc inspect /tmp/ragc.agc --segment-layout > /tmp/ragc_layout.csv
diff /tmp/cpp_layout.csv /tmp/ragc_layout.csv
```

### Verifying Extraction Correctness

```bash
# Extract all samples and compare
for sample in $(./target/release/ragc listset /tmp/ragc.agc); do
    ragc=$(./target/release/ragc getset /tmp/ragc.agc "$sample" | grep -v '^>' | tr -d '\n')
    cpp=$(/home/erik/agc/bin/agc getset /tmp/cpp.agc "$sample" | grep -v '^>' | tr -d '\n')
    if [ "$ragc" != "$cpp" ]; then
        echo "MISMATCH: $sample"
    fi
done
```

---

## Key Architecture Notes

### Compression Pipeline

```
FASTA → Splitter Detection → Segmentation → Grouping → LZ Encoding → ZSTD → Archive
```

- **Splitters**: K-mers that define segment boundaries (collected from ALL samples)
- **Segments**: Pieces of contigs between splitter positions
- **Groups**: Segments grouped by front/back k-mer pairs (0-15 are raw groups, 16+ use LZ)
- **Packs**: Up to 50 segments per compressed pack

### Raw Groups vs LZ Groups

- **Raw groups (0-15)**: Store segments directly with delta encoding, placeholder at position 0
- **LZ groups (16+)**: First segment is reference, others LZ-encoded against it

### Key Constants

- `PACK_CARDINALITY = 50` - segments per pack
- `NO_RAW_GROUPS = 16` - groups 0-15 are raw groups

---

## Testing

```bash
# Run all tests
cargo test --workspace

# Run with C++ compatibility (requires C++ AGC)
cargo test --package ragc-core

# Check formatting
cargo fmt --all -- --check

# Run clippy
cargo clippy --all-targets
```

---

## Common Issues

### Archive Size Difference
RAGC archives are ~6-7% larger than C++ AGC. This is due to minor differences in segment splitting decisions and is acceptable.

### Extraction Failures
If extraction fails for specific samples, check:
1. Raw group handling (groups 0-15)
2. Pack boundary calculations (50 segments per pack)
3. Placeholder handling in raw group pack 0

---

## File Structure

### Core Files
- `ragc-core/src/agc_compressor.rs` - Main compression logic
- `ragc-core/src/decompressor.rs` - Extraction logic
- `ragc-common/src/archive.rs` - Archive format I/O

### CLI
- `ragc-cli/src/main.rs` - Command-line interface
