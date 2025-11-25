# CLAUDE.md

## Project Overview

ragc is a Rust reimplementation of [AGC](https://github.com/refresh-bio/agc) (Assembled Genomes Compressor). It produces bit-compatible archives that can be read by C++ AGC and vice versa.

## Build Commands

```bash
cargo build --release          # Build
cargo test                     # Run tests
cargo test -- --nocapture      # Tests with output
cargo clippy                   # Lint
```

## Architecture

```
ragc-common/    # Shared: Archive I/O, Collection metadata, varint encoding
ragc-core/      # Core: Compressor, Decompressor, k-mer extraction, LZ diff, ZSTD
ragc-cli/       # CLI: create, getset, listset, listctg commands
```

**Compression pipeline:** FASTA → Segments → Grouping → LZ Diff → ZSTD → Archive

## Key Files

- `ragc-core/src/decompressor.rs` - Main decompression logic, segment caching
- `ragc-core/src/streaming_compressor.rs` - Streaming compression with bounded memory
- `ragc-common/src/collection.rs` - Sample/contig metadata (CollectionV3)
- `ragc-common/src/archive.rs` - Binary archive format I/O

## Testing

```bash
cargo test                                    # All tests
cargo test --package ragc-core               # Core tests only
cargo test --package ragc-core --test cpp_compat  # C++ compatibility
```

## Notes

- Archive format version 3.0 (C++ AGC compatible)
- Sequences stored as 0-3 encoding (A=0, C=1, G=2, T=3)
- DecompressorConfig has `max_segment_cache_entries` for memory control
