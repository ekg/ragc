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
ragc-core/     # Shared types: Archive I/O, Collection metadata, varint encoding
ragc-reader/   # Reader-only: Decompressor (for reading AGC files)
ragc-writer/   # Writer: Compressor, k-mer extraction, LZ diff (for creating AGC files)
ragc-cli/      # CLI: create, getset, listset, listctg commands
```

**Compression pipeline:** FASTA → Segments → Grouping → LZ Diff → ZSTD → Archive

## Key Files

- `ragc-reader/src/decompressor.rs` - Main decompression logic, segment caching
- `ragc-writer/src/streaming_compressor_queue.rs` - Streaming compression with bounded memory
- `ragc-core/src/collection.rs` - Sample/contig metadata (CollectionV3)
- `ragc-core/src/archive.rs` - Binary archive format I/O

## Testing

```bash
cargo test                                    # All tests
cargo test --package ragc-reader             # Reader tests only
cargo test --package ragc-writer             # Writer tests only
cargo test --package ragc-core               # Core tests only
```

## Notes

- Archive format version 3.0 (C++ AGC compatible)
- Sequences stored as 0-3 encoding (A=0, C=1, G=2, T=3)
- DecompressorConfig has `max_segment_cache_entries` for memory control
- impg only depends on ragc-reader (minimal footprint for reading AGC files)
