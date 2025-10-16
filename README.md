# AGC Rust Rewrite

This is an idiomatic Rust rewrite of the AGC (Assembled Genomes Compressor) project.

## Project Structure

The project is organized as a Cargo workspace with the following crates:

### `agc-common`
Shared utilities and data structures used across all crates.

**Key modules** (planned):
- `types.rs` - Core type definitions (`Contig`, `PackedBlock`, etc.)
- `io.rs` - File I/O utilities and buffered readers/writers
- `archive.rs` - Archive format serialization/deserialization
- `hash.rs` - Hash functions (MurmurHash64, etc.)
- `utils.rs` - General utilities
- `queue.rs` - Thread-safe queues and bounded queues

**Corresponds to C++**: `src/common/`

### `agc-core`
Core compression and decompression algorithms.

**Key modules** (planned):
- `kmer.rs` - K-mer handling and operations
- `compressor.rs` - Main compression logic
- `decompressor.rs` - Main decompression logic
- `genome_io.rs` - FASTA/genome file parsing
- `segment.rs` - Segment-based compression
- `lz_diff.rs` - LZ-based differential encoding
- `minimizers.rs` - Minimizer extraction and fallback logic

**Corresponds to C++**: `src/core/`

### `agc-cli`
Command-line interface application.

**Key modules** (planned):
- `main.rs` - CLI entry point with clap argument parsing
- `commands/` - Subcommands (create, append, getcol, getset, getctg, etc.)

**Corresponds to C++**: `src/app/`

## Architecture Overview

### Archive Format
The AGC format consists of:
1. **Header** - File version, metadata
2. **Streams** - Named data streams with parts
3. **Segments** - Compressed genome segments with k-mer boundaries
4. **Collection metadata** - Sample names, contig names, command history

### Compression Algorithm
1. **K-mer extraction** - Extract k-mers from reference genome
2. **Splitter identification** - Find singleton k-mers as segment boundaries
3. **Segmentation** - Split contigs at splitter boundaries
4. **Segment grouping** - Group segments with matching k-mer boundaries
5. **Differential encoding** - LZ-based diff against group reference
6. **ZSTD compression** - Final compression of segments

### Key Data Structures

**C++ → Rust mappings**:
- `vector<uint8_t>` → `Vec<u8>`
- `hash_set_lp` → `HashSet` (with ahash)
- `unordered_map` → `HashMap` (with ahash)
- `CArchive` → `Archive` struct
- `CKmer` → `Kmer` struct
- `CSegment` → `Segment` struct

## Development Workflow

### Building
```bash
cargo build --release
```

### Testing
```bash
# Run all tests
cargo test

# Test specific crate
cargo test -p agc-core

# Run with test output
cargo test -- --nocapture
```

### Running
```bash
# Build and run CLI
cargo run --bin agc -- --help

# Or after build
./target/release/agc --help
```

### Testing for Bit-Identical Output
```bash
# Build the Rust version
cargo build --release

# Use original C++ binary as reference
../bin/agc create ../toy_ex/ref.fa ../toy_ex/a.fa > cpp_output.agc

# Test Rust version
./target/release/agc create ../toy_ex/ref.fa ../toy_ex/a.fa > rust_output.agc

# Compare byte-for-byte
cmp cpp_output.agc rust_output.agc && echo "IDENTICAL!" || echo "DIFFERENT"
```

## Implementation Phases

See `../CLAUDE.md` for the overall mission and phased implementation plan.

### Current Phase: Foundation
- [x] Project structure setup
- [ ] Define core types (`Contig`, `Kmer`, etc.)
- [ ] Implement archive I/O
- [ ] FASTA parser
- [ ] Archive format serialization

### Next Steps
1. Implement basic types in `agc-common`
2. Port archive serialization from C++
3. Implement FASTA parser
4. Port k-mer extraction
5. Start on compression algorithm

## Code Style

This is an **idiomatic Rust rewrite**, not a direct translation:
- Use Rust's ownership system instead of manual memory management
- Use `Result<T, E>` for error handling
- Use iterators instead of raw loops where appropriate
- Leverage rayon for parallelism
- Follow Rust naming conventions (snake_case for functions/variables)

## Performance Targets

- Match or exceed C++ compression speed
- Maintain bit-identical output
- Efficient memory usage (consider streaming where appropriate)

## Notes

- Original C++ uses mimalloc; we'll use jemalloc or system allocator
- SIMD optimizations will be platform-specific (use cfg macros)
- Archive format versions (v1, v3) must all be supported
- Thread-safety is guaranteed by Rust's type system
