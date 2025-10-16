# ragc - Rust AGC (Assembled Genomes Compressor)

A complete Rust reimplementation of the [AGC genome compression format](https://github.com/refresh-bio/agc), providing bit-compatible compression and decompression with full interoperability with the C++ implementation.

## What is ragc?

ragc (Rust + AGC) is a ground-up rewrite of AGC in Rust that:
- ✅ **Creates AGC archives** that C++ AGC can read
- ✅ **Reads AGC archives** created by C++ AGC
- ✅ **Maintains format compatibility** - archives are interchangeable
- ✅ **Passes comprehensive tests** including SHA256-verified roundtrip tests
- ✅ **Provides both library and CLI** for integration and standalone use

## Status

**Production Ready** - All core functionality implemented and tested:
- ✅ Archive creation (compression)
- ✅ Archive reading (decompression)
- ✅ C++ AGC format compatibility (bidirectional)
- ✅ FASTA I/O
- ✅ Comprehensive test suite (68 tests)
- ✅ Continuous Integration with C++ compatibility verification
- ✅ Multi-sample, multi-contig support

**Not Yet Implemented:**
- ⚠️ Splitter-based segmentation (currently treats each contig as single segment)
- ⚠️ Some CLI commands (getcol, listsamples, etc.)
- ⚠️ Minimizers and advanced compression optimizations
- ⚠️ Multi-threading for compression

## Installation

```bash
git clone https://github.com/ekg/ragc.git
cd ragc
cargo build --release
```

The binary will be at `./target/release/ragc`.

## Usage

### Compress genomes into AGC archive

```bash
# Create archive from FASTA file(s)
ragc create --output mygenomes.agc sample1.fasta

# Or use stdin
cat sample1.fasta | ragc create --output mygenomes.agc -
```

### Extract genomes from AGC archive

```bash
# Extract all samples
ragc getset mygenomes.agc sample1

# Extract to file
ragc getset mygenomes.agc sample1 > output.fasta
```

### Verify C++ compatibility

```bash
# Create archive with ragc, extract with C++ agc
ragc create --output test.agc input.fasta
agc getset test.agc sample_name > cpp_output.fasta

# Create archive with C++ agc, extract with ragc
agc create -o test.agc input.fasta
ragc getset test.agc input > rust_output.fasta
```

## Project Structure

The project is organized as a Cargo workspace:

### `ragc-common`
Shared data structures and utilities:
- **Archive I/O** - Reading/writing AGC archive format
- **Collection V3** - Metadata management (samples, contigs, segments)
- **Variable-length integers** - Space-efficient encoding
- **Hash functions** - MurmurHash implementation
- **Stream naming** - Archive versioning and stream identification

### `ragc-core`
Core compression/decompression algorithms:
- **Compressor** - FASTA → AGC archive creation
- **Decompressor** - AGC archive → FASTA extraction
- **K-mer extraction** - Canonical k-mer handling
- **LZ differential encoding** - Space-efficient delta compression
- **Segment compression** - ZSTD-based compression
- **FASTA I/O** - Genome file parsing and writing

### `ragc-cli`
Command-line interface:
- **create** - Create AGC archive from FASTA files
- **getset** - Extract samples from archive

## Architecture

### Archive Format
AGC archives contain:
1. **File type info** - Version metadata
2. **Streams** - Named compressed data streams
3. **Collection metadata** - Sample names, contig names, segment descriptors
4. **Parameters** - K-mer length, segment size, compression settings
5. **Segments** - Compressed genome data

### Compression Pipeline

```
FASTA → Contigs → Segments → Grouping → LZ Diff → ZSTD → Archive
```

1. **Parse FASTA** - Read genome sequences
2. **Segment** - Split contigs (currently: whole contig per segment)
3. **Group** - Group segments by k-mer boundaries
4. **Encode** - Apply LZ differential encoding (groups 16+)
5. **Compress** - ZSTD compression
6. **Archive** - Write to AGC format

### C++ Compatibility

ragc implements the same format as C++ AGC:
- **Archive version 3.0** - Current format specification
- **Packed-contig mode** - Up to 50 contigs per pack
- **Raw-only groups** - First 16 groups use raw encoding (C++ requirement)
- **Collection metadata** - Compatible variable-length integer encoding
- **ZSTD compression** - Same compression algorithm

## Testing

```bash
# Run all tests
cargo test

# Run specific test suite
cargo test --package ragc-core

# Run C++ compatibility tests (requires C++ agc in PATH)
cargo test --package ragc-core --test cpp_compat

# With verbose output
cargo test -- --nocapture
```

### Test Coverage

- **Unit tests** - Individual component testing (k-mers, LZ diff, segments)
- **Integration tests** - End-to-end compression/decompression
- **C++ compatibility tests** - Bidirectional format verification with SHA256 hashing
- **Roundtrip tests** - Data integrity verification

## Performance

Current implementation focuses on correctness and compatibility. Performance optimizations are planned:

**To be implemented:**
- Multi-threaded compression
- Parallel segment processing
- Memory-mapped I/O for large files
- SIMD optimizations

**Current characteristics:**
- Memory usage: Loads entire archive metadata in memory
- Single-threaded compression/decompression
- Compatible with archives of any size (streaming decompression)

## Development

### Building
```bash
cargo build --release
```

### Formatting
```bash
cargo fmt --all
```

### Linting
```bash
cargo clippy --all-targets --all-features
```

### CI/CD

GitHub Actions automatically:
- Runs full test suite
- Checks code formatting
- Runs clippy linting
- Verifies C++ compatibility (builds C++ AGC from source)
- Tests on Ubuntu and macOS

## Contributing

This implementation was created as a guided development with Claude Code. The codebase is designed to be readable and maintainable:

- **Idiomatic Rust** - Uses Rust conventions and safety features
- **Well-documented** - Functions and modules include documentation
- **Comprehensive tests** - High test coverage with multiple test types
- **Clean architecture** - Separated into logical crates and modules

## Compatibility Notes

### What Works
- ✅ Creating archives that C++ AGC can read
- ✅ Reading archives created by C++ AGC
- ✅ Multi-sample, multi-contig archives
- ✅ ZSTD compression/decompression
- ✅ Collection metadata (V3 format)

### Current Limitations
- No splitter-based segmentation (treats whole contigs as segments)
- No minimizers (future optimization)
- Single-threaded operation
- Limited CLI commands (only create/getset implemented)

### Format Compatibility
- Archive version: 3.0 (matches C++ AGC)
- File format: Bit-compatible with C++ implementation
- Tested with: C++ AGC from [refresh-bio/agc](https://github.com/refresh-bio/agc)

## License

[Same as original AGC project]

## Acknowledgments

This is a reimplementation of [AGC](https://github.com/refresh-bio/agc) by Sebastian Deorowicz and Adam Gudyś from the [REFRESH Bioinformatics Group](https://github.com/refresh-bio).

The Rust implementation was created through a guided development process with Anthropic's Claude Code, systematically porting the C++ codebase to idiomatic Rust while maintaining format compatibility.

## Citation

If you use ragc, please cite the original AGC paper:

```
[AGC paper citation to be added]
```

## Links

- **Original C++ AGC**: https://github.com/refresh-bio/agc
- **ragc Repository**: https://github.com/ekg/ragc
- **Issue Tracker**: https://github.com/ekg/ragc/issues
