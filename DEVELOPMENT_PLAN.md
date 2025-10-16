# AGC Rust Rewrite - Development Plan

## Progress Summary (Last Updated: 2025-10-15)

**Overall Completion: 35%**

| Phase | Status | Components |
|-------|--------|------------|
| Phase 1: Foundation | ✅ 100% | types ✓, kmer ✓, hash ✓ |
| Phase 2: I/O & Serialization | 🟡 80% | varint ✓, FASTA ✓, archive ⚠️ |
| Phase 3: Core Algorithms | 🟡 70% | kmer_extract ✓, splitters ✓, segment ✓, lz_diff ⏸, compress ⏸ |
| Phase 4: Integration | ⏸ 0% | Not started |
| Phase 5: Decompression | ⏸ 0% | Not started |

**All completed components validated as bit-identical with C++ implementation.**

---

## Strategy: Test Harness Validation

We will rewrite AGC in Rust using a **bottom-up, component-by-component approach** with validation at each step.

### Core Principle
For each component, we:
1. Create a **C++ test harness** that exercises the function/module with known inputs
2. Write the **Rust equivalent**
3. **Verify identical outputs** (byte-for-byte or value-for-value)
4. Move to the next component

This avoids FFI complexity while ensuring correctness at every step.

---

## Phase 1: Foundation & Pure Functions ✅ 100% COMPLETE

### 1.1 Core Types & Constants ✅ COMPLETE
**Location**: `agc-common/src/types.rs`

- ✅ Version constants (AGC_VER_MAJOR, AGC_VER_MINOR, etc.)
- ✅ Type aliases (Contig = Vec<u8>, PackedBlock = Vec<u8>)
- ✅ Basic enums (file versions, orientations, Base enum)

**C++ Reference**: `src/common/defs.h`

**Validation**: ✓ Verified - Simple equality checks

---

### 1.2 K-mer Operations ✅ COMPLETE
**Location**: `agc-core/src/kmer.rs`

**Functions implemented**:
- ✅ 2-bit packed k-mer representation (MSB-first)
- ✅ `Kmer::insert()` - Insert base and maintain direct/RC forms
- ✅ `reverse_complement()` - RC of k-mer
- ✅ Canonical k-mer mode (minimum of kmer and RC)
- ✅ Full k-mer operations (reset, is_full, data access)

**C++ Reference**: `src/core/kmer.h`

**Test Harness**: ✓ Created `src/test-harness/test_kmer.cpp`

**Validation**: ✓ PASSED
```bash
./test-harness/test_kmer > rust-agc/test-data/kmer_reference.txt
cargo run --bin test-kmer > rust-agc/test-data/kmer_rust.txt
diff rust-agc/test-data/kmer_{reference,rust}.txt  # 0 differences
```

---

### 1.3 Hash Functions ✅ COMPLETE
**Location**: `agc-common/src/hash.rs`

**Functions implemented**:
- ✅ `MurMur64Hash::hash(key: u64) -> u64` - MurmurHash64
- ✅ `MurMurPair64Hash::hash(k1: u64, k2: u64) -> u64` - Hash pair
- ✅ `MurMurMix32Hash::hash(key: u32) -> u32` - MurmurHash32
- ✅ `MurMurStringsHash::hash(s: &str) -> u64` - Hash string

**C++ Reference**: `src/common/utils.h` (MurMur64Hash, MurMurPair64Hash, MurMurStringsHash)

**Test Harness**: ✓ Created `src/test-harness/test_hash.cpp`

**Validation**: ✓ PASSED - All hash variants produce identical output

---

## Phase 2: I/O & Serialization 🟡 80% COMPLETE

### 2.1 FASTA Parser ✅ COMPLETE
**Location**: `agc-core/src/genome_io.rs`

**Functions implemented**:
- ✅ `GenomeIO::new(reader)` - Create parser from BufReader
- ✅ `read_contig()` - Read raw contig (preserves original bases)
- ✅ `read_contig_converted()` - Convert IUPAC codes to 0-4 encoding
- ✅ Handle multi-line sequences, comments, whitespace
- ✅ Support all IUPAC nucleotide codes (N, R, Y, S, W, K, M, B, D, H, V, U)

**C++ Reference**: `src/core/genome_io.{cpp,h}`

**Test Harness**: ✓ Created `src/test-harness/test_fasta_parser.cpp`
- Uses real `genome_io.cpp` from AGC codebase
- Tests with multiple FASTA formats

**Validation**: ✓ PASSED
```bash
diff rust-agc/test-data/fasta_{cpp,rust}.txt  # 0 differences
```

---

### 2.2 Archive Format I/O
**Location**: `agc-common/src/archive.rs`

**Core Archive operations**:
- [ ] `Archive::create()` - Create new archive
- [ ] `Archive::open(path)` - Open existing archive
- [ ] `Archive::register_stream(name) -> stream_id` - Register stream
- [ ] `Archive::add_part(stream_id, data, metadata)` - Add data part
- [ ] `Archive::get_part(stream_id, part_id) -> (data, metadata)` - Read part
- [ ] Variable-length integer encoding/decoding
- [ ] Serialize/deserialize stream metadata

**C++ Reference**: `src/common/archive.{cpp,h}`

**Test Harness**: Create small AGC files with known content
```cpp
// test-harness/test_archive_write.cpp
// Creates minimal archive with known streams and data
```

**Validation**:
- Byte-for-byte comparison of written archives
- Read back and verify all data matches

---

### 2.3 Variable-Length Encoding ✅ COMPLETE
**Location**: `agc-common/src/varint.rs`

**Functions implemented**:
- ✅ `write_varint(writer, value: u64)` - Variable-length write
- ✅ `read_varint(reader) -> (u64, usize)` - Variable-length read
- ✅ Encoding: first byte = length, followed by big-endian value
- ✅ Special case: value=0 encoded as single 0x00 byte

**C++ Reference**: `src/common/archive.h` (write/read template methods)

**Test Harness**: ✓ Created `src/test-harness/test_varint.cpp`
- Tests: 0, 1, 127, 128, 255, 256, 65535, 65536, u64::MAX

**Validation**: ✓ PASSED
```bash
diff rust-agc/test-data/varint_{cpp,rust}.txt  # 0 differences
```

---

## Phase 3: Core Algorithms 🟡 70% COMPLETE

### 3.1 K-mer Extraction & Enumeration ✅ COMPLETE
**Location**: `agc-core/src/kmer_extract.rs`

**Functions implemented**:
- ✅ `enumerate_kmers(contig, k) -> Vec<u64>` - Extract canonical k-mers from contig
- ✅ `remove_non_singletons(vec, virtual_begin)` - Filter to singleton k-mers only
- ✅ Handle special cases (N bases reset k-mer state, short sequences)
- ✅ Efficient in-place singleton filtering

**C++ Reference**: `src/core/agc_compressor.cpp` (enumerate_kmers, remove_non_singletons)

**Test Harness**: ✓ Created `src/test-harness/test_kmer_extract.cpp`

**Validation**: ✓ PASSED
```bash
diff rust-agc/test-data/kmer_extract_{cpp,rust}.txt  # 0 differences
```

---

### 3.2 Splitter Identification ✅ COMPLETE
**Location**: `agc-core/src/splitters.rs`

**Functions implemented**:
- ✅ `determine_splitters(contigs, k) -> HashSet<u64>` - Find singleton k-mers across all contigs
- ✅ Aggregates k-mers from all reference contigs
- ✅ Sorts and filters to singleton k-mers only
- ✅ Returns as HashSet for O(1) lookup during segmentation

**C++ Reference**: `src/core/agc_compressor.cpp` (determine_splitters)

**Test Harness**: ✓ Created `src/test-harness/test_splitters.cpp`

**Validation**: ✓ PASSED
```bash
diff rust-agc/test-data/splitters_{cpp,rust}.txt  # 0 differences
```

---

### 3.3 Segmentation ✅ COMPLETE
**Location**: `agc-core/src/segment.rs`

**Functions implemented**:
- ✅ `split_at_splitters(contig, splitters, k) -> Vec<Segment>` - Split contig at splitter positions
- ✅ `Segment` struct with data, front_kmer, back_kmer fields
- ✅ Handle edge cases:
  - No splitters (entire contig as single segment)
  - Multiple consecutive splitters
  - Splitters at start/end
  - Short contigs (< k bases)
  - N bases (reset k-mer state)

**C++ Reference**: `src/common/segment.{cpp,h}`

**Test Harness**: ✓ Created `src/test-harness/test_segment.cpp`
- 7 comprehensive test cases covering all edge cases

**Validation**: ✓ PASSED
```bash
diff rust-agc/test-data/segment_{cpp,rust}.txt  # 0 differences
```

---

### 3.4 LZ Diff Encoding ⏸ NOT STARTED
**Location**: `agc-core/src/lz_diff.rs`

**Priority**: 🔴 HIGH - Required for compression pipeline

**Functions to implement**:
- ⏸ `lz_diff_encode(target: &[u8], reference: &[u8]) -> Vec<u8>` - LZ diff encoding
- ⏸ `lz_diff_decode(encoded: &[u8], reference: &[u8]) -> Vec<u8>` - Decode back to original
- ⏸ Match exact encoding scheme (literal/match opcodes, offset encoding)

**C++ Reference**: `src/common/lz_diff.{cpp,h}`

**Test Harness**: Will need to create `test_lz_diff.cpp` with encode/decode pairs

**Estimated Time**: 3-5 days (complex algorithm, needs careful bit-level validation)

---

### 3.5 Segment Compression ⏸ NOT STARTED
**Location**: `agc-core/src/segment_compression.rs`

**Priority**: 🔴 HIGH - Required for compression pipeline

**Functions to implement**:
- ⏸ `compress_segment(segment: &[u8]) -> Vec<u8>` - ZSTD compression
- ⏸ `decompress_segment(compressed: &[u8]) -> Vec<u8>` - ZSTD decompression
- ⏸ Match compression parameters (level, dictionary usage)

**C++ Reference**: Uses ZSTD library directly

**Dependencies**: `zstd` crate (already in Cargo.toml)

**Test Harness**: Simple compress/decompress roundtrip tests

**Estimated Time**: 1 day (straightforward library usage)

**Note**: Compressed output might vary slightly between ZSTD versions, but decompressed output must be byte-identical

---

## Next Steps: Path to End-to-End Compression

To reach a working end-to-end compressor, complete these components in order:

### Step 1: LZ Diff Encoding (3-5 days)
- Most complex remaining algorithm
- Critical for segment matching between samples
- Needs careful bit-level validation

### Step 2: Segment Compression (1 day)
- Straightforward ZSTD wrapper
- Test with compression/decompression roundtrips

### Step 3: Archive I/O Completion (3-4 days)
- Complete stream registration and data writing
- Implement metadata serialization
- Test with minimal AGC files

### Step 4: Collection Metadata (2-3 days)
- Sample and contig name registries
- Version-specific serialization (v1 vs v3)
- Command line history tracking

### Step 5: Compressor Core (5-7 days)
- Integrate all components into compression pipeline
- Segment grouping and matching logic
- Reference vs sample handling

### Step 6: CLI Integration (2-3 days)
- Command-line argument parsing
- `create` command implementation
- Error handling and user feedback

**Total estimated time to working compressor**: 3-4 weeks

---

## Phase 4: Integration & Compression Pipeline ⏸ NOT STARTED

### 4.1 Compressor Core
**Location**: `agc-core/src/compressor.rs`

**Components**:
- [ ] `Compressor` struct
- [ ] `compress_sample(sample_fasta) -> compressed_segments` - Single sample
- [ ] Segment grouping logic
- [ ] Metadata tracking (sample names, contig names)

**C++ Reference**: `src/core/agc_compressor.{cpp,h}`

**Validation**: Compare compressed segments for toy_ex samples

---

### 4.2 Collection Metadata
**Location**: `agc-common/src/collection.rs`

**Functions**:
- [ ] Serialize/deserialize collection metadata (v1, v3 formats)
- [ ] Sample name registry
- [ ] Contig name registry
- [ ] Command line history

**C++ Reference**: `src/common/collection_v1.{cpp,h}`, `collection_v3.{cpp,h}`

**Test Harness**: Create minimal collection, serialize, compare bytes

**Validation**: Serialized format must be byte-identical

---

### 4.3 End-to-End Compression
**Location**: `agc-cli/src/commands/create.rs`

**Full pipeline**:
- [ ] Parse command-line arguments
- [ ] Read FASTA files
- [ ] Determine splitters from reference
- [ ] Compress all samples
- [ ] Write archive
- [ ] Store metadata

**Validation**:
```bash
# C++ version
../bin/agc create toy_ex/ref.fa toy_ex/a.fa toy_ex/b.fa > cpp.agc

# Rust version
cargo run --release --bin agc create toy_ex/ref.fa toy_ex/a.fa toy_ex/b.fa > rust.agc

# Compare
cmp cpp.agc rust.agc && echo "SUCCESS: Bit-for-bit identical!"
```

---

## Phase 5: Decompression & Other Commands ⏸ NOT STARTED

### 5.1 Decompressor
**Location**: `agc-core/src/decompressor.rs`

- [ ] Read archive metadata
- [ ] Extract samples (getset)
- [ ] Extract contigs (getctg)
- [ ] Extract collection (getcol)

**Validation**: Decompressed FASTA must match original inputs byte-for-byte

---

### 5.2 Additional Commands
- [ ] `append` - Add samples to existing archive
- [ ] `listset` - List sample names
- [ ] `listctg` - List contig names
- [ ] `info` - Show archive statistics

**Validation**: Output must match C++ version exactly

---

## Development Workflow

### For Each Component:

1. **Read C++ source** - Understand the algorithm
2. **Create C++ test harness** (if needed) - Small program that exercises the function
3. **Write Rust implementation** - Idiomatic Rust
4. **Create Rust test binary** - Matching the harness
5. **Validate outputs match** - Diff or byte comparison
6. **Write unit tests** - Property-based tests with proptest
7. **Document** - Inline docs explaining algorithms
8. **Commit** - One component at a time

### Test Harness Location
```
src/test-harness/
├── test_kmer.cpp
├── test_hash.cpp
├── test_fasta_parser.cpp
├── test_archive_io.cpp
├── test_splitters.cpp
└── Makefile
```

### Rust Test Binaries
```
rust-agc/
├── agc-core/
│   └── tests/
│       ├── test_kmer.rs
│       ├── test_hash.rs
│       └── ...
```

---

## Success Metrics

- [ ] **Phase 1-3**: All component tests pass (100% output match)
- [ ] **Phase 4**: Toy example compresses to bit-identical archive
- [ ] **Phase 5**: Toy example decompresses to identical FASTA
- [ ] **Final**: Human genome dataset compresses to identical archive

---

## Timeline Estimate

- **Phase 1**: ~1-2 weeks (foundation)
- **Phase 2**: ~1-2 weeks (I/O)
- **Phase 3**: ~2-3 weeks (core algorithms)
- **Phase 4**: ~2-3 weeks (integration)
- **Phase 5**: ~1-2 weeks (decompression)

**Total**: ~7-12 weeks for feature parity with validated correctness

---

## Notes

- Use `--release` builds for performance testing
- Keep C++ test harnesses minimal and focused
- Document any deviations from C++ implementation
- Use git commits to track progress on each component
- Parallel work possible: I/O while algorithms are being developed
