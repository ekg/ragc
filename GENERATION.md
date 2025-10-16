# How ragc Was Generated

**ragc** (Rust AGC) was created through a guided implementation process to achieve bit-for-bit compatible interaction with the AGC format.

## Implementation Approach

This implementation was developed through:

1. **Format Analysis**: Deep study of the C++ AGC reference implementation v3.2.1 to understand:
   - Binary archive structure (streams, parts, metadata)
   - Packed-contig mode with separator bytes
   - LZ diff encoding for delta compression
   - Version-specific stream naming (x{N}d for v3+)
   - Raw-only groups constraint (groups 0-15)
   - ZSTD compression with marker bytes

2. **Guided Development**: Systematic implementation guided by compatibility requirements:
   - Archive I/O matching C++ byte layout
   - Collection metadata with custom varint encoding
   - K-mer extraction and segment grouping
   - LZ diff algorithm compatible with C++
   - Proper handling of the 16 raw-only groups

3. **Iterative Testing**: Continuous validation against C++ implementation:
   - Cross-implementation roundtrip testing
   - Binary format comparison
   - Edge case identification and fixing

## Capabilities

ragc is capable of both:
- ✅ **Creating** AGC archives that C++ can read
- ✅ **Reading** AGC archives that C++ creates

This bidirectional compatibility was achieved through careful attention to:
- Exact binary format matching
- C++ convention adherence (raw-only groups, separator metadata)
- Compatible compression parameters
- Proper stream naming and structure

## Key Implementation Decisions

1. **Multi-crate architecture**: Separated into `ragc-common`, `ragc-core`, and `ragc-cli` for modularity
2. **Pure Rust**: No FFI bindings - complete rewrite for memory safety and portability
3. **Format fidelity**: Prioritized exact compatibility over optimization
4. **C++ conventions**: Respects C++ hardcoded constants like `no_raw_groups=16`

## Verification

Full compatibility verified through:
- 60+ unit tests
- Cross-implementation archive creation/extraction
- Binary format validation
- Continuous integration with C++ AGC

## Result

A production-ready, format-compatible Rust implementation that can seamlessly interoperate with the C++ reference implementation while providing Rust's safety guarantees and ecosystem benefits.
