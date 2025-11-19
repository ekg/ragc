# C++ AGC Fork Plan for Rust Integration

## Goal
Enable progressive replacement of C++ AGC components with Rust implementations while maintaining byte-identical output.

## Current Status
- ✅ Rust FFI for splitter detection complete and verified (matches C++ AGC exactly)
- ✅ C++ wrapper class (RustSplitterDetector) ready
- ❌ Cannot access C++ AGC private members from external class

## Fork Strategy

### Option 1: Minimal Patch Approach (RECOMMENDED)
Create minimal modifications to C++ AGC source:

**Files to modify:**
1. `/home/erik/agc/src/core/agc_compressor.h`
   - Change `private:` → `protected:` for:
     - `v_candidate_kmers`
     - `v_duplicated_kmers`
     - `v_candidate_kmers_offset`
     - `concatenated_genomes`, `adaptive_compression`
     - `fallback_frac`, `fallback_filter`
     - `out_archive`
   - Add new overload:
     ```cpp
     bool CreateWithPrecomputedSplitters(
         const string& _file_name,
         const vector<uint64_t>& singletons,
         const vector<uint64_t>& duplicates,
         ... // other parameters
     );
     ```

2. `/home/erik/agc/src/core/agc_compressor.cpp`
   - Implement CreateWithPrecomputedSplitters() that:
     - Skips determine_splitters() call
     - Sets v_candidate_kmers = singletons
     - Sets v_duplicated_kmers = duplicates
     - Continues with rest of Create() logic

**Benefits:**
- Minimal changes to C++ AGC
- Easy to maintain/update
- Clear separation between original and modified code

**Implementation:**
1. Create patch file documenting exact changes
2. Apply to local C++ AGC copy
3. Test byte-identical archives
4. Document divergence for future C++ AGC updates

### Option 2: Copy and Modify
Copy entire C++ AGC source into ragc repo and modify directly.

**Benefits:**
- Complete control
- No dependency on external C++ AGC

**Drawbacks:**
- Large code duplication
- Harder to track upstream C++ AGC changes
- More code to maintain

## Recommended Next Steps

1. **Create patch for C++ AGC** (1 hour)
   - Modify agc_compressor.h to make members protected
   - Add CreateWithPrecomputedSplitters() method
   - Test compilation

2. **Update FFI wrapper** (30 min)
   - Modify agc_compress.cpp to call CreateWithPrecomputedSplitters()
   - Pass Rust-computed splitters

3. **Test byte-identical output** (30 min)
   - Create archive with Rust splitters
   - Compare with native C++ AGC archive
   - Verify SHA256 match

4. **Document and commit** (15 min)
   - Record exact patch applied
   - Commit with clear message about fork point

## Testing Protocol

```bash
# 1. Native C++ AGC (baseline)
/home/erik/agc/bin/agc create -o baseline.agc -k 21 -s 10000 -l 20 -t 1 test.fa
sha256sum baseline.agc > baseline.sha

# 2. With modified C++ AGC using Rust splitters
./ragc-with-rust-splitters create -o rust_splitters.agc -k 21 -s 10000 -l 20 -t 1 test.fa
sha256sum rust_splitters.agc > rust.sha

# 3. MUST be identical
diff baseline.sha rust.sha || echo "ERROR: Archives differ!"
```

## Success Criteria
- ✅ Archives are byte-identical (SHA256 match)
- ✅ Both RAGC and C++ AGC can read archives
- ✅ Extraction produces identical sequences
- ✅ Patch is minimal and well-documented
