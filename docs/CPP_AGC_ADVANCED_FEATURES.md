# C++ AGC Advanced Features - RAGC is Missing These!

## MAJOR DISCOVERY: C++ AGC Has Three Splitter Mechanisms

### Mechanism 1: Base Reference Splitters (Default Mode)
**What RAGC Currently Implements**

- Extract k-mers from FIRST genome only
- Find singletons (k-mers appearing once)
- Two-pass algorithm:
  - Pass 1: Enumerate all k-mers, sort, keep singletons → candidates
  - Pass 2: Scan reference again to find which candidates are ACTUALLY used with distance >= segment_size
- Result: Fixed splitter set for all samples

**Limitation**: Splitters work well for reference but may be poor for divergent samples.

---

### Mechanism 2: Adaptive Splitters (Optional `-a` Flag)
**RAGC IS MISSING THIS!**

**CLI Flag**: `agc create -a ...`

**When**: Enabled with `-a` (adaptive mode) flag

**How It Works**:
```cpp
// agc_compressor.cpp:2033-2039
if (adaptive_compression &&
    contig_processing_stage == contig_processing_stage_t::all_contigs &&
    split_kmer == CKmer(kmer_length, kmer_mode_t::canonical))
{
    if(contig.size() >= segment_size)
        find_new_splitters(contig, thread_id);

    return false;
}
```

**Algorithm** (`find_new_splitters`, lines 2046-2077):
```cpp
void CAGCCompressor::find_new_splitters(contig_t& ctg, uint32_t thread_id)
{
    vector<uint64_t> v_contig_kmers;
    vector<uint64_t> v_tmp;

    // 1. Extract all k-mers from this contig
    enumerate_kmers(ctg, v_contig_kmers);
    sort(v_contig_kmers.begin(), v_contig_kmers.end());
    remove_non_singletons(v_contig_kmers, 0);

    v_tmp.resize(v_contig_kmers.size());

    // 2. Exclude k-mers in reference genome - singletons
    auto p_end = set_difference(v_contig_kmers.begin(), v_contig_kmers.end(),
        v_candidate_kmers.begin() + v_candidate_kmers_offset, v_candidate_kmers.end(),
        v_tmp.begin());

    v_tmp.erase(p_end, v_tmp.end());

    // 3. Exclude k-mers in reference genome - duplicated
    p_end = set_difference(v_tmp.begin(), v_tmp.end(),
        v_duplicated_kmers.begin(), v_duplicated_kmers.end(),
        v_contig_kmers.begin());

    v_contig_kmers.erase(p_end, v_contig_kmers.end());

    // 4. Add NEW splitters to thread-local set
    auto& th_v_splitters = vv_splitters[thread_id];
    th_v_splitters.insert(th_v_splitters.end(), v_contig_kmers.begin(), v_contig_kmers.end());
}
```

**Synchronization** (lines 2186-2191):
```cpp
// After each sample, sync workers and add new splitters
pq_contigs_desc->EmplaceManyNoCost(make_tuple(
    adaptive_compression ? contig_processing_stage_t::new_splitters : contig_processing_stage_t::registration,
    "", "", contig_t()), sample_priority, no_workers);
```

**Worker Sync** (lines 1154-1177):
```cpp
if (get<0>(task) == contig_processing_stage_t::new_splitters)
{
    bar.arrive_and_wait();  // Synchronization barrier

    auto bloom_insert = [&] {
        // Add new splitters to global set
        for (auto& v : vv_splitters)
        {
            for (auto& x : v)
            {
                hs_splitters.insert_fast(x);
            }
            v.clear();
        }
    };

    bloom_insert();  // All workers add their new splitters
    bar.arrive_and_wait();  // Sync again before continuing
    continue;
}
```

**Flow**:
1. Process sample N with current splitters
2. For each contig that can't be segmented (no terminal splitter):
   - Find singleton k-mers in that contig
   - Exclude k-mers already in reference (both singletons and duplicates)
   - Add remaining k-mers as NEW splitters (sample-specific)
3. After sample N completes, synchronization barrier
4. All workers merge their new splitters into global set
5. Process sample N+1 with EXPANDED splitter set

**Advantages**:
- ✅ Pangenome-aware: Each sample contributes its own unique k-mers
- ✅ Prevents huge unsegmentable contigs
- ✅ Better compression for diverse datasets
- ✅ Gradually expanding splitter set

**Disadvantages**:
- ⚠️ Order-dependent: Later samples benefit from earlier samples' splitters
- ⚠️ Slightly slower (more synchronization)
- ⚠️ Larger splitter set (more metadata)

---

### Mechanism 3: Fallback Minimizers
**RAGC IS MISSING THIS TOO!**

**When**: Segment doesn't have terminal splitter k-mer (novel sequence region)

**Purpose**: Even without exact splitter matches, use minimizers for LZ-encoding

**How It Works**:

**1. During Splitter Finding** (lines 788-802):
```cpp
// Track minimizer k-mers between splitters
if (fallback_filter(kmer.data()) && kmer.data_dir() != kmer.data_rc())
    fallback_kmers_in_segment.emplace_back(kmer.data(), kmer.is_dir_oriented());

if (current_len >= segment_size)
{
    uint64_t d = kmer.data();

    if (binary_search(v_begin, v_end, d))  // Found splitter
    {
        v_splitters.emplace_back(d);

        // Save fallback mappings: (prev_splitter, current_splitter) → [minimizers]
        for (auto& x : fallback_kmers_in_segment)
            v_fallbacks.emplace_back(array<uint64_t, 4>{
                prev_splitter, d, x.first, (uint64_t)x.second
            });

        fallback_kmers_in_segment.clear();
        prev_splitter = d;
    }
}
```

**2. Fallback Mapping Storage** (lines 383-422):
```cpp
// Map: minimizer_kmer → [(splitter1, splitter2), ...]
unordered_map<uint64_t, vector<pair<uint64_t, uint64_t>>> map_fallback_minimizers;

void CAGCCompressor::add_fallback_mapping(
    uint64_t splitter1, uint64_t splitter2,
    vector<pair<uint64_t, bool>>& cand_fallback_kmers)
{
    pair<uint64_t, uint64_t> sp_pair_dir{ splitter1, splitter2 };
    pair<uint64_t, uint64_t> sp_pair_rc{ splitter2, splitter1 };

    for (auto x : cand_fallback_kmers)
    {
        auto& mfm_kd = map_fallback_minimizers[x.first];
        auto& to_add = x.second ? sp_pair_dir : sp_pair_rc;

        // Add (splitter1, splitter2) to this minimizer's candidate list
        if (find(mfm_kd.begin(), mfm_kd.end(), to_add) == mfm_kd.end())
            mfm_kd.emplace_back(to_add);
    }
}
```

**3. Using Fallbacks During Matching** (lines 1715-1724):
```cpp
if (!kmer_front.is_full() && !kmer_back.is_full())
{
    // No terminal splitters present - try fallback minimizers
    if (fallback_filter)
    {
        tie(pk, store_rc) = find_cand_segment_using_fallback_minimizers(segment, 1);

        if(pk != pk_empty && store_rc)
            reverse_complement_copy(segment, segment_rc);
    }
    else
        pk = pk_empty;
}
```

**Algorithm**:
- Between each pair of splitters, track ALL minimizer k-mers
- Store mapping: `minimizer_kmer → list of (splitter_before, splitter_after) pairs`
- When segment has no terminal splitters:
  - Extract minimizers from segment
  - Look up which splitter pairs contain those minimizers
  - Use best match for LZ-encoding (even without exact splitter match!)

**Advantages**:
- ✅ Can compress novel sequences without exact splitter matches
- ✅ Uses minimizers (well-studied, robust)
- ✅ Maintains some compression even for highly divergent regions

**Disadvantages**:
- ⚠️ More complex metadata (minimizer → splitter pair mappings)
- ⚠️ Slower matching (must search minimizer space)
- ⚠️ Less precise than exact splitter matches

---

## Comparison: What RAGC is Missing

| Feature | C++ AGC | RAGC | Impact |
|---------|---------|------|--------|
| **Base splitters (reference)** | ✅ Yes | ✅ Yes | Both have this |
| **Adaptive splitters (-a flag)** | ✅ Yes | ❌ NO | **RAGC can't handle divergent samples well** |
| **Fallback minimizers** | ✅ Yes | ❌ NO | **RAGC fails on novel sequences** |
| **Order dependency** | Yes (in `-a` mode) | No | Adaptive mode depends on sample order |
| **Pangenome-aware** | Yes (in `-a` mode) | No | **RAGC is reference-centric only** |

---

## Testing: Does `-a` Mode Improve Compression?

Let's test yeast235 with and without adaptive mode:

```bash
# Default mode (reference-only splitters)
time ~/agc/bin/agc create -c -v 1 -o yeast235_default.agc ~/scrapy/yeast235.fa.gz

# Adaptive mode (expanding splitter set)
time ~/agc/bin/agc create -a -c -v 1 -o yeast235_adaptive.agc ~/scrapy/yeast235.fa.gz

# Compare sizes
ls -lh yeast235_default.agc yeast235_adaptive.agc
```

**Hypothesis**:
- Adaptive mode should produce SMALLER archive (better compression)
- Because yeast strains have divergent regions that benefit from sample-specific splitters
- Trade-off: Slightly larger metadata (more splitters), but much better segment compression

---

## Recommendations for RAGC

### Priority 1: Implement Adaptive Splitters (High Impact)

**Why**: This is THE feature for pangenome compression.

**Implementation**:
```rust
pub struct StreamingCompressorConfig {
    // ... existing fields ...
    pub adaptive_mode: bool,
}

impl StreamingCompressor {
    /// Find new splitters for contigs that can't be segmented
    fn find_new_splitters(&mut self, contig: &Contig) -> HashSet<u64> {
        // 1. Extract singleton k-mers from this contig
        let contig_singletons = find_singleton_kmers(&[contig.clone()], self.config.kmer_length);

        // 2. Exclude k-mers from reference genome
        let mut new_splitters = HashSet::new();
        for kmer in contig_singletons {
            if !self.reference_kmers.contains(&kmer) {
                new_splitters.insert(kmer);
            }
        }

        new_splitters
    }

    /// Process sample with adaptive splitters
    pub fn add_sample_adaptive(&mut self, sample_name: &str, contigs: Vec<(String, Contig)>) -> Result<()> {
        let mut sample_new_splitters = HashSet::new();

        for (contig_name, contig) in contigs {
            self.collection.register_sample_contig(sample_name, &contig_name)?;

            // Try to segment with current splitters
            let segments = split_at_splitters_with_size(&contig, &self.splitters,
                self.config.kmer_length, self.config.segment_size);

            // Check if we got poor segmentation (huge segments)
            let max_segment_size = self.config.segment_size * 10;
            let needs_new_splitters = segments.iter().any(|s| s.data.len() > max_segment_size);

            if needs_new_splitters {
                // Find new splitters for this contig
                let new_splitters = self.find_new_splitters(&contig);
                sample_new_splitters.extend(new_splitters);
            }

            // Add segments (will use updated splitters for next contig)
            for (seg_idx, segment) in segments.iter().enumerate() {
                self.add_segment_with_kmers(sample_name, &contig_name, seg_idx,
                    segment.data.clone(), segment.front_kmer, segment.back_kmer, false)?;
            }
        }

        // After sample completes, add new splitters to global set
        if self.config.adaptive_mode && !sample_new_splitters.is_empty() {
            eprintln!("Adaptive mode: Adding {} new splitters from sample {}",
                sample_new_splitters.len(), sample_name);
            self.splitters.extend(sample_new_splitters);
        }

        Ok(())
    }
}
```

### Priority 2: Implement Fallback Minimizers (Medium Impact)

**Why**: Handles novel sequences gracefully.

**Complexity**: Higher - requires minimizer extraction, mapping storage, and matching algorithm.

**Recommendation**: Implement after adaptive splitters are working.

### Priority 3: Make Adaptive Mode Default (Low Priority)

**Question**: Should `-a` be default (like C++ AGC defaults to false)?

**Arguments for default=true**:
- Better compression for pangenomes
- More robust to sample order
- Main use case is diverse genomes

**Arguments for default=false**:
- Matches C++ AGC behavior
- Faster (no synchronization overhead)
- Simpler for small datasets

---

## Next Steps

1. **Test C++ AGC adaptive mode on yeast235**:
   ```bash
   time ~/agc/bin/agc create -a -c -v 1 -o /tmp/yeast235_adaptive.agc ~/scrapy/yeast235.fa.gz
   ls -lh /tmp/yeast235_adaptive.agc /tmp/yeast235_default.agc
   ```

2. **Implement adaptive mode in RAGC**:
   - Add `adaptive_mode` flag to config
   - Implement `find_new_splitters()`
   - Add synchronization between samples
   - Test on yeast235

3. **Benchmark compression improvement**:
   - Compare archive sizes
   - Measure compression ratio improvement
   - Analyze splitter set growth

4. **Document algorithm in RAGC**:
   - Explain three mechanisms
   - Add CLI flag `-a` / `--adaptive`
   - Update README with pangenome features
