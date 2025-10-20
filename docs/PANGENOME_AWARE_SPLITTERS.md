# Pangenome-Aware Splitter Improvements for AGC

## Current Algorithm Limitations

### Issue 1: Reference-Only Singleton Detection
**Problem**: Singletons are computed from the first genome only.

**Example**:
```
Genome 1 (reference): AAACCCGGGTTT → k-mer "CCC" appears once (singleton ✓)
Genome 2:             AAACCCCCCGGG → k-mer "CCC" appears 5 times
Genome 3:             AAAGGGTTTTAA → k-mer "CCC" appears 0 times
```

When we use "CCC" as a splitter:
- Genome 1: Splits once (good!)
- Genome 2: Splits 5 times in a row (terrible segmentation!)
- Genome 3: Never splits (huge segments!)

**Impact**: Splitters that work well for reference may work poorly for other genomes.

### Issue 2: No Pangenome Conservation Check
**Problem**: We don't check if splitter k-mers are conserved across samples.

**Better Approach**: Use k-mers that appear in MOST or ALL genomes, ideally at similar positions.

### Issue 3: Distance Guarantee Only in Reference
**Problem**: Splitters are spaced >= segment_size in reference, but might cluster in other genomes.

**Example**:
```
Reference: ----S1--------S2--------S3----  (nice spacing)
Sample 2:  --------S1S2----------S3------  (S1 and S2 adjacent!)
Sample 3:  ------------------S1S2S3------  (all clustered!)
```

## Proposed Improvements

### Improvement 1: Multi-Genome Splitter Discovery (Pangenome-Aware)

**Idea**: Determine splitters by analyzing a SAMPLE of genomes, not just the first one.

**Algorithm**:
```rust
/// Pangenome-aware splitter determination
///
/// Uses first N genomes (or all if N < threshold) to find CONSERVED splitters
pub fn determine_splitters_pangenome_aware(
    all_samples: &[(String, Vec<Contig>)],
    k: usize,
    segment_size: usize,
    num_reference_samples: usize,
) -> HashSet<u64> {
    let num_refs = num_reference_samples.min(all_samples.len()).max(1);

    // Step 1: Find candidate k-mers from EACH reference genome
    let mut per_genome_candidates: Vec<HashSet<u64>> = Vec::new();

    for i in 0..num_refs {
        let candidates = find_singleton_kmers(&all_samples[i].1, k);
        per_genome_candidates.push(candidates);
    }

    // Step 2: Find INTERSECTION or high-frequency k-mers
    // Option A: Strict - must appear as singleton in ALL refs
    let conserved_candidates = find_intersection(&per_genome_candidates);

    // Option B: Relaxed - must appear as singleton in >= 50% of refs
    // let conserved_candidates = find_high_frequency(&per_genome_candidates, 0.5);

    // Step 3: Find actually-used splitters across ALL reference genomes
    let mut all_splitters = HashSet::new();

    for i in 0..num_refs {
        let genome_splitters = find_actual_splitters(
            &all_samples[i].1,
            &conserved_candidates,
            k,
            segment_size
        );
        all_splitters.extend(genome_splitters);
    }

    all_splitters
}
```

**Advantages**:
- Splitters work well across multiple genomes
- Better segmentation quality for all samples
- More robust to structural variation

**Disadvantages**:
- Slower (must analyze multiple genomes)
- Fewer splitters (intersection is smaller)
- May need to use more genomes as "reference set"

### Improvement 2: Minimizer-Based Splitters (Guaranteed Spacing)

**Idea**: Use minimizers or syncmers instead of singleton k-mers.

**Why Better for Pangenomes**:
- Minimizers are **sequence-intrinsic**: Same sequence → same minimizers
- Automatically spaced (approximately every w bases for window size w)
- Work consistently across all genomes at conserved regions

**Algorithm**:
```rust
/// Use minimizers as splitters (guarantees spacing in ALL genomes)
pub fn determine_splitters_minimizer_based(
    reference_contigs: &[Contig],
    k: usize,
    window_size: usize,
) -> HashSet<u64> {
    let mut splitters = HashSet::new();

    for contig in reference_contigs {
        let minimizers = compute_minimizers(contig, k, window_size);
        splitters.extend(minimizers);
    }

    splitters
}

fn compute_minimizers(sequence: &[u8], k: usize, w: usize) -> Vec<u64> {
    let mut minimizers = Vec::new();
    let mut last_minimizer_pos = 0;

    for window_start in 0..=sequence.len().saturating_sub(w + k - 1) {
        let window_end = window_start + w + k - 1;
        let window = &sequence[window_start..window_end];

        // Find minimum k-mer in this window
        let (min_kmer, min_pos) = find_minimum_kmer(window, k);

        let abs_pos = window_start + min_pos;
        if abs_pos > last_minimizer_pos {
            minimizers.push(min_kmer);
            last_minimizer_pos = abs_pos;
        }
    }

    minimizers
}
```

**Advantages**:
- Guaranteed approximate spacing in ALL genomes
- Fast (single pass)
- Well-studied in bioinformatics (used in minimap2, etc.)

**Disadvantages**:
- Different splitter set than C++ AGC (compatibility issue)
- May have too many or too few splitters depending on window size
- Doesn't account for genome-specific features

### Improvement 3: Adaptive Splitters (Per-Sample Refinement)

**Idea**: Start with reference splitters, but allow each sample to ADD new splitters where needed.

**Algorithm**:
```rust
/// Adaptive splitter determination
/// - Use reference splitters as base
/// - Add sample-specific splitters where segments are too large
pub fn determine_splitters_adaptive(
    reference_contigs: &[Contig],
    sample_contig: &Contig,
    base_splitters: &HashSet<u64>,
    k: usize,
    segment_size: usize,
    max_segment_size: usize,
) -> HashSet<u64> {
    let mut adaptive_splitters = base_splitters.clone();

    // Segment using base splitters
    let segments = split_at_splitters(sample_contig, base_splitters, k);

    // Find segments that are too large
    for segment in segments {
        if segment.len() > max_segment_size {
            // Add new splitters from THIS segment's singleton k-mers
            let segment_singletons = find_singleton_kmers(&[segment.clone()], k);
            let new_splitters = find_actual_splitters_in_contig(
                &segment,
                &segment_singletons,
                k,
                segment_size
            );
            adaptive_splitters.extend(new_splitters);
        }
    }

    adaptive_splitters
}
```

**Advantages**:
- Prevents pathologically large segments
- Backward compatible with reference-based approach
- Each sample gets good segmentation

**Disadvantages**:
- Different samples may have different splitter sets (complexity!)
- Larger archive metadata (must store per-sample splitters)
- C++ AGC doesn't support this

### Improvement 4: Core-Genome Splitters (MSA-Aware)

**Idea**: If you have a multiple sequence alignment (MSA) or know conserved regions, use splitters from those regions ONLY.

**Algorithm**:
```rust
/// Use k-mers from CONSERVED regions as splitters
///
/// Assumes you have alignment or can identify conserved regions
pub fn determine_splitters_core_genome(
    conserved_regions: &[(usize, usize)], // (start, end) positions
    reference_contigs: &[Contig],
    k: usize,
    segment_size: usize,
) -> HashSet<u64> {
    let mut candidates = HashSet::new();

    // Only extract k-mers from conserved regions
    for (start, end) in conserved_regions {
        for contig in reference_contigs {
            if *end <= contig.len() {
                let region = &contig[*start..*end];
                let region_kmers = enumerate_kmers(region, k);
                candidates.extend(region_kmers);
            }
        }
    }

    // Find which candidates work as splitters
    find_actual_splitters(reference_contigs, &candidates, k, segment_size)
}
```

**Advantages**:
- Splitters are in conserved regions → appear in all genomes
- Best compression for pangenomes with known core genome

**Disadvantages**:
- Requires pre-computed alignment or conservation info
- Doesn't work for highly variable genomes
- Complex to implement

## Comparative Analysis

| Approach | Speed | Compression Quality | Pangenome-Aware | C++ AGC Compatible |
|----------|-------|--------------------|-----------------|--------------------|
| **Current (reference-only)** | Fast | Poor for diverse pangenomes | ❌ No | ✅ Yes |
| **Multi-genome intersection** | Medium | Better | ✅ Yes | ❌ No |
| **Minimizer-based** | Very Fast | Good | ✅ Yes | ❌ No |
| **Adaptive per-sample** | Slow | Best | ✅ Yes | ❌ No |
| **Core-genome** | Fast | Best (if core exists) | ✅ Yes | ❌ No |

## Recommendation: Hybrid Approach

**Best of both worlds**:

1. **For compatibility mode**: Keep current algorithm (first-genome singletons)
2. **Add new mode**: Multi-genome splitter discovery
   - Use first 5-10 genomes as "reference set"
   - Find k-mers that appear as singletons in >= 70% of reference set
   - Ensures splitters work across multiple genomes

**Implementation**:
```rust
pub struct SplitterConfig {
    pub mode: SplitterMode,
    pub num_reference_samples: usize,
    pub conservation_threshold: f64,
}

pub enum SplitterMode {
    /// C++ AGC compatible: first genome only
    Reference,

    /// Pangenome-aware: multiple reference genomes
    MultiGenome { conservation_threshold: f64 },

    /// Minimizer-based: guaranteed spacing
    Minimizer { window_size: usize },

    /// Adaptive: add splitters where needed
    Adaptive { max_segment_size: usize },
}
```

## Testing Strategy

To validate improvements:

1. **Compression ratio**: Measure archive size for yeast235 dataset
2. **Segmentation quality**: Measure segment size distribution across ALL samples (not just reference)
3. **Decompression speed**: Ensure no performance regression
4. **Cross-genome consistency**: Measure how many segments align across genomes

Would you like me to implement any of these approaches?
