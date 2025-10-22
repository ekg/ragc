# RAGC Performance Fix: Refactor to C++ AGC Architecture

## Problem Identified
Current architecture has massive mutex contention (67% system time):
- We parallelize at CONTIG level
- Multiple workers compete for same 770 groups  
- Each segment acquisition requires mutex lock
- Result: 66s for 10 genomes (should be ~10s)

## Root Cause
C++ AGC uses DIFFERENT parallelization strategy:
- C++ AGC: Parallelize at GROUP level (each worker gets exclusive groups)
- RAGC: Parallelize at CONTIG level (workers compete for shared groups)

## Solution: Match C++ AGC Architecture

### Phase 1: Single-threaded Segmentation
```
for each contig:
    segments = split_at_splitters(contig)
    all_segments.extend(segments)
```

### Phase 2: Single-threaded Grouping  
```
groups = HashMap<SegmentGroupKey, Vec<Segment>>::new()
for segment in all_segments:
    key = (kmer_front, kmer_back)
    groups[key].push(segment)
```

### Phase 3: Parallel Group Processing
```
// Partition groups across workers
groups.par_chunks(N).for_each(|group_chunk| {
    for (key, segments) in group_chunk:
        // Process all segments - NO MUTEX needed!
        for segment in segments:
            process_segment(segment)
})
```

## Implementation Location
File: `ragc-core/src/compressor_streaming.rs`
Function: `add_contigs_parallel_streaming` (line ~600)

## Key Changes
1. Remove channel-based contig distribution
2. Collect all segments first (single-threaded)
3. Group by key (single-threaded)
4. Use rayon par_iter on groups (parallel, no contention!)
5. Remove ALL mutex locks on group writers
6. Update terminator tracking to work with new flow

## Expected Result
- 10 genomes: ~10 seconds (vs current 66s)
- System time: <5% (vs current 67%)
- 235 genomes: proportionally faster

## Status
- [x] Problem identified
- [x] Root cause found  
- [x] Solution designed
- [ ] Implementation in progress
- [ ] Testing
