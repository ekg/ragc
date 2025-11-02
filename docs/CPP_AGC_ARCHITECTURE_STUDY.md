# C++ AGC Multi-Sample Architecture - Complete Study

**Date**: 2025-11-02
**Goal**: Understand C++ AGC's worker-based parallelization to fix RAGC

## Critical Discovery

**The fundamental architecture**: C++ AGC uses **buffering + batch registration**

## Data Flow

### Phase 1: Parallel Segmentation (Workers)

```
Worker Threads (all in parallel):
├─ Pull contig from priority queue
├─ Segment contig (split_at_splitters)
├─ For each segment:
│  ├─ Determine k-mer pair (front, back)
│  ├─ Check if group exists in map_segments
│  │
│  ├─ If group EXISTS:
│  │  └─ buffered_seg_part.add_known(group_id, segment)
│  │     └─ Adds to vector at vl_seg_part[group_id]
│  │
│  └─ If group UNKNOWN:
│     └─ buffered_seg_part.add_new(kmer_pair, segment)
│        └─ Adds to set s_seg_part
│
└─ Continue until registration token received
```

### Phase 2: Registration Barrier (Thread 0 Only)

```
Thread 0 (while others wait at barrier):
├─ buffered_seg_part.sort_known(n_t)
│  └─ Sort each group's vector by (sample, contig, seg_part_no)
│
├─ buffered_seg_part.process_new()
│  ├─ For each segment in s_seg_part (unknown set):
│  │  ├─ Assign new group_id
│  │  ├─ Add to map_segments[kmer_pair] = group_id
│  │  └─ Move to vl_seg_part[group_id]
│  │
│  └─ Return count of new groups
│
├─ Register new stream IDs in archive
│
└─ buffered_seg_part.distribute_segments()
   └─ Distribute segments from set to vectors
```

### Phase 3: Parallel Compression (ALL Threads)

```
ALL Worker Threads (in parallel):
├─ Each thread calls store_segments(zstd_cctx, zstd_dctx)
│
├─ Loop: get_vec_id() returns next group_id block
│  │
│  ├─ For each group_id in block:
│  │  ├─ While segments available:
│  │  │  ├─ Pop segment from vl_seg_part[group_id]
│  │  │  │
│  │  │  ├─ If first segment in group:
│  │  │  │  ├─ Create CSegment object
│  │  │  │  ├─ Register in map_segments (if new)
│  │  │  │  └─ Update map_segments_terminators + SORT
│  │  │  │
│  │  │  ├─ Compress and add to group
│  │  │  └─ Record metadata
│  │  │
│  │  └─ Next segment
│  │
│  └─ Next block
│
└─ Barrier: wait for all threads to finish
```

## Key Data Structures

### CBufferedSegPart

```cpp
class CBufferedSegPart {
    // Known groups: indexed by group_id
    vector<list_seg_part_t> vl_seg_part;

    // Unknown groups: temporary set
    set<kk_seg_part_t> s_seg_part;

    // Thread-safe atomic for distributing work
    atomic<int32_t> a_v_part_id;

    // Add segment to known group
    void add_known(group_id, kmer1, kmer2, ...);

    // Add segment to unknown set
    void add_new(kmer1, kmer2, ...);

    // Sort known segments
    void sort_known(n_threads);

    // Process unknowns → assign groups
    uint32_t process_new();

    // Distribute segments to vectors
    void distribute_segments(...);

    // Get next group_id block (thread-safe)
    int32_t get_vec_id();
};
```

### Worker Thread Loop

```cpp
while (true) {
    task = pq_contigs_desc_working->PopLarge();

    if (task == registration_token) {
        // BARRIER 1: Wait for all workers
        bar.arrive_and_wait();

        // Thread 0 only: register segments
        if (thread_id == 0)
            register_segments(n_t);

        // BARRIER 2: Wait for registration
        bar.arrive_and_wait();

        // ALL threads: compress segments
        store_segments(zstd_cctx, zstd_dctx);

        // BARRIER 3: Wait for completion
        bar.arrive_and_wait();

        // Thread 0 only: cleanup
        if (thread_id == 0)
            buffered_seg_part.clear();

        // BARRIER 4: Final sync
        bar.arrive_and_wait();

        continue;
    }

    // Process contig (segment + add to buffer)
    compress_contig(task, ...);
}
```

## Why This Works

### 1. **Buffering Prevents Data Loss**
- Workers add segments to buffer (not processing immediately)
- Buffer holds ALL segments until registration
- NO segments are lost because buffer isn't cleared until after compression

### 2. **Batch Registration**
- Thread 0 processes ALL unknown segments at once
- All groups are registered before ANY compression starts
- map_segments is complete when store_segments runs

### 3. **Parallel Compression is Safe**
- Each thread gets different group_id blocks via atomic get_vec_id()
- No race conditions on group creation (already done in registration)
- Each group processes its own segments independently

### 4. **Inline Splitting Works**
- Workers can check map_segments during segmentation
- If match found, segment goes to known group immediately
- If no match, segment goes to unknown set for later processing

## What RAGC Got Wrong

### ❌ Immediate Processing Approach (BROKEN)

```rust
// Thread 0 registration - WRONG!
for seg in pending_segments.drain() {
    let gid = get_or_create_group(seg.key);
    compress_and_write_immediately(seg);  // ← TOO EARLY!
}

// Problem: Workers still adding to pending_segments!
// Later segments never get processed → DATA LOSS
```

### ✅ Buffered Batch Processing (CORRECT)

```rust
// Workers: Add to buffer
for seg in all_segments {
    if group_exists(seg.key) {
        buffer.add_known(group_id, seg);
    } else {
        buffer.add_new(seg.key, seg);
    }
}

// Thread 0: Register ALL at once
for seg in buffer.unknown_segments() {
    let gid = create_group(seg.key);
    buffer.move_to_known(gid, seg);
}

// ALL threads: Compress from buffer
for block in buffer.get_next_block() {
    for seg in block.segments {
        compress_and_write(seg);
    }
}
```

## Implementation Plan for RAGC

### Option 1: Match C++ AGC Exactly

1. Create `SegmentBuffer` struct:
   - `Vec<Vec<UncompressedSegment>>` for known groups
   - `Vec<UncompressedSegment>` for unknown segments
   - Atomic counter for work distribution

2. Worker thread changes:
   - Add segments to buffer (not channels)
   - Don't compress immediately

3. Registration phase:
   - Thread 0 only
   - Process ALL unknown segments
   - Assign group IDs

4. Compression phase:
   - ALL threads participate
   - Each gets group_id blocks atomically
   - Compress from buffer

### Option 2: Optimize Current RAGC

Keep existing Rayon architecture but fix the data flow:

1. **Collect ALL segments first**
   ```rust
   let all_segments: Vec<UncompressedSegment> = files
       .par_iter()
       .flat_map(|file| segment_file(file))
       .collect();  // ← Buffer everything
   ```

2. **Group segments by key**
   ```rust
   let mut groups: HashMap<SegmentGroupKey, Vec<Segment>> = HashMap::new();
   for seg in all_segments {
       groups.entry(seg.key).or_default().push(seg);
   }
   ```

3. **Compress groups in parallel**
   ```rust
   groups.par_iter().for_each(|(key, segments)| {
       let gid = get_or_create_group(key);
       for seg in segments {
           compress_and_add(gid, seg);
       }
   });
   ```

## Key Insights

1. **Never process incrementally during parallel segmentation**
   - All segmentation must complete before registration
   - Registration must complete before compression
   - Clear phase boundaries prevent data loss

2. **Registration is a barrier operation**
   - All workers must finish segmentation
   - One thread processes unknowns
   - All workers wait for completion

3. **Compression can be parallel because groups are independent**
   - Each group has its own CSegment object
   - No shared state between groups during compression
   - Work distribution via atomic counter

4. **The buffer is the key abstraction**
   - Decouples segmentation from compression
   - Allows workers to run independently
   - Provides synchronization point for registration

## Next Steps

1. ✅ Study complete - architecture understood
2. Choose implementation approach (Option 1 or Option 2)
3. Implement segment buffer
4. Test with multi-sample dataset
5. Verify 0% data loss
6. Measure performance improvement
