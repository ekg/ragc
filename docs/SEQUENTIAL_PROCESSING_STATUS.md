# Sequential Processing Implementation - Current Status

**Date**: 2025-10-30
**Status**: In Progress (60% complete)

---

## What's Been Done

### ✅ Completed Tasks

1. **Root Cause Identified**
   - Duplicate `seg_part_no` values from Phase 1 buffering + Phase 4 splits
   - Documented in `BUG_INVESTIGATION_CHECKLIST.md`

2. **Design Document Created**
   - `SEQUENTIAL_PROCESSING_DESIGN.md` with complete implementation plan
   - Matches C++ AGC architecture

3. **Writer Thread Simplified**
   - Lines 1238-1290 updated
   - Removed broken renumbering logic
   - Direct seg_part_no usage (will be correct from sequential processing)

4. **Data Structures Updated**
   - Lines 1207-1223: Changed to use HashMap instead of buffer
   - Removed `segment_buffer: Vec<BufferedSegment>`
   - Added sequential processing structures

---

## What Needs To Be Done

### 🚧 Task: Replace Phase 2/3/4 with Sequential Loop

**Location**: Lines 1292-1940 in `ragc-core/src/compressor_streaming.rs`

**Current Code** (BROKEN - uses buffering):
```rust
// Line 1292-1940: Old phase-based processing
// Phase 2: Sort all segments
// Phase 3: Process NEW segments
// Phase 4: Process KNOWN segments (splits cause duplicate part_no here!)
```

**Need To Replace With**:
```rust
// Sequential contig processing loop
let (pack_tx_clone, next_group_id_clone) = (pack_tx.clone(), next_group_id.clone());

while let Some((sample_name, contig_name, sequence)) = iterator.next_contig()? {
    if sequence.is_empty() {
        continue;
    }

    // Log sample changes
    if current_sample != sample_name {
        if !current_sample.is_empty() && self.config.verbosity > 0 {
            println!("Finished sample: {}", current_sample);
        }
        current_sample = sample_name.clone();
        if self.config.verbosity > 0 {
            println!("Processing sample: {}", sample_name);
        }
    }

    total_contigs += 1;
    total_bases += sequence.len();

    // Segment contig
    let segments = split_at_splitters_with_size(
        &sequence,
        splitters,
        self.config.kmer_length as usize,
        self.config.segment_size as usize,
    );

    // Process each segment IMMEDIATELY (sequential!)
    let mut seg_part_no = 0;

    for segment in segments {
        let kmer_front = segment.front_kmer;
        let kmer_back = segment.back_kmer;
        let segment_data = segment.data;

        // Create normalized key
        let (key, _) = SegmentGroupKey::new_normalized(kmer_front, kmer_back);

        // Classify: NEW or KNOWN?
        let is_new = !known_groups.contains_key(&key);

        if is_new {
            // NEW SEGMENT: Create group
            let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
            known_groups.insert(key.clone(), gid);

            // Add terminators
            if kmer_front != MISSING_KMER && kmer_back != MISSING_KMER {
                group_terminators.entry(kmer_front).or_default().push(kmer_back);
                if kmer_front != kmer_back {
                    group_terminators.entry(kmer_back).or_default().push(kmer_front);
                }
            }

            // Create group and add segment
            let stream_id = gid as usize;
            let ref_stream_id = if gid >= 16 { Some(10000 + gid as usize) } else { None };
            let mut group_writer = GroupWriter::new(gid, stream_id, ref_stream_id);

            let seg_info = Self::prepare_segment_info(
                &self.config,
                &sample_name,
                &contig_name,
                seg_part_no,
                segment_data,
                Some(kmer_front),
                Some(kmer_back),
            )?;

            if let Some(pack) = group_writer.add_segment(seg_info.segment, &self.config)? {
                pack_tx_clone.send(convert_to_compressed_pack(pack)?)?;
            }

            groups.insert(key, (gid, group_writer));

            seg_part_no += 1;  // Increment by 1 for normal segment
            total_segments += 1;

        } else {
            // KNOWN SEGMENT: Try to split
            let group_id = *known_groups.get(&key).unwrap();

            // Check if we can split using existing split logic
            let can_split = kmer_front != MISSING_KMER
                         && kmer_back != MISSING_KMER
                         && segment_data.len() >= 2 * self.config.segment_size as usize;

            let mut did_split = false;

            if can_split {
                if let Some(middle_kmer) = self.find_shared_kmer(kmer_front, kmer_back, &group_terminators) {
                    // Both target groups exist, calculate split position
                    if let Some(split_pos) = Self::calculate_split_position(
                        kmer_front,
                        kmer_back,
                        middle_kmer,
                        &segment_data,
                        &groups,
                        &self.config,
                    ) {
                        // SPLIT THE SEGMENT!
                        let kmer_len = self.config.kmer_length as usize;
                        let seg2_start_pos = split_pos.saturating_sub(kmer_len / 2);

                        let seg1_data = segment_data[..seg2_start_pos + kmer_len].to_vec();
                        let seg2_data = segment_data[seg2_start_pos..].to_vec();

                        // Create first split segment (seg_part_no)
                        let seg1_key = SegmentGroupKey::new_normalized(kmer_front, middle_kmer).0;
                        let seg1_info = Self::prepare_segment_info(
                            &self.config,
                            &sample_name,
                            &contig_name,
                            seg_part_no,
                            seg1_data,
                            Some(kmer_front),
                            Some(middle_kmer),
                        )?;

                        // Create second split segment (seg_part_no + 1)
                        let seg2_key = SegmentGroupKey::new_normalized(middle_kmer, kmer_back).0;
                        let seg2_info = Self::prepare_segment_info(
                            &self.config,
                            &sample_name,
                            &contig_name,
                            seg_part_no + 1,
                            seg2_data,
                            Some(middle_kmer),
                            Some(kmer_back),
                        )?;

                        // Add both segments to their respective groups
                        for (seg_key, seg_info_prepared) in [(seg1_key, seg1_info), (seg2_key, seg2_info)] {
                            let mut group_entry = groups.entry(seg_key.clone()).or_insert_with(|| {
                                let gid = next_group_id.fetch_add(1, Ordering::SeqCst);
                                known_groups.insert(seg_key.clone(), gid);
                                let stream_id = gid as usize;
                                let ref_stream_id = if gid >= 16 { Some(10000 + gid as usize) } else { None };
                                (gid, GroupWriter::new(gid, stream_id, ref_stream_id))
                            });

                            if let Some(pack) = group_entry.1.add_segment(seg_info_prepared.segment, &self.config)? {
                                pack_tx_clone.send(convert_to_compressed_pack(pack)?)?;
                            }
                        }

                        seg_part_no += 2;  // Increment by 2 for split
                        total_segments += 2;
                        did_split = true;
                    }
                }
            }

            if !did_split {
                // No split - add as normal KNOWN segment
                let seg_info = Self::prepare_segment_info(
                    &self.config,
                    &sample_name,
                    &contig_name,
                    seg_part_no,
                    segment_data,
                    Some(kmer_front),
                    Some(kmer_back),
                )?;

                let mut group_entry = groups.get_mut(&key).unwrap();
                if let Some(pack) = group_entry.1.add_segment(seg_info.segment, &self.config)? {
                    pack_tx_clone.send(convert_to_compressed_pack(pack)?)?;
                }

                seg_part_no += 1;  // Increment by 1
                total_segments += 1;
            }
        }
    }
}

// Flush remaining groups
for (_, (_, mut group_writer)) in groups {
    if let Some(pack) = group_writer.prepare_pack(&self.config)? {
        pack_tx.send(convert_to_compressed_pack(PackToWrite::Uncompressed(pack))?)?;
    }
}

drop(pack_tx);  // Signal writer thread to finish

// Wait for writer
writer_handle.join().unwrap()?;

// Restore archive and collection
self.archive = Arc::try_unwrap(archive)?.into_inner().unwrap();
self.collection = Arc::try_unwrap(collection)?.into_inner().unwrap();

if self.config.verbosity > 0 {
    println!("Compression complete: {} contigs, {} bases, {} segments",
             total_contigs, total_bases, total_segments);
}

Ok(())
```

---

## Key Points

1. **seg_part_no increments correctly**:
   - Normal segment: `seg_part_no += 1`
   - Split segment: `seg_part_no += 2`
   - NO PRE-ASSIGNED VALUES = NO CONFLICTS!

2. **No buffering**:
   - Segments processed immediately
   - Lower memory usage

3. **Matches C++ AGC**:
   - Sequential contig processing
   - Inline NEW/KNOWN classification
   - Inline split handling

---

## Testing Plan

After implementation:
1. Build: `cargo build --release`
2. Test 2 samples: `bash /tmp/test_two_samples.sh`
3. Verify AAB#0 extracts to 12,299,861 bytes (not 8.5MB)
4. Check for duplicate seg_part_no (should be NONE)
5. Test full 10-sample dataset

---

## Estimated Time

- Code replacement: 30 minutes
- Build + fix compile errors: 15 minutes
- Testing + debugging: 30 minutes
- **Total: ~75 minutes**

The implementation is straightforward - just replacing buffered/phased processing with a sequential loop.
