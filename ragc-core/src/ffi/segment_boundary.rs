// FFI helper for segment boundary calculations
// Simple arithmetic but critical for correct segment boundaries

/// Calculate split position for a segment
///
/// Matches C++ AGC's split position calculation (agc_compressor.cpp:2044):
/// ```cpp
/// split_pos = pos + 1 - kmer_length;
/// ```
///
/// This determines where the next segment starts after finding a splitter.
/// The -kmer_length creates overlap so each segment includes the full
/// k-mer at its boundaries.
///
/// # Arguments
/// * `pos` - Current position in contig (where splitter k-mer ends)
/// * `kmer_length` - Length of k-mers
///
/// # Returns
/// Position where next segment should start
#[no_mangle]
pub extern "C" fn ragc_calculate_split_position(
    pos: u64,
    kmer_length: u32,
) -> u64 {
    // Match C++ logic exactly
    // pos + 1 gives us the position after the splitter k-mer
    // Subtracting kmer_length creates the k-byte overlap
    (pos + 1).saturating_sub(kmer_length as u64)
}

/// Calculate segment end position
///
/// When a splitter is found at position `pos`, the segment ends at `pos + 1`
/// (just after the splitter k-mer completes).
///
/// Matches C++ AGC logic in get_part() call (agc_compressor.cpp:2037):
/// ```cpp
/// get_part(contig, split_pos, pos + 1 - split_pos)
/// ```
#[no_mangle]
pub extern "C" fn ragc_calculate_segment_end(
    pos: u64,
) -> u64 {
    pos + 1
}

/// Calculate segment length
///
/// Given start and end positions, calculate the segment length.
/// This matches the length parameter in get_part().
///
/// # Arguments
/// * `split_pos` - Start of segment
/// * `segment_end` - End of segment (pos + 1)
///
/// # Returns
/// Length of segment
#[no_mangle]
pub extern "C" fn ragc_calculate_segment_length(
    split_pos: u64,
    segment_end: u64,
) -> u64 {
    segment_end.saturating_sub(split_pos)
}

/// Complete boundary calculation for a segment split
///
/// Given current position and parameters, calculate all boundary values.
/// This combines the operations from compress_contig() into one call.
///
/// # Returns
/// (segment_end, new_split_pos, segment_length)
#[repr(C)]
pub struct SegmentBoundary {
    pub segment_end: u64,
    pub new_split_pos: u64,
    pub segment_length: u64,
}

#[no_mangle]
pub extern "C" fn ragc_calculate_segment_boundary(
    current_pos: u64,
    current_split_pos: u64,
    kmer_length: u32,
) -> SegmentBoundary {
    let segment_end = ragc_calculate_segment_end(current_pos);
    let segment_length = ragc_calculate_segment_length(current_split_pos, segment_end);
    let new_split_pos = ragc_calculate_split_position(current_pos, kmer_length);

    SegmentBoundary {
        segment_end,
        new_split_pos,
        segment_length,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_position() {
        // Example: pos=100, k=21
        // split_pos = 100 + 1 - 21 = 80
        assert_eq!(ragc_calculate_split_position(100, 21), 80);

        // Edge case: pos < k
        assert_eq!(ragc_calculate_split_position(10, 21), 0); // saturating_sub

        // k=3
        assert_eq!(ragc_calculate_split_position(50, 3), 48);
    }

    #[test]
    fn test_segment_end() {
        assert_eq!(ragc_calculate_segment_end(100), 101);
        assert_eq!(ragc_calculate_segment_end(0), 1);
        assert_eq!(ragc_calculate_segment_end(999), 1000);
    }

    #[test]
    fn test_segment_length() {
        // split_pos=0, segment_end=101
        assert_eq!(ragc_calculate_segment_length(0, 101), 101);

        // split_pos=80, segment_end=101
        assert_eq!(ragc_calculate_segment_length(80, 101), 21);

        // Edge case: end < start (shouldn't happen, but saturating)
        assert_eq!(ragc_calculate_segment_length(100, 50), 0);
    }

    #[test]
    fn test_complete_boundary() {
        // Scenario: Found splitter at pos=100, k=21, current split_pos=0
        let boundary = ragc_calculate_segment_boundary(100, 0, 21);

        assert_eq!(boundary.segment_end, 101);      // pos + 1
        assert_eq!(boundary.new_split_pos, 80);     // pos + 1 - k
        assert_eq!(boundary.segment_length, 101);   // segment_end - split_pos
    }

    #[test]
    fn test_complete_boundary_mid_contig() {
        // Scenario: Second segment, found splitter at pos=250, k=21, current split_pos=80
        let boundary = ragc_calculate_segment_boundary(250, 80, 21);

        assert_eq!(boundary.segment_end, 251);      // pos + 1
        assert_eq!(boundary.new_split_pos, 230);    // pos + 1 - k
        assert_eq!(boundary.segment_length, 171);   // 251 - 80
    }

    #[test]
    fn test_boundary_with_small_k() {
        // k=3 (small k-mer)
        let boundary = ragc_calculate_segment_boundary(100, 0, 3);

        assert_eq!(boundary.segment_end, 101);
        assert_eq!(boundary.new_split_pos, 98);  // pos + 1 - 3
        assert_eq!(boundary.segment_length, 101);
    }
}
