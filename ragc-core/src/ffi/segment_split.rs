// FFI helper for segment splitting
// Used when a segment needs to be split into two parts with overlap

/// Split a segment into two parts with overlap
///
/// Matches C++ AGC's segment splitting logic (agc_compressor.cpp:1434-1438):
/// ```cpp
/// uint32_t seg2_start_pos = left_size - kmer_length / 2;
/// segment2.assign(segment.begin() + seg2_start_pos, segment.end());
/// segment.resize((size_t)seg2_start_pos + kmer_length);
/// ```
///
/// # Arguments
/// * `left_size` - Size of the left part (where split occurs)
/// * `kmer_length` - Length of k-mers (for overlap calculation)
///
/// # Returns
/// Position where second segment should start
///
/// # Safety
/// This function performs safe arithmetic with overflow protection
#[no_mangle]
pub extern "C" fn ragc_calculate_seg2_start_pos(
    left_size: u32,
    kmer_length: u32,
) -> u32 {
    // Match C++ AGC logic exactly
    // seg2_start_pos = left_size - kmer_length / 2
    left_size.saturating_sub(kmer_length / 2)
}

/// Calculate the new size for the first segment after split
///
/// Matches C++ AGC: segment.resize((size_t)seg2_start_pos + kmer_length);
///
/// # Arguments
/// * `seg2_start_pos` - Starting position of second segment
/// * `kmer_length` - Length of k-mers (for overlap)
///
/// # Returns
/// New size for first segment
#[no_mangle]
pub extern "C" fn ragc_calculate_segment1_size(
    seg2_start_pos: u32,
    kmer_length: u32,
) -> u32 {
    seg2_start_pos.saturating_add(kmer_length)
}

/// Perform complete segment split calculation
///
/// Returns all values needed to split a segment into two overlapping parts.
///
/// # Arguments
/// * `left_size` - Size of the left part (where split occurs)
/// * `kmer_length` - Length of k-mers
///
/// # Returns
/// Struct containing seg2_start_pos and segment1_new_size
#[repr(C)]
pub struct SegmentSplitInfo {
    pub seg2_start_pos: u32,
    pub segment1_new_size: u32,
}

#[no_mangle]
pub extern "C" fn ragc_calculate_segment_split(
    left_size: u32,
    kmer_length: u32,
) -> SegmentSplitInfo {
    let seg2_start_pos = left_size.saturating_sub(kmer_length / 2);
    let segment1_new_size = seg2_start_pos + kmer_length;

    SegmentSplitInfo {
        seg2_start_pos,
        segment1_new_size,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seg2_start_pos() {
        // Standard case: k=21
        assert_eq!(ragc_calculate_seg2_start_pos(100, 21), 90); // 100 - 21/2 = 100 - 10 = 90

        // k=25
        assert_eq!(ragc_calculate_seg2_start_pos(200, 25), 188); // 200 - 25/2 = 200 - 12 = 188

        // Small left_size
        assert_eq!(ragc_calculate_seg2_start_pos(5, 21), 0); // 5 - 10 = saturates to 0

        // Zero left_size
        assert_eq!(ragc_calculate_seg2_start_pos(0, 21), 0);
    }

    #[test]
    fn test_segment1_size() {
        // Standard case: k=21, seg2_start=90
        assert_eq!(ragc_calculate_segment1_size(90, 21), 111); // 90 + 21 = 111

        // k=25, seg2_start=188
        assert_eq!(ragc_calculate_segment1_size(188, 25), 213); // 188 + 25 = 213

        // Zero start
        assert_eq!(ragc_calculate_segment1_size(0, 21), 21);
    }

    #[test]
    fn test_complete_split() {
        // Standard case: left_size=100, k=21
        let info = ragc_calculate_segment_split(100, 21);
        assert_eq!(info.seg2_start_pos, 90); // 100 - 10
        assert_eq!(info.segment1_new_size, 111); // 90 + 21

        // k=25
        let info = ragc_calculate_segment_split(200, 25);
        assert_eq!(info.seg2_start_pos, 188); // 200 - 12
        assert_eq!(info.segment1_new_size, 213); // 188 + 25

        // Edge case: small left_size
        let info = ragc_calculate_segment_split(5, 21);
        assert_eq!(info.seg2_start_pos, 0);
        assert_eq!(info.segment1_new_size, 21);
    }

    #[test]
    fn test_overlap_size() {
        // The overlap should be kmer_length
        // segment1 ends at seg2_start_pos + kmer_length
        // segment2 starts at seg2_start_pos
        // So overlap = kmer_length

        let info = ragc_calculate_segment_split(100, 21);
        let overlap = info.segment1_new_size - info.seg2_start_pos;
        assert_eq!(overlap, 21); // Overlap equals kmer_length

        let info = ragc_calculate_segment_split(200, 25);
        let overlap = info.segment1_new_size - info.seg2_start_pos;
        assert_eq!(overlap, 25);
    }

    #[test]
    fn test_realistic_values() {
        // Realistic scenario: segment of 500bp, split at position 250, k=21
        let info = ragc_calculate_segment_split(250, 21);
        assert_eq!(info.seg2_start_pos, 240); // 250 - 10
        assert_eq!(info.segment1_new_size, 261); // 240 + 21

        // segment1: [0..261) = 261 bases
        // segment2: [240..500) = 260 bases
        // overlap: [240..261) = 21 bases ✓
    }
}
