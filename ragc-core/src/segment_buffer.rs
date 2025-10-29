use crate::compressor_streaming::SegmentGroupKey;

/// Represents a segment that has been compressed but not yet written to the archive.
/// Matches C++ AGC's buffered_seg_part structure.
#[derive(Debug)]
pub enum BufferedSegment {
    /// Segment for a NEW group that doesn't exist yet
    New {
        key: SegmentGroupKey,
        data: Vec<u8>,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
        is_rev_comp: bool,
    },
    /// Segment for a KNOWN group that already exists
    Known {
        group_id: u32,
        data: Vec<u8>,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
        is_rev_comp: bool,
    },
}

/// Buffers segments during contig processing, to be flushed at synchronization points.
/// This matches C++ AGC's CBufferedSegmentsPart behavior.
#[derive(Debug, Default)]
pub struct SegmentBuffer {
    new_segments: Vec<BufferedSegment>,
    known_segments: Vec<BufferedSegment>,
}

impl SegmentBuffer {
    pub fn new() -> Self {
        Self {
            new_segments: Vec::new(),
            known_segments: Vec::new(),
        }
    }

    /// Add a segment for a NEW group (group doesn't exist yet)
    pub fn add_new(
        &mut self,
        key: SegmentGroupKey,
        data: Vec<u8>,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
        is_rev_comp: bool,
    ) {
        self.new_segments.push(BufferedSegment::New {
            key,
            data,
            sample_name,
            contig_name,
            seg_part_no,
            is_rev_comp,
        });
    }

    /// Add a segment for a KNOWN group (group already exists)
    pub fn add_known(
        &mut self,
        group_id: u32,
        data: Vec<u8>,
        sample_name: String,
        contig_name: String,
        seg_part_no: u32,
        is_rev_comp: bool,
    ) {
        self.known_segments.push(BufferedSegment::Known {
            group_id,
            data,
            sample_name,
            contig_name,
            seg_part_no,
            is_rev_comp,
        });
    }

    /// Drain all NEW segments
    pub fn drain_new(&mut self) -> impl Iterator<Item = BufferedSegment> + '_ {
        self.new_segments.drain(..)
    }

    /// Drain all KNOWN segments
    pub fn drain_known(&mut self) -> impl Iterator<Item = BufferedSegment> + '_ {
        self.known_segments.drain(..)
    }

    /// Check if buffer is empty
    pub fn is_empty(&self) -> bool {
        self.new_segments.is_empty() && self.known_segments.is_empty()
    }

    /// Get counts for reporting
    pub fn counts(&self) -> (usize, usize) {
        (self.new_segments.len(), self.known_segments.len())
    }

    /// Clear all buffered segments
    pub fn clear(&mut self) {
        self.new_segments.clear();
        self.known_segments.clear();
    }
}
