// Buffered segment storage matching C++ AGC's CBufferedSegPart
// Reference: agc_compressor.h lines 27-536

use std::collections::BTreeSet;
use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::{Arc, Mutex};

/// Segment part data
///
/// Matches C++ AGC's seg_part_t (agc_compressor.h:29-120) and kk_seg_part_t (lines 124-165)
///
/// Fields:
/// - `kmer1`, `kmer2`: First and last k-mers of segment
/// - `sample_name`, `contig_name`: Origin of segment
/// - `seg_data`: Compressed segment bytes
/// - `is_rev_comp`: Whether segment is reverse complemented
/// - `seg_part_no`: Part number within contig
///
/// Ordering: By (sample_name, contig_name, seg_part_no) for deterministic storage
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SegmentPart {
    pub kmer1: u64,
    pub kmer2: u64,
    pub sample_name: String,
    pub contig_name: String,
    pub seg_data: Vec<u8>,
    pub is_rev_comp: bool,
    pub seg_part_no: u32,
}

impl PartialOrd for SegmentPart {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SegmentPart {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Match C++ AGC ordering: (sample_name, contig_name, seg_part_no)
        (&self.sample_name, &self.contig_name, self.seg_part_no).cmp(&(
            &other.sample_name,
            &other.contig_name,
            other.seg_part_no,
        ))
    }
}

/// Thread-safe list of segments for one group
///
/// Matches C++ AGC's list_seg_part_t (agc_compressor.h:169-296)
struct SegmentPartList {
    /// Vector of segments (grows during compression, read during storage)
    parts: Mutex<Vec<SegmentPart>>,

    /// Virtual begin index for pop operations (avoid actually removing elements)
    virt_begin: Mutex<usize>,
}

impl SegmentPartList {
    fn new() -> Self {
        SegmentPartList {
            parts: Mutex::new(Vec::new()),
            virt_begin: Mutex::new(0),
        }
    }

    /// Add segment (thread-safe)
    ///
    /// Matches C++ AGC's list_seg_part_t::emplace (agc_compressor.h:224-228)
    fn emplace(&self, part: SegmentPart) {
        let mut parts = self.parts.lock().unwrap();
        parts.push(part);
    }

    /// Sort segments by (sample, contig, part_no)
    ///
    /// Matches C++ AGC's list_seg_part_t::sort (agc_compressor.h:235-238)
    fn sort(&self) {
        let mut parts = self.parts.lock().unwrap();
        parts.sort();
    }

    /// Pop next segment (uses virt_begin to avoid actual removal)
    ///
    /// Matches C++ AGC's list_seg_part_t::pop (agc_compressor.h:251-265)
    fn pop(&self) -> Option<SegmentPart> {
        let mut virt_begin = self.virt_begin.lock().unwrap();
        let parts = self.parts.lock().unwrap();

        if *virt_begin >= parts.len() {
            drop(parts);  // Release lock before clearing
            let mut parts = self.parts.lock().unwrap();
            *virt_begin = 0;
            parts.clear();
            return None;
        }

        let part = parts[*virt_begin].clone();
        *virt_begin += 1;

        Some(part)
    }

    /// Check if list is empty (from virt_begin perspective)
    ///
    /// Matches C++ AGC's list_seg_part_t::empty (agc_compressor.h:246-249)
    fn is_empty(&self) -> bool {
        let virt_begin = self.virt_begin.lock().unwrap();
        let parts = self.parts.lock().unwrap();
        *virt_begin >= parts.len()
    }

    /// Clear list and reset virt_begin
    ///
    /// Matches C++ AGC's list_seg_part_t::clear (agc_compressor.h:240-244)
    fn clear(&self) {
        let mut parts = self.parts.lock().unwrap();
        let mut virt_begin = self.virt_begin.lock().unwrap();
        parts.clear();
        *virt_begin = 0;
    }

    fn size(&self) -> usize {
        let parts = self.parts.lock().unwrap();
        parts.len()
    }
}

/// Buffered segments storage
///
/// Matches C++ AGC's CBufferedSegPart (agc_compressor.h:27-536)
///
/// **Architecture**:
/// - `vl_seg_part`: Vector of thread-safe lists, indexed by group_id (KNOWN segments)
/// - `s_seg_part`: BTreeSet of NEW segments (not yet assigned group_id)
///
/// **Workflow**:
/// 1. During compression: Workers call `add_known()` or `add_new()`
/// 2. At registration barrier: Main thread calls `sort_known()`, `process_new()`, `distribute_segments()`
/// 3. During storage: Workers call `get_vec_id()` and `get_part()` to read segments
/// 4. After storage: Main thread calls `clear()`
pub struct BufferedSegments {
    /// KNOWN segments indexed by group_id
    /// Matches C++ AGC's vector<list_seg_part_t> vl_seg_part (line 298)
    vl_seg_part: Vec<SegmentPartList>,

    /// NEW segments (no group_id yet)
    /// Matches C++ AGC's set<kk_seg_part_t> s_seg_part (line 300)
    s_seg_part: Mutex<BTreeSet<SegmentPart>>,

    /// Atomic counter for reading segments (starts at size-1, decrements to 0)
    /// Matches C++ AGC's atomic<int32_t> a_v_part_id (line 303)
    a_v_part_id: AtomicI32,

    /// Mutex for resizing vl_seg_part
    /// Matches C++ AGC's mutex mtx (line 301)
    resize_mtx: Mutex<()>,
}

impl BufferedSegments {
    /// Create new buffered segments storage
    ///
    /// Matches C++ AGC's CBufferedSegPart::CBufferedSegPart (agc_compressor.h:308-311)
    ///
    /// # Arguments
    /// * `no_raw_groups` - Initial number of group IDs
    pub fn new(no_raw_groups: usize) -> Self {
        let mut vl_seg_part = Vec::with_capacity(no_raw_groups);
        for _ in 0..no_raw_groups {
            vl_seg_part.push(SegmentPartList::new());
        }

        BufferedSegments {
            vl_seg_part,
            s_seg_part: Mutex::new(BTreeSet::new()),
            a_v_part_id: AtomicI32::new(0),
            resize_mtx: Mutex::new(()),
        }
    }

    /// Add segment to KNOWN group
    ///
    /// Matches C++ AGC's CBufferedSegPart::add_known (agc_compressor.h:320-324)
    ///
    /// Thread-safe: `SegmentPartList::emplace()` has internal mutex
    pub fn add_known(
        &self,
        group_id: u32,
        kmer1: u64,
        kmer2: u64,
        sample_name: String,
        contig_name: String,
        seg_data: Vec<u8>,
        is_rev_comp: bool,
        seg_part_no: u32,
    ) {
        self.vl_seg_part[group_id as usize].emplace(SegmentPart {
            kmer1,
            kmer2,
            sample_name,
            contig_name,
            seg_data,
            is_rev_comp,
            seg_part_no,
        });
    }

    /// Add NEW segment (not yet assigned group_id)
    ///
    /// Matches C++ AGC's CBufferedSegPart::add_new (agc_compressor.h:326-331)
    ///
    /// Thread-safe: Locks `s_seg_part` mutex
    pub fn add_new(
        &self,
        kmer1: u64,
        kmer2: u64,
        sample_name: String,
        contig_name: String,
        seg_data: Vec<u8>,
        is_rev_comp: bool,
        seg_part_no: u32,
    ) {
        let mut s_seg_part = self.s_seg_part.lock().unwrap();
        s_seg_part.insert(SegmentPart {
            kmer1,
            kmer2,
            sample_name,
            contig_name,
            seg_data,
            is_rev_comp,
            seg_part_no,
        });
    }

    /// Sort all KNOWN segments in parallel
    ///
    /// Matches C++ AGC's CBufferedSegPart::sort_known (agc_compressor.h:333-377)
    ///
    /// **Note**: In C++ AGC, this uses std::async for parallelism. For now, we'll use
    /// a simple sequential implementation. Parallelism can be added later with rayon.
    ///
    /// # Arguments
    /// * `_num_threads` - Number of threads (unused in sequential version)
    pub fn sort_known(&self, _num_threads: usize) {
        // TODO: Implement parallel sorting with rayon
        // For now, sequential sorting
        for list in &self.vl_seg_part {
            list.sort();
        }
    }

    /// Process NEW segments: assign group IDs and move to KNOWN
    ///
    /// Matches C++ AGC's CBufferedSegPart::process_new (agc_compressor.h:384-415)
    ///
    /// Returns: Number of NEW groups created
    pub fn process_new(&mut self) -> u32 {
        let _lock = self.resize_mtx.lock().unwrap();
        let mut s_seg_part = self.s_seg_part.lock().unwrap();

        if s_seg_part.is_empty() {
            return 0;
        }

        // Assign group IDs to unique (kmer1, kmer2) pairs
        let mut m_kmers = std::collections::HashMap::new();
        let mut group_id = self.vl_seg_part.len() as u32;

        for part in s_seg_part.iter() {
            let key = (part.kmer1, part.kmer2);
            if !m_kmers.contains_key(&key) {
                m_kmers.insert(key, group_id);
                group_id += 1;
            }
        }

        let no_new = group_id - self.vl_seg_part.len() as u32;

        // Resize vl_seg_part to accommodate new groups
        let new_size = group_id as usize;
        if self.vl_seg_part.capacity() < new_size {
            self.vl_seg_part.reserve((new_size as f64 * 1.2) as usize - self.vl_seg_part.len());
        }
        while self.vl_seg_part.len() < new_size {
            self.vl_seg_part.push(SegmentPartList::new());
        }

        // Move NEW segments to KNOWN groups
        for part in s_seg_part.iter() {
            let key = (part.kmer1, part.kmer2);
            let group_id = m_kmers[&key] as usize;

            self.vl_seg_part[group_id].emplace(part.clone());
        }

        s_seg_part.clear();

        no_new
    }

    /// Get the number of NEW segments (not yet assigned group_id)
    ///
    /// **For testing** - count segments in s_seg_part
    pub fn get_num_new(&self) -> usize {
        let s_seg_part = self.s_seg_part.lock().unwrap();
        s_seg_part.len()
    }

    /// Distribute segments from src_id to range [dest_from, dest_to)
    ///
    /// Matches C++ AGC's CBufferedSegPart::distribute_segments (agc_compressor.h:417-435)
    ///
    /// **Pattern**: Round-robin distribution
    ///
    /// Example: `distribute_segments(0, 0, num_workers)`
    /// - Distributes group 0 segments among workers 0..num_workers
    pub fn distribute_segments(&self, src_id: u32, dest_id_from: u32, dest_id_to: u32) {
        let src_id = src_id as usize;
        let no_in_src = self.vl_seg_part[src_id].size();
        let mut dest_id_curr = dest_id_from;

        for _ in 0..no_in_src {
            if dest_id_curr != src_id as u32 {
                if let Some(part) = self.vl_seg_part[src_id].pop() {
                    self.vl_seg_part[dest_id_curr as usize].emplace(part);
                }
            }

            dest_id_curr += 1;
            if dest_id_curr == dest_id_to {
                dest_id_curr = dest_id_from;
            }
        }
    }

    /// Clear all buffered segments
    ///
    /// Matches C++ AGC's CBufferedSegPart::clear (agc_compressor.h:461-507)
    ///
    /// # Arguments
    /// * `_num_threads` - Number of threads (unused in sequential version)
    pub fn clear(&mut self, _num_threads: usize) {
        // TODO: Implement parallel clearing with rayon
        // For now, sequential clearing
        let _lock = self.resize_mtx.lock().unwrap();

        let mut s_seg_part = self.s_seg_part.lock().unwrap();
        s_seg_part.clear();
        drop(s_seg_part);

        for list in &self.vl_seg_part {
            list.clear();
        }
    }

    /// Restart reading from highest group_id
    ///
    /// Matches C++ AGC's CBufferedSegPart::restart_read_vec (agc_compressor.h:509-514)
    pub fn restart_read_vec(&self) {
        let _lock = self.resize_mtx.lock().unwrap();
        self.a_v_part_id
            .store((self.vl_seg_part.len() - 1) as i32, Ordering::SeqCst);
    }

    /// Atomically get next group_id to process (decrements from size-1 to 0)
    ///
    /// Matches C++ AGC's CBufferedSegPart::get_vec_id (agc_compressor.h:516-520)
    ///
    /// Returns: group_id to process, or negative if done
    pub fn get_vec_id(&self) -> i32 {
        self.a_v_part_id.fetch_sub(1, Ordering::SeqCst)
    }

    /// Check if group is empty
    ///
    /// Matches C++ AGC's CBufferedSegPart::is_empty_part (agc_compressor.h:527-530)
    pub fn is_empty_part(&self, group_id: i32) -> bool {
        if group_id < 0 || group_id as usize >= self.vl_seg_part.len() {
            return true;
        }
        self.vl_seg_part[group_id as usize].is_empty()
    }

    /// Pop next segment from group
    ///
    /// Matches C++ AGC's CBufferedSegPart::get_part (agc_compressor.h:532-535)
    ///
    /// Returns: (kmer1, kmer2, sample_name, contig_name, seg_data, is_rev_comp, seg_part_no)
    pub fn get_part(
        &self,
        group_id: i32,
    ) -> Option<(u64, u64, String, String, Vec<u8>, bool, u32)> {
        if group_id < 0 || group_id as usize >= self.vl_seg_part.len() {
            return None;
        }

        self.vl_seg_part[group_id as usize].pop().map(|part| {
            (
                part.kmer1,
                part.kmer2,
                part.sample_name,
                part.contig_name,
                part.seg_data,
                part.is_rev_comp,
                part.seg_part_no,
            )
        })
    }

    /// Get current number of groups
    pub fn get_no_parts(&self) -> usize {
        self.vl_seg_part.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_part_ordering() {
        let part1 = SegmentPart {
            kmer1: 100,
            kmer2: 200,
            sample_name: "sample1".to_string(),
            contig_name: "chr1".to_string(),
            seg_data: vec![0, 1, 2],
            is_rev_comp: false,
            seg_part_no: 0,
        };

        let part2 = SegmentPart {
            kmer1: 100,
            kmer2: 200,
            sample_name: "sample1".to_string(),
            contig_name: "chr1".to_string(),
            seg_data: vec![3, 4, 5],
            is_rev_comp: false,
            seg_part_no: 1,
        };

        assert!(part1 < part2);  // Ordered by seg_part_no
    }

    #[test]
    fn test_buffered_segments_add_known() {
        let buf = BufferedSegments::new(10);

        buf.add_known(
            5,
            100,
            200,
            "sample1".to_string(),
            "chr1".to_string(),
            vec![0, 1, 2, 3],
            false,
            0,
        );

        assert!(!buf.is_empty_part(5));
    }

    #[test]
    fn test_buffered_segments_add_new_and_process() {
        let mut buf = BufferedSegments::new(10);

        // Add NEW segment
        buf.add_new(
            300,
            400,
            "sample1".to_string(),
            "chr1".to_string(),
            vec![4, 5, 6, 7],
            false,
            0,
        );

        // Process NEW segments
        let no_new = buf.process_new();
        assert_eq!(no_new, 1);  // One new group created

        // New group should be at index 10 (after initial 10)
        assert_eq!(buf.get_no_parts(), 11);
        assert!(!buf.is_empty_part(10));
    }

    #[test]
    fn test_buffered_segments_get_vec_id() {
        let buf = BufferedSegments::new(5);

        buf.restart_read_vec();

        // Should return 4, 3, 2, 1, 0, then negative
        assert_eq!(buf.get_vec_id(), 4);
        assert_eq!(buf.get_vec_id(), 3);
        assert_eq!(buf.get_vec_id(), 2);
        assert_eq!(buf.get_vec_id(), 1);
        assert_eq!(buf.get_vec_id(), 0);
        assert!(buf.get_vec_id() < 0);
    }

    #[test]
    fn test_buffered_segments_get_part() {
        let buf = BufferedSegments::new(10);

        buf.add_known(
            5,
            100,
            200,
            "sample1".to_string(),
            "chr1".to_string(),
            vec![0, 1, 2, 3],
            false,
            0,
        );

        let part = buf.get_part(5);
        assert!(part.is_some());

        let (kmer1, kmer2, sample, contig, data, is_rev, part_no) = part.unwrap();
        assert_eq!(kmer1, 100);
        assert_eq!(kmer2, 200);
        assert_eq!(sample, "sample1");
        assert_eq!(contig, "chr1");
        assert_eq!(data, vec![0, 1, 2, 3]);
        assert_eq!(is_rev, false);
        assert_eq!(part_no, 0);

        // Second call should return None (list empty)
        assert!(buf.get_part(5).is_none());
    }

    #[test]
    fn test_buffered_segments_sort() {
        let buf = BufferedSegments::new(1);

        // Add segments in wrong order
        buf.add_known(
            0,
            100,
            200,
            "sample1".to_string(),
            "chr1".to_string(),
            vec![2],
            false,
            2,
        );
        buf.add_known(
            0,
            100,
            200,
            "sample1".to_string(),
            "chr1".to_string(),
            vec![0],
            false,
            0,
        );
        buf.add_known(
            0,
            100,
            200,
            "sample1".to_string(),
            "chr1".to_string(),
            vec![1],
            false,
            1,
        );

        buf.sort_known(1);

        // Should pop in sorted order: 0, 1, 2
        let part0 = buf.get_part(0).unwrap();
        assert_eq!(part0.6, 0);  // seg_part_no

        let part1 = buf.get_part(0).unwrap();
        assert_eq!(part1.6, 1);

        let part2 = buf.get_part(0).unwrap();
        assert_eq!(part2.6, 2);
    }
}
