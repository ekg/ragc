// Collection metadata structures
// Manages sample names, contig names, and segment descriptors

use crate::Archive;
use anyhow::{Context, Result};
use std::collections::HashMap;

/// A segment descriptor identifying a compressed segment
///
/// # Architecture
///
/// Segments are organized into groups by terminal k-mers (front and back).
/// Each group has two streams:
/// - **ref_stream** (xGr): Contains reference segment (in_group_id = 0)
/// - **delta_stream** (xGd): Contains LZ-encoded segments (in_group_id = 1, 2, 3, ...)
///
/// For groups >= 16 (LZ-enabled):
/// - Reference is stored raw in ref_stream
/// - Other segments are LZ-encoded against reference
///
/// For groups 0-15 (raw-only):
/// - No LZ encoding, all segments stored raw
///
/// # Critical: in_group_id Indexing
///
/// - **in_group_id = 0**: ALWAYS the reference segment (in ref_stream)
/// - **in_group_id >= 1**: Delta segments (in delta_stream)
///
/// IMPORTANT: Do NOT store reference in delta_stream! It must be:
/// 1. Written to ref_stream
/// 2. Registered with in_group_id = 0
/// 3. Not included in delta pack
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SegmentDesc {
    pub group_id: u32,
    pub in_group_id: u32, // 0 = reference, 1+ = deltas
    pub is_rev_comp: bool,
    pub raw_length: u32,
}

impl SegmentDesc {
    pub fn new(group_id: u32, in_group_id: u32, is_rev_comp: bool, raw_length: u32) -> Self {
        SegmentDesc {
            group_id,
            in_group_id,
            is_rev_comp,
            raw_length,
        }
    }

    pub fn empty() -> Self {
        SegmentDesc {
            group_id: u32::MAX,
            in_group_id: u32::MAX,
            is_rev_comp: false,
            raw_length: 0,
        }
    }
}

/// A contig with its name and segments
#[derive(Debug, Clone)]
pub struct ContigDesc {
    pub name: String,
    pub segments: Vec<SegmentDesc>,
}

impl ContigDesc {
    pub fn new(name: String) -> Self {
        ContigDesc {
            name,
            segments: Vec::new(),
        }
    }
}

/// A sample with its name and contigs
#[derive(Debug, Clone)]
pub struct SampleDesc {
    pub name: String,
    pub contigs: Vec<ContigDesc>,
}

impl SampleDesc {
    pub fn new(name: String) -> Self {
        SampleDesc {
            name,
            contigs: Vec::new(),
        }
    }
}

/// Custom variable-length integer encoding for Collection metadata
/// Different from standard varint - uses prefix bits to indicate length
///
/// Format:
/// - 0-127:         1 byte  (0xxxxxxx)
/// - 128-16511:     2 bytes (10xxxxxx xxxxxxxx)
/// - 16512-2113663: 3 bytes (110xxxxx xxxxxxxx xxxxxxxx)
/// - etc.
pub struct CollectionVarInt;

impl CollectionVarInt {
    const THR_1: u32 = 1 << 7;
    const THR_2: u32 = Self::THR_1 + (1 << 14);
    const THR_3: u32 = Self::THR_2 + (1 << 21);
    const THR_4: u32 = Self::THR_3 + (1 << 28);

    const PREF_1: u8 = 0;
    const PREF_2: u8 = 0b10000000u8;
    const PREF_3: u8 = 0b11000000u8;
    const PREF_4: u8 = 0b11100000u8;
    const PREF_5: u8 = 0b11110000u8;

    const MASK_1: u8 = 0b10000000u8;
    const MASK_2: u8 = 0b11000000u8;
    const MASK_3: u8 = 0b11100000u8;
    const MASK_4: u8 = 0b11110000u8;

    pub fn encode(data: &mut Vec<u8>, num: u32) {
        if num < Self::THR_1 {
            data.push(Self::PREF_1 + num as u8);
        } else if num < Self::THR_2 {
            let num = num - Self::THR_1;
            data.push(Self::PREF_2 + (num >> 8) as u8);
            data.push((num & 0xff) as u8);
        } else if num < Self::THR_3 {
            let num = num - Self::THR_2;
            data.push(Self::PREF_3 + (num >> 16) as u8);
            data.push(((num >> 8) & 0xff) as u8);
            data.push((num & 0xff) as u8);
        } else if num < Self::THR_4 {
            let num = num - Self::THR_3;
            data.push(Self::PREF_4 + (num >> 24) as u8);
            data.push(((num >> 16) & 0xff) as u8);
            data.push(((num >> 8) & 0xff) as u8);
            data.push((num & 0xff) as u8);
        } else {
            let num = num - Self::THR_4;
            data.push(Self::PREF_5);
            data.push(((num >> 24) & 0xff) as u8);
            data.push(((num >> 16) & 0xff) as u8);
            data.push(((num >> 8) & 0xff) as u8);
            data.push((num & 0xff) as u8);
        }
    }

    pub fn decode(ptr: &mut &[u8]) -> Result<u32> {
        if ptr.is_empty() {
            anyhow::bail!("Unexpected end of data while decoding varint");
        }

        let first = ptr[0];

        if (first & Self::MASK_1) == Self::PREF_1 {
            let num = (first - Self::PREF_1) as u32;
            *ptr = &ptr[1..];
            Ok(num)
        } else if (first & Self::MASK_2) == Self::PREF_2 {
            if ptr.len() < 2 {
                anyhow::bail!("Unexpected end of data while decoding 2-byte varint");
            }
            let num =
                ((ptr[0] as u32) << 8) + ptr[1] as u32 + Self::THR_1 - ((Self::PREF_2 as u32) << 8);
            *ptr = &ptr[2..];
            Ok(num)
        } else if (first & Self::MASK_3) == Self::PREF_3 {
            if ptr.len() < 3 {
                anyhow::bail!("Unexpected end of data while decoding 3-byte varint");
            }
            let num =
                ((ptr[0] as u32) << 16) + ((ptr[1] as u32) << 8) + ptr[2] as u32 + Self::THR_2
                    - ((Self::PREF_3 as u32) << 16);
            *ptr = &ptr[3..];
            Ok(num)
        } else if (first & Self::MASK_4) == Self::PREF_4 {
            if ptr.len() < 4 {
                anyhow::bail!("Unexpected end of data while decoding 4-byte varint");
            }
            let num = ((ptr[0] as u32) << 24)
                + ((ptr[1] as u32) << 16)
                + ((ptr[2] as u32) << 8)
                + ptr[3] as u32
                + Self::THR_3
                - ((Self::PREF_4 as u32) << 24);
            *ptr = &ptr[4..];
            Ok(num)
        } else {
            // PREF_5
            if ptr.len() < 5 {
                anyhow::bail!("Unexpected end of data while decoding 5-byte varint");
            }
            let mut num = ptr[1] as u32;
            num <<= 8;
            num += ptr[2] as u32;
            num <<= 8;
            num += ptr[3] as u32;
            num <<= 8;
            num += ptr[4] as u32;
            num += Self::THR_4;
            *ptr = &ptr[5..];
            Ok(num)
        }
    }

    pub fn encode_string(data: &mut Vec<u8>, s: &str) {
        data.extend_from_slice(s.as_bytes());
        data.push(0);
    }

    pub fn decode_string(ptr: &mut &[u8]) -> Result<String> {
        let end = ptr
            .iter()
            .position(|&b| b == 0)
            .context("Null terminator not found in string")?;

        let s = String::from_utf8(ptr[..end].to_vec()).context("Invalid UTF-8 in string")?;

        *ptr = &ptr[end + 1..];
        Ok(s)
    }
}

/// Zigzag encoding for signed differences
pub fn zigzag_encode_i64(x: i64) -> u64 {
    if x >= 0 {
        (2 * x) as u64
    } else {
        (2 * (-x) - 1) as u64
    }
}

pub fn zigzag_decode_i64(x: u64) -> i64 {
    if x & 1 != 0 {
        -(x.div_ceil(2) as i64)
    } else {
        (x / 2) as i64
    }
}

/// Zigzag encoding with prediction
pub fn zigzag_encode(x_curr: u64, x_prev: u64) -> u64 {
    if x_curr < x_prev {
        2 * (x_prev - x_curr) - 1
    } else if x_curr < 2 * x_prev {
        2 * (x_curr - x_prev)
    } else {
        x_curr
    }
}

pub fn zigzag_decode(x_val: u64, x_prev: u64) -> u64 {
    if x_val >= 2 * x_prev {
        x_val
    } else if x_val & 1 != 0 {
        (2 * x_prev - x_val) / 2
    } else {
        (x_val + 2 * x_prev) / 2
    }
}

/// Collection V3 - manages metadata for AGC archives
pub struct CollectionV3 {
    sample_desc: Vec<SampleDesc>,
    sample_ids: HashMap<String, usize>,

    // Archive stream IDs
    collection_samples_id: Option<usize>,
    collection_contigs_id: Option<usize>,
    collection_details_id: Option<usize>,

    // Configuration
    batch_size: usize,
    segment_size: u32,
    kmer_length: u32,

    // State tracking
    prev_sample_name: String,
    #[allow(dead_code)]
    placing_sample_name: String,
    #[allow(dead_code)]
    placing_sample_id: usize,
    no_samples_in_last_batch: usize,

    // For in_group_id delta encoding
    in_group_ids: Vec<i32>,
}

impl Default for CollectionV3 {
    fn default() -> Self {
        Self::new()
    }
}

impl CollectionV3 {
    pub fn new() -> Self {
        CollectionV3 {
            sample_desc: Vec::new(),
            sample_ids: HashMap::new(),
            collection_samples_id: None,
            collection_contigs_id: None,
            collection_details_id: None,
            batch_size: 1 << 20, // 1MB batches
            segment_size: 0,
            kmer_length: 0,
            prev_sample_name: String::new(),
            placing_sample_name: String::new(),
            placing_sample_id: 0,
            no_samples_in_last_batch: 0,
            in_group_ids: Vec::new(),
        }
    }

    pub fn set_config(&mut self, segment_size: u32, kmer_length: u32, batch_size: Option<usize>) {
        self.segment_size = segment_size;
        self.kmer_length = kmer_length;
        if let Some(bs) = batch_size {
            self.batch_size = bs;
        }
    }

    pub fn prepare_for_compression(&mut self, archive: &mut Archive) -> Result<()> {
        self.collection_samples_id = Some(archive.register_stream("collection-samples"));
        self.collection_contigs_id = Some(archive.register_stream("collection-contigs"));
        self.collection_details_id = Some(archive.register_stream("collection-details"));
        Ok(())
    }

    pub fn prepare_for_decompression(&mut self, archive: &Archive) -> Result<()> {
        self.collection_samples_id = archive.get_stream_id("collection-samples");
        self.collection_contigs_id = archive.get_stream_id("collection-contigs");
        self.collection_details_id = archive.get_stream_id("collection-details");

        if self.collection_samples_id.is_none() {
            anyhow::bail!("collection-samples stream not found in archive");
        }
        if self.collection_contigs_id.is_none() {
            anyhow::bail!("collection-contigs stream not found in archive");
        }
        if self.collection_details_id.is_none() {
            anyhow::bail!("collection-details stream not found in archive");
        }

        Ok(())
    }

    /// Register a sample and contig (idempotent - adds contig to existing sample if needed)
    pub fn register_sample_contig(&mut self, sample_name: &str, contig_name: &str) -> Result<bool> {
        let mut stored_sample_name = sample_name.to_string();

        if sample_name.is_empty() {
            stored_sample_name = Self::extract_contig_name(contig_name);
        }

        // Get or create sample
        let sample_id = if let Some(&id) = self.sample_ids.get(&stored_sample_name) {
            // Sample already exists
            id
        } else {
            // New sample
            let id = self.sample_ids.len();
            self.sample_ids.insert(stored_sample_name.clone(), id);
            self.sample_desc
                .push(SampleDesc::new(stored_sample_name.clone()));
            id
        };

        self.prev_sample_name = stored_sample_name;

        // Add contig to the sample (avoid duplicates)
        let sample = &mut self.sample_desc[sample_id];
        if !sample.contigs.iter().any(|c| c.name == contig_name) {
            sample
                .contigs
                .push(ContigDesc::new(contig_name.to_string()));
            Ok(true)
        } else {
            Ok(false) // Contig already registered
        }
    }

    /// Add segment placement information
    #[allow(clippy::too_many_arguments)]
    pub fn add_segment_placed(
        &mut self,
        sample_name: &str,
        contig_name: &str,
        place: usize,
        group_id: u32,
        in_group_id: u32,
        is_rev_comp: bool,
        raw_length: u32,
    ) -> Result<()> {
        let mut stored_sample_name = sample_name.to_string();

        if sample_name.is_empty() {
            stored_sample_name = Self::extract_contig_name(contig_name);
        }

        // BUGFIX: Always lookup sample_id fresh instead of caching
        // The cache was causing segments to be registered to the wrong sample
        // when GroupWriters buffer segments across sample boundaries
        let sample_id = *self
            .sample_ids
            .get(&stored_sample_name)
            .context(format!("Sample not found: {stored_sample_name}"))?;

        if std::env::var("RAGC_DEBUG_REGISTER").is_ok() {
            eprintln!("REGISTER: sample='{stored_sample_name}' (id={sample_id}), contig='{contig_name}', place={place}");
        }

        let sample = &mut self.sample_desc[sample_id];
        for contig in &mut sample.contigs {
            if contig.name == contig_name {
                if place >= contig.segments.len() {
                    contig.segments.resize(place + 1, SegmentDesc::empty());
                }
                contig.segments[place] =
                    SegmentDesc::new(group_id, in_group_id, is_rev_comp, raw_length);
                return Ok(());
            }
        }

        anyhow::bail!("Contig {contig_name} not found in sample {stored_sample_name}");
    }

    /// Get list of samples
    pub fn get_samples_list(&self, sorted: bool) -> Vec<String> {
        let mut samples: Vec<String> = self.sample_desc.iter().map(|s| s.name.clone()).collect();
        if sorted {
            samples.sort();
        }
        samples
    }

    /// Get number of samples
    pub fn get_no_samples(&self) -> usize {
        self.sample_desc.len()
    }

    /// Get number of contigs in a sample
    pub fn get_no_contigs(&self, sample_name: &str) -> Option<usize> {
        self.sample_ids
            .get(sample_name)
            .map(|&id| self.sample_desc[id].contigs.len())
    }

    /// Get contig list for a sample
    pub fn get_contig_list(&self, sample_name: &str) -> Option<Vec<String>> {
        self.sample_ids.get(sample_name).map(|&id| {
            self.sample_desc[id]
                .contigs
                .iter()
                .map(|c| c.name.clone())
                .collect()
        })
    }

    /// Get sample descriptor with all contigs and segments
    pub fn get_sample_desc(&self, sample_name: &str) -> Option<Vec<(String, Vec<SegmentDesc>)>> {
        self.sample_ids.get(sample_name).map(|&id| {
            self.sample_desc[id]
                .contigs
                .iter()
                .map(|c| (c.name.clone(), c.segments.clone()))
                .collect()
        })
    }

    /// Get contig descriptor for a specific contig in a sample
    pub fn get_contig_desc(
        &self,
        sample_name: &str,
        contig_name: &str,
    ) -> Option<Vec<SegmentDesc>> {
        self.sample_ids.get(sample_name).and_then(|&id| {
            self.sample_desc[id]
                .contigs
                .iter()
                .find(|c| c.name == contig_name)
                .map(|c| c.segments.clone())
        })
    }

    /// Get configuration parameters (segment_size, kmer_length)
    pub fn get_params(&self) -> (u32, u32) {
        (self.segment_size, self.kmer_length)
    }

    /// Extract short contig name (first word)
    fn extract_contig_name(s: &str) -> String {
        s.split_whitespace().next().unwrap_or(s).to_string()
    }

    // Serialization/Deserialization methods

    /// Serialize sample names
    fn serialize_sample_names(&self) -> Vec<u8> {
        let mut data = Vec::new();
        CollectionVarInt::encode(&mut data, self.sample_desc.len() as u32);

        for sample in &self.sample_desc {
            CollectionVarInt::encode_string(&mut data, &sample.name);
        }

        data
    }

    /// Deserialize sample names
    fn deserialize_sample_names(&mut self, data: &[u8]) -> Result<()> {
        let mut ptr = data;

        let no_samples = CollectionVarInt::decode(&mut ptr)?;

        self.sample_desc.clear();
        self.sample_ids.clear();

        for i in 0..no_samples {
            let name = CollectionVarInt::decode_string(&mut ptr)?;
            self.sample_ids.insert(name.clone(), i as usize);
            self.sample_desc.push(SampleDesc::new(name));
        }

        Ok(())
    }

    /// Split a string by spaces
    fn split_string(s: &str) -> Vec<String> {
        s.split(' ').map(|s| s.to_string()).collect()
    }

    /// Encode contig name using delta from previous name (returns raw bytes)
    fn encode_split(prev_split: &[String], curr_split: &[String]) -> Vec<u8> {
        let mut enc = Vec::new();

        for i in 0..curr_split.len() {
            if prev_split[i] == curr_split[i] {
                enc.push((-127i8) as u8); // Same component marker (0x81 = 129)
            } else if prev_split[i].len() != curr_split[i].len() {
                enc.extend_from_slice(curr_split[i].as_bytes());
            } else {
                // Encode with run-length for identical characters
                let mut cnt: i8 = 0;
                let p_bytes = prev_split[i].as_bytes();
                let c_bytes = curr_split[i].as_bytes();

                for j in 0..c_bytes.len() {
                    if p_bytes[j] == c_bytes[j] {
                        if cnt == 100 {
                            enc.push((-cnt) as u8); // Repetition marker
                            cnt = 1;
                        } else {
                            cnt += 1;
                        }
                    } else {
                        if cnt > 0 {
                            enc.push((-cnt) as u8); // Repetition marker
                            cnt = 0;
                        }
                        enc.push(c_bytes[j]);
                    }
                }

                if cnt > 0 {
                    enc.push((-cnt) as u8); // Repetition marker
                }
            }

            enc.push(b' ');
        }

        // Remove final space
        if enc.last() == Some(&b' ') {
            enc.pop();
        }

        enc
    }

    /// Decode contig name using delta from previous name (from raw bytes)
    #[allow(clippy::ptr_arg)]
    fn decode_split(prev_split: &mut Vec<String>, curr_split_bytes: &[u8]) -> String {
        // Split by space
        let parts: Vec<&[u8]> = curr_split_bytes.split(|&b| b == b' ').collect();
        let mut dec_parts = Vec::new();

        for i in 0..parts.len() {
            if parts[i].len() == 1 && (parts[i][0] as i8) == -127 {
                // Same component marker
                dec_parts.push(prev_split[i].clone());
            } else {
                let mut cmp = Vec::new();
                let p_bytes = prev_split[i].as_bytes();
                let mut p_idx = 0;

                for &byte in parts[i] {
                    let c = byte as i8;
                    if c >= 0 {
                        cmp.push(byte);
                        p_idx += 1;
                    } else {
                        // Repetition count
                        let count = -c as usize;
                        cmp.extend_from_slice(&p_bytes[p_idx..p_idx + count]);
                        p_idx += count;
                    }
                }

                let component = String::from_utf8(cmp).expect("Invalid UTF-8 in component");
                prev_split[i] = component.clone();
                dec_parts.push(component);
            }
        }

        dec_parts.join(" ")
    }

    /// Serialize contig names for a batch of samples
    fn serialize_contig_names(&self, id_from: usize, id_to: usize) -> Vec<u8> {
        let mut data = Vec::new();

        CollectionVarInt::encode(&mut data, (id_to - id_from) as u32);

        if std::env::var("RAGC_DEBUG_CONTIG_NAMES").is_ok() {
            eprintln!("SERIALIZE_CONTIG_NAMES: samples {id_from}..{id_to}");
        }

        for (sample_idx, sample) in self.sample_desc[id_from..id_to].iter().enumerate() {
            CollectionVarInt::encode(&mut data, sample.contigs.len() as u32);

            if std::env::var("RAGC_DEBUG_CONTIG_NAMES").is_ok() {
                eprintln!(
                    "  Sample {} ({}): {} contigs",
                    id_from + sample_idx,
                    sample.name,
                    sample.contigs.len()
                );
            }

            let mut prev_split = Vec::new();

            for (contig_idx, contig) in sample.contigs.iter().enumerate() {
                let curr_split = Self::split_string(&contig.name);
                let before_len = data.len();

                if curr_split.len() != prev_split.len() {
                    CollectionVarInt::encode_string(&mut data, &contig.name);

                    if std::env::var("RAGC_DEBUG_CONTIG_NAMES").is_ok() {
                        eprintln!(
                            "    Contig {}: '{}' (FULL) -> {} bytes",
                            contig_idx,
                            contig.name,
                            data.len() - before_len
                        );
                    }
                } else {
                    let enc_bytes = Self::encode_split(&prev_split, &curr_split);
                    // Write as raw bytes with null terminator
                    data.extend_from_slice(&enc_bytes);
                    data.push(0);

                    if std::env::var("RAGC_DEBUG_CONTIG_NAMES").is_ok() {
                        eprintln!("    Contig {}: '{}' (DELTA prev='{}') -> enc_bytes={:?}, {} bytes total",
                            contig_idx, contig.name, prev_split.join(" "), enc_bytes, data.len() - before_len);
                    }
                }

                prev_split = curr_split;
            }
        }

        if std::env::var("RAGC_DEBUG_CONTIG_NAMES").is_ok() {
            eprintln!("  Total serialized size: {} bytes", data.len());
        }

        data
    }

    /// Decode raw bytes string (handles non-UTF8 bytes)
    fn decode_bytes_string(ptr: &mut &[u8]) -> Result<Vec<u8>> {
        let end = ptr
            .iter()
            .position(|&b| b == 0)
            .context("Null terminator not found in bytes string")?;

        let bytes = ptr[..end].to_vec();
        *ptr = &ptr[end + 1..];
        Ok(bytes)
    }

    /// Deserialize contig names for a batch of samples
    fn deserialize_contig_names(&mut self, data: &[u8], i_sample: usize) -> Result<()> {
        let mut ptr = data;

        let no_samples_in_curr_batch = CollectionVarInt::decode(&mut ptr)? as usize;

        for i in 0..no_samples_in_curr_batch {
            let no_contigs = CollectionVarInt::decode(&mut ptr)? as usize;

            let curr_sample = &mut self.sample_desc[i_sample + i];
            curr_sample.contigs.clear();
            curr_sample.contigs.reserve(no_contigs);

            let mut prev_split = Vec::new();

            for _ in 0..no_contigs {
                let enc_bytes = Self::decode_bytes_string(&mut ptr)?;

                // Try to decode as UTF-8 first
                let name = if let Ok(enc_str) = std::str::from_utf8(&enc_bytes) {
                    let curr_split = Self::split_string(enc_str);

                    if prev_split.is_empty() || curr_split.len() != prev_split.len() {
                        // Plain text (first contig or different structure)
                        prev_split = curr_split;
                        enc_str.to_string()
                    } else {
                        // Encoded with delta
                        let decoded = Self::decode_split(&mut prev_split, &enc_bytes);
                        prev_split = Self::split_string(&decoded);
                        decoded
                    }
                } else {
                    // Contains non-UTF8 bytes, must be encoded (should not happen if prev_split is empty)
                    if prev_split.is_empty() {
                        anyhow::bail!("Cannot decode non-UTF8 bytes without previous split");
                    }
                    let decoded = Self::decode_split(&mut prev_split, &enc_bytes);
                    prev_split = Self::split_string(&decoded);
                    decoded
                };

                curr_sample.contigs.push(ContigDesc::new(name));
            }
        }

        self.no_samples_in_last_batch = no_samples_in_curr_batch;

        Ok(())
    }

    /// Clear in_group_ids tracking
    fn clear_in_group_ids(&mut self) {
        self.in_group_ids.clear();
    }

    /// Get in_group_id for a group
    fn get_in_group_id(&self, pos: usize) -> i32 {
        self.in_group_ids.get(pos).copied().unwrap_or(-1)
    }

    /// Set in_group_id for a group
    fn set_in_group_id(&mut self, pos: usize, val: i32) {
        if pos >= self.in_group_ids.len() {
            self.in_group_ids
                .resize((pos as f64 * 1.2) as usize + 1, -1);
        }
        self.in_group_ids[pos] = val;
    }

    /// Serialize contig details (5-stream format)
    fn serialize_contig_details(&mut self, id_from: usize, id_to: usize) -> [Vec<u8>; 5] {
        let mut v_data: [Vec<u8>; 5] = Default::default();

        CollectionVarInt::encode(&mut v_data[0], (id_to - id_from) as u32);

        self.clear_in_group_ids();

        let pred_raw_length = self.segment_size + self.kmer_length;

        // First pass: compute encoded values
        // IMPORTANT: Update in_group_ids IMMEDIATELY after encoding each segment (matching C++ AGC)
        let mut encoded_segments = Vec::new();

        for sample_idx in id_from..id_to {
            let mut sample_encoded = Vec::new();

            for contig_idx in 0..self.sample_desc[sample_idx].contigs.len() {
                let mut contig_encoded = Vec::new();

                for seg_idx in 0..self.sample_desc[sample_idx].contigs[contig_idx]
                    .segments
                    .len()
                {
                    let seg = &self.sample_desc[sample_idx].contigs[contig_idx].segments[seg_idx];
                    let prev_in_group_id = self.get_in_group_id(seg.group_id as usize);

                    let e_group_id = seg.group_id;
                    let e_in_group_id = if prev_in_group_id == -1 {
                        seg.in_group_id
                    } else if seg.in_group_id == 0 {
                        0
                    } else if seg.in_group_id as i32 == prev_in_group_id + 1 {
                        1
                    } else {
                        zigzag_encode(seg.in_group_id as u64, (prev_in_group_id + 1) as u64) as u32
                            + 1
                    };

                    let e_raw_length =
                        zigzag_encode(seg.raw_length as u64, pred_raw_length as u64) as u32;

                    // DEBUG: Print encoded values during serialization
                    if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                        eprintln!(
                            "ENCODE: group={}, in_group_id={}, prev={}, e_in_group_id={}",
                            seg.group_id, seg.in_group_id, prev_in_group_id, e_in_group_id
                        );
                    }

                    contig_encoded.push((e_group_id, e_in_group_id, e_raw_length, seg.is_rev_comp));

                    // Apply update IMMEDIATELY (matching C++ AGC behavior at collection_v3.cpp:581-582)
                    // NOTE: C++ AGC condition is: c_in_group_id > prev_in_group_id && c_in_group_id > 0
                    // We must match this EXACTLY even though it seems wrong!
                    if seg.in_group_id as i32 > prev_in_group_id && seg.in_group_id > 0 {
                        self.set_in_group_id(seg.group_id as usize, seg.in_group_id as i32);
                    }
                }

                sample_encoded.push(contig_encoded);
            }

            encoded_segments.push(sample_encoded);
        }

        // Second pass: write to output streams
        if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
            eprintln!(
                "SERIALIZE: Writing {} samples to collection",
                encoded_segments.len()
            );
        }

        for (sample_idx, sample_encoded) in encoded_segments.iter().enumerate() {
            let before_len = v_data[0].len();
            CollectionVarInt::encode(&mut v_data[0], sample_encoded.len() as u32);
            let after_len = v_data[0].len();

            if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                eprintln!("  Sample {}: {} contigs", sample_idx, sample_encoded.len());
                eprintln!(
                    "    Wrote bytes [{}..{}]: {:?}",
                    before_len,
                    after_len,
                    &v_data[0][before_len..after_len]
                );
            }

            for (contig_idx, contig_encoded) in sample_encoded.iter().enumerate() {
                let before_len = v_data[0].len();
                CollectionVarInt::encode(&mut v_data[0], contig_encoded.len() as u32);
                let after_len = v_data[0].len();

                if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                    eprintln!(
                        "    Contig {}: {} segments",
                        contig_idx,
                        contig_encoded.len()
                    );
                    eprintln!(
                        "      Wrote bytes [{}..{}]: {:?}",
                        before_len,
                        after_len,
                        &v_data[0][before_len..after_len]
                    );
                }

                for &(e_group_id, e_in_group_id, e_raw_length, is_rev_comp) in contig_encoded {
                    CollectionVarInt::encode(&mut v_data[1], e_group_id);
                    CollectionVarInt::encode(&mut v_data[2], e_in_group_id);
                    CollectionVarInt::encode(&mut v_data[3], e_raw_length);
                    CollectionVarInt::encode(&mut v_data[4], is_rev_comp as u32);
                }
            }

            if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                let total_segments: usize = sample_encoded.iter().map(|c| c.len()).sum();
                eprintln!("  Sample {sample_idx} total: {total_segments} segments");
            }
        }

        if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
            eprintln!(
                "Stream sizes: [0]={}, [1]={}, [2]={}, [3]={}, [4]={}",
                v_data[0].len(),
                v_data[1].len(),
                v_data[2].len(),
                v_data[3].len(),
                v_data[4].len()
            );
        }

        v_data
    }

    /// Deserialize contig details (5-stream format)
    fn deserialize_contig_details(&mut self, v_data: &[Vec<u8>; 5], i_sample: usize) -> Result<()> {
        let mut ptr0 = v_data[0].as_slice();

        let no_samples_in_curr_batch = CollectionVarInt::decode(&mut ptr0)? as usize;

        if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
            eprintln!("DESERIALIZE: Reading {no_samples_in_curr_batch} samples from collection");
        }

        // First pass: decode structure and counts
        let mut structure = Vec::new();
        for sample_idx in 0..no_samples_in_curr_batch {
            let no_contigs = CollectionVarInt::decode(&mut ptr0)? as usize;
            let mut contig_seg_counts = Vec::new();

            if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                eprintln!("  Sample {sample_idx}: {no_contigs} contigs");
            }

            for contig_idx in 0..no_contigs {
                let no_segments = CollectionVarInt::decode(&mut ptr0)? as usize;
                contig_seg_counts.push(no_segments);

                if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                    eprintln!("    Contig {contig_idx}: {no_segments} segments");
                }
            }

            structure.push(contig_seg_counts);
        }

        // Decode detail streams
        let mut v_det: [Vec<u32>; 5] = Default::default();
        let mut no_items = 0;

        for counts in &structure {
            for &count in counts {
                no_items += count;
            }
        }

        for i in 1..5 {
            v_det[i].reserve(no_items);
            let mut ptr = v_data[i].as_slice();

            for _ in 0..no_items {
                v_det[i].push(CollectionVarInt::decode(&mut ptr)?);
            }
        }

        // Reconstruct segments - first pass: decode all values
        self.clear_in_group_ids();

        let pred_raw_length = self.segment_size + self.kmer_length;
        let mut item_idx = 0;

        // Collect all decoded segments first
        // IMPORTANT: Update in_group_ids IMMEDIATELY after decoding each segment (matching C++ AGC)
        let mut decoded_segments = Vec::new();

        for contig_seg_counts in &structure {
            let mut sample_segs = Vec::new();

            for &no_segments in contig_seg_counts {
                let mut contig_segs = Vec::new();

                for _ in 0..no_segments {
                    let c_group_id = v_det[1][item_idx];
                    let prev_in_group_id = self.get_in_group_id(c_group_id as usize);

                    let e_in_group_id = v_det[2][item_idx];
                    let c_in_group_id = if prev_in_group_id == -1 {
                        e_in_group_id
                    } else if e_in_group_id == 0 {
                        0
                    } else if e_in_group_id == 1 {
                        (prev_in_group_id + 1) as u32
                    } else {
                        zigzag_decode(e_in_group_id as u64 - 1, (prev_in_group_id + 1) as u64)
                            as u32
                    };

                    // DEBUG: Print encoded vs decoded values for investigation
                    if std::env::var("RAGC_DEBUG_COLLECTION").is_ok() {
                        eprintln!("DECODE: item={item_idx}, group={c_group_id}, e_in_group_id={e_in_group_id}, prev={prev_in_group_id}, c_in_group_id={c_in_group_id}");
                    }

                    let c_raw_length =
                        zigzag_decode(v_det[3][item_idx] as u64, pred_raw_length as u64) as u32;
                    let c_is_rev_comp = v_det[4][item_idx] != 0;

                    contig_segs.push(SegmentDesc::new(
                        c_group_id,
                        c_in_group_id,
                        c_is_rev_comp,
                        c_raw_length,
                    ));

                    // Apply update IMMEDIATELY (matching C++ AGC behavior at collection_v3.cpp:674-675)
                    // NOTE: C++ AGC condition is: c_in_group_id > prev_in_group_id && c_in_group_id > 0
                    // We must match this EXACTLY even though it causes issues with multiple samples per group!
                    if c_in_group_id as i32 > prev_in_group_id && c_in_group_id > 0 {
                        self.set_in_group_id(c_group_id as usize, c_in_group_id as i32);
                    }

                    item_idx += 1;
                }

                sample_segs.push(contig_segs);
            }

            decoded_segments.push(sample_segs);
        }

        // Second pass: assign to sample_desc
        for (i, sample_segs) in decoded_segments.into_iter().enumerate() {
            let curr_sample = &mut self.sample_desc[i_sample + i];

            for (j, contig_segs) in sample_segs.into_iter().enumerate() {
                curr_sample.contigs[j].segments = contig_segs;
            }
        }

        Ok(())
    }

    /// Store sample names batch (with ZSTD compression)
    pub fn store_batch_sample_names(&mut self, archive: &mut Archive) -> Result<()> {
        let v_tmp = self.serialize_sample_names();

        let v_data = zstd::encode_all(&v_tmp[..], 19).context("Failed to compress sample names")?;

        let stream_id = self
            .collection_samples_id
            .context("collection-samples stream not registered")?;

        archive.add_part(stream_id, &v_data, v_tmp.len() as u64)?;

        Ok(())
    }

    /// Load sample names batch (with ZSTD decompression)
    pub fn load_batch_sample_names(&mut self, archive: &mut Archive) -> Result<()> {
        let stream_id = self
            .collection_samples_id
            .context("collection-samples stream not found")?;

        let (v_tmp, raw_size) = archive
            .get_part(stream_id)?
            .context("No sample names batch found")?;

        let v_data = zstd::decode_all(&v_tmp[..]).context("Failed to decompress sample names")?;

        if v_data.len() != raw_size as usize {
            anyhow::bail!(
                "Decompressed size mismatch: expected {}, got {}",
                raw_size,
                v_data.len()
            );
        }

        self.deserialize_sample_names(&v_data)?;

        Ok(())
    }

    /// Load a batch of contigs (names + details)
    #[allow(clippy::needless_range_loop)]
    pub fn load_contig_batch(&mut self, archive: &mut Archive, id_batch: usize) -> Result<()> {
        // Load contig names
        let contig_stream_id = self
            .collection_contigs_id
            .context("collection-contigs stream not found")?;

        let (v_tmp_names, raw_size_names) = archive.get_part_by_id(contig_stream_id, id_batch)?;
        let v_data_names =
            zstd::decode_all(&v_tmp_names[..]).context("Failed to decompress contig names")?;

        if v_data_names.len() != raw_size_names as usize {
            anyhow::bail!(
                "Decompressed size mismatch for contig names: expected {}, got {}",
                raw_size_names,
                v_data_names.len()
            );
        }

        self.deserialize_contig_names(&v_data_names, id_batch * self.batch_size)?;

        // Load contig details
        let details_stream_id = self
            .collection_details_id
            .context("collection-details stream not found")?;

        let (v_stream, _) = archive.get_part_by_id(details_stream_id, id_batch)?;

        // Unpack the 5 streams
        let mut ptr = v_stream.as_slice();
        let mut sizes: [(u32, u32); 5] = Default::default();

        for i in 0..5 {
            sizes[i].0 = CollectionVarInt::decode(&mut ptr)?; // raw size
            sizes[i].1 = CollectionVarInt::decode(&mut ptr)?; // compressed size
        }

        let mut v_data_details_compressed: [Vec<u8>; 5] = Default::default();
        for i in 0..5 {
            v_data_details_compressed[i] = ptr[..sizes[i].1 as usize].to_vec();
            ptr = &ptr[sizes[i].1 as usize..];
        }

        // Decompress each stream
        let mut v_data_details: [Vec<u8>; 5] = Default::default();
        for i in 0..5 {
            v_data_details[i] = zstd::decode_all(&v_data_details_compressed[i][..])
                .context(format!("Failed to decompress contig details stream {i}"))?;

            if v_data_details[i].len() != sizes[i].0 as usize {
                anyhow::bail!(
                    "Decompressed size mismatch for details stream {}: expected {}, got {}",
                    i,
                    sizes[i].0,
                    v_data_details[i].len()
                );
            }
        }

        self.deserialize_contig_details(&v_data_details, id_batch * self.batch_size)?;

        Ok(())
    }

    /// Store a batch of contigs (names + details)
    #[allow(clippy::needless_range_loop)]
    pub fn store_contig_batch(
        &mut self,
        archive: &mut Archive,
        id_from: usize,
        id_to: usize,
    ) -> Result<()> {
        // Store contig names
        let v_tmp_names = self.serialize_contig_names(id_from, id_to);
        let v_data_names =
            zstd::encode_all(&v_tmp_names[..], 18).context("Failed to compress contig names")?;

        let contig_stream_id = self
            .collection_contigs_id
            .context("collection-contigs stream not registered")?;

        archive.add_part(contig_stream_id, &v_data_names, v_tmp_names.len() as u64)?;

        // Store contig details (5 streams)
        let v_data_details_raw = self.serialize_contig_details(id_from, id_to);

        let mut v_data_details_compressed: [Vec<u8>; 5] = Default::default();
        for i in 0..5 {
            v_data_details_compressed[i] = zstd::encode_all(&v_data_details_raw[i][..], 19)
                .context(format!("Failed to compress contig details stream {i}"))?;
        }

        // Pack into single stream
        let mut v_stream = Vec::new();
        for i in 0..5 {
            CollectionVarInt::encode(&mut v_stream, v_data_details_raw[i].len() as u32);
            CollectionVarInt::encode(&mut v_stream, v_data_details_compressed[i].len() as u32);
        }
        for i in 0..5 {
            v_stream.extend_from_slice(&v_data_details_compressed[i]);
        }

        let details_stream_id = self
            .collection_details_id
            .context("collection-details stream not registered")?;

        archive.add_part(details_stream_id, &v_stream, 0)?;

        // Clear contigs from memory (for large datasets)
        for sample in &mut self.sample_desc[id_from..id_to] {
            sample.contigs.clear();
            sample.contigs.shrink_to_fit();
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collection_varint_1byte() {
        let mut data = Vec::new();
        CollectionVarInt::encode(&mut data, 0);
        CollectionVarInt::encode(&mut data, 42);
        CollectionVarInt::encode(&mut data, 127);

        let mut ptr = data.as_slice();
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 0);
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 42);
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 127);
    }

    #[test]
    fn test_collection_varint_2byte() {
        let mut data = Vec::new();
        CollectionVarInt::encode(&mut data, 128);
        CollectionVarInt::encode(&mut data, 255);
        CollectionVarInt::encode(&mut data, 16511);

        let mut ptr = data.as_slice();
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 128);
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 255);
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 16511);
    }

    #[test]
    fn test_collection_varint_3byte() {
        let mut data = Vec::new();
        CollectionVarInt::encode(&mut data, 16512);
        CollectionVarInt::encode(&mut data, 65536);
        CollectionVarInt::encode(&mut data, 2113663);

        let mut ptr = data.as_slice();
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 16512);
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 65536);
        assert_eq!(CollectionVarInt::decode(&mut ptr).unwrap(), 2113663);
    }

    #[test]
    fn test_collection_varint_string() {
        let mut data = Vec::new();
        CollectionVarInt::encode_string(&mut data, "Hello");
        CollectionVarInt::encode_string(&mut data, "World");

        let mut ptr = data.as_slice();
        assert_eq!(CollectionVarInt::decode_string(&mut ptr).unwrap(), "Hello");
        assert_eq!(CollectionVarInt::decode_string(&mut ptr).unwrap(), "World");
    }

    #[test]
    fn test_zigzag_i64() {
        assert_eq!(zigzag_encode_i64(0), 0);
        assert_eq!(zigzag_encode_i64(1), 2);
        assert_eq!(zigzag_encode_i64(-1), 1);
        assert_eq!(zigzag_encode_i64(2), 4);
        assert_eq!(zigzag_encode_i64(-2), 3);

        assert_eq!(zigzag_decode_i64(0), 0);
        assert_eq!(zigzag_decode_i64(2), 1);
        assert_eq!(zigzag_decode_i64(1), -1);
        assert_eq!(zigzag_decode_i64(4), 2);
        assert_eq!(zigzag_decode_i64(3), -2);
    }

    #[test]
    fn test_zigzag_predictive() {
        // Test predictive zigzag encoding
        assert_eq!(zigzag_encode(100, 100), 0); // Same as previous
        assert_eq!(zigzag_encode(101, 100), 2); // Slightly above
        assert_eq!(zigzag_encode(99, 100), 1); // Slightly below
        assert_eq!(zigzag_encode(200, 100), 200); // Far away

        assert_eq!(zigzag_decode(0, 100), 100);
        assert_eq!(zigzag_decode(2, 100), 101);
        assert_eq!(zigzag_decode(1, 100), 99);
        assert_eq!(zigzag_decode(200, 100), 200);
    }

    #[test]
    fn test_segment_desc() {
        let seg = SegmentDesc::new(42, 7, true, 1000);
        assert_eq!(seg.group_id, 42);
        assert_eq!(seg.in_group_id, 7);
        assert!(seg.is_rev_comp);
        assert_eq!(seg.raw_length, 1000);
    }

    #[test]
    fn test_collection_register() {
        let mut coll = CollectionV3::new();

        assert!(coll.register_sample_contig("sample1", "contig1").unwrap());
        assert!(coll.register_sample_contig("sample1", "contig2").unwrap());
        assert!(coll.register_sample_contig("sample2", "contig1").unwrap());

        assert_eq!(coll.get_no_samples(), 2);

        let samples = coll.get_samples_list(false);
        assert_eq!(samples, vec!["sample1", "sample2"]);
    }

    /// Test that in_group_id delta encoding/decoding works correctly
    /// This specifically tests the bug where deserialize updates in_group_id
    /// during the loop instead of after (like serialize does)
    #[test]
    fn test_in_group_id_delta_encoding_roundtrip() {
        let mut coll = CollectionV3::new();
        coll.set_config(60000, 31, None);

        // Register samples and contigs
        coll.register_sample_contig("sample1", "contig1").unwrap();
        coll.register_sample_contig("sample1", "contig2").unwrap();
        coll.register_sample_contig("sample2", "contig1").unwrap();

        // Add multiple segments from the same group (group 93) with consecutive in_group_ids
        // This reproduces the real-world scenario where multiple contigs have segments in the same group
        let test_segments = vec![
            ("sample1", "contig1", 0, 93, 0, false, 61000), // group 93, in_group_id 0
            ("sample1", "contig1", 1, 93, 1, false, 61000), // group 93, in_group_id 1
            ("sample1", "contig2", 0, 93, 2, false, 61000), // group 93, in_group_id 2
            ("sample1", "contig2", 1, 93, 3, false, 61000), // group 93, in_group_id 3
            ("sample2", "contig1", 0, 93, 4, false, 61000), // group 93, in_group_id 4
            ("sample2", "contig1", 1, 93, 5, false, 61000), // group 93, in_group_id 5
            ("sample2", "contig1", 2, 93, 6, false, 61000), // group 93, in_group_id 6
        ];

        for (sample, contig, place, group_id, in_group_id, is_rev_comp, raw_len) in &test_segments {
            coll.add_segment_placed(
                sample,
                contig,
                *place,
                *group_id,
                *in_group_id,
                *is_rev_comp,
                *raw_len,
            )
            .unwrap();
        }

        // Serialize the collection details
        let serialized = coll.serialize_contig_details(0, 2);

        // Create a new collection and deserialize
        let mut coll2 = CollectionV3::new();
        coll2.set_config(60000, 31, None);

        // Must set up the same sample structure
        coll2.register_sample_contig("sample1", "contig1").unwrap();
        coll2.register_sample_contig("sample1", "contig2").unwrap();
        coll2.register_sample_contig("sample2", "contig1").unwrap();

        // Deserialize
        coll2.deserialize_contig_details(&serialized, 0).unwrap();

        // Verify all in_group_ids match original values
        let sample1_contig1 = coll2.get_contig_desc("sample1", "contig1").unwrap();
        assert_eq!(
            sample1_contig1.len(),
            2,
            "sample1/contig1 should have 2 segments"
        );
        assert_eq!(
            sample1_contig1[0].in_group_id, 0,
            "sample1/contig1 segment 0 should have in_group_id=0"
        );
        assert_eq!(
            sample1_contig1[1].in_group_id, 1,
            "sample1/contig1 segment 1 should have in_group_id=1"
        );

        let sample1_contig2 = coll2.get_contig_desc("sample1", "contig2").unwrap();
        assert_eq!(
            sample1_contig2.len(),
            2,
            "sample1/contig2 should have 2 segments"
        );
        assert_eq!(
            sample1_contig2[0].in_group_id, 2,
            "sample1/contig2 segment 0 should have in_group_id=2"
        );
        assert_eq!(
            sample1_contig2[1].in_group_id, 3,
            "sample1/contig2 segment 1 should have in_group_id=3"
        );

        let sample2_contig1 = coll2.get_contig_desc("sample2", "contig1").unwrap();
        assert_eq!(
            sample2_contig1.len(),
            3,
            "sample2/contig1 should have 3 segments"
        );
        assert_eq!(
            sample2_contig1[0].in_group_id, 4,
            "sample2/contig1 segment 0 should have in_group_id=4"
        );
        assert_eq!(
            sample2_contig1[1].in_group_id, 5,
            "sample2/contig1 segment 1 should have in_group_id=5"
        );
        assert_eq!(
            sample2_contig1[2].in_group_id, 6,
            "sample2/contig1 segment 2 should have in_group_id=6"
        );
    }
}
