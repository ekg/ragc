// LZ Diff Encoding
// Encodes a target sequence as differences from a reference sequence

use ragc_common::{hash::MurMur64Hash, types::Contig};
use std::collections::HashMap;

/// Constants for LZ diff encoding
const N_CODE: u8 = 4;
const N_RUN_STARTER_CODE: u8 = 30;
const MIN_NRUN_LEN: u32 = 4;
const MAX_NO_TRIES: usize = 64;
const HASHING_STEP: usize = 4; // USE_SPARSE_HT mode

/// LZ Diff encoder/decoder (V2 implementation)
pub struct LZDiff {
    reference: Vec<u8>,
    ht: HashMap<u64, Vec<u32>>, // Hash table: kmer_hash -> list of positions
    min_match_len: u32,
    key_len: u32,
    key_mask: u64,
}

impl LZDiff {
    /// Create a new LZ diff encoder with the given minimum match length
    pub fn new(min_match_len: u32) -> Self {
        let key_len = min_match_len - (HASHING_STEP as u32) + 1;
        let key_mask = if key_len >= 32 {
            !0u64
        } else {
            (1u64 << (2 * key_len)) - 1
        };

        LZDiff {
            reference: Vec::new(),
            ht: HashMap::new(),
            min_match_len,
            key_len,
            key_mask,
        }
    }

    /// Prepare the encoder with a reference sequence
    pub fn prepare(&mut self, reference: &Contig) {
        self.reference = reference.clone();
        // Add padding for key_len
        self.reference
            .resize(self.reference.len() + self.key_len as usize, 31);
        self.build_index();
    }

    /// Build hash table index for k-mers in reference
    fn build_index(&mut self) {
        self.ht.clear();
        let ref_len = self.reference.len();

        let mut i = 0;
        while i + (self.key_len as usize) < ref_len {
            if let Some(code) = self.get_code(&self.reference[i..]) {
                let hash = MurMur64Hash::hash(code);
                // Store i / HASHING_STEP (like C++ implementation)
                self.ht
                    .entry(hash)
                    .or_default()
                    .push((i / HASHING_STEP) as u32);
            }
            i += HASHING_STEP;
        }
    }

    /// Extract k-mer code from sequence
    #[allow(clippy::needless_range_loop)]
    fn get_code(&self, seq: &[u8]) -> Option<u64> {
        let mut code = 0u64;
        for i in 0..(self.key_len as usize) {
            if seq[i] > 3 {
                return None; // Invalid base (N or other)
            }
            code = (code << 2) | (seq[i] as u64);
        }
        Some(code)
    }

    /// Extract k-mer code using sliding window optimization
    fn get_code_skip1(&self, prev_code: u64, seq: &[u8]) -> Option<u64> {
        let last_base_idx = (self.key_len as usize) - 1;
        if seq[last_base_idx] > 3 {
            return None;
        }
        let code = ((prev_code << 2) & self.key_mask) | (seq[last_base_idx] as u64);
        Some(code)
    }

    /// Check for N-run (at least 3 consecutive N bases)
    fn get_nrun_len(&self, seq: &[u8], max_len: usize) -> u32 {
        if seq.len() < 3 || seq[0] != N_CODE || seq[1] != N_CODE || seq[2] != N_CODE {
            return 0;
        }

        let mut len = 3;
        while len < max_len && seq[len] == N_CODE {
            len += 1;
        }
        len as u32
    }

    /// Encode a literal base
    fn encode_literal(&self, base: u8, encoded: &mut Vec<u8>) {
        encoded.push(b'A' + base);
    }

    /// Encode an N-run
    fn encode_nrun(&self, len: u32, encoded: &mut Vec<u8>) {
        encoded.push(N_RUN_STARTER_CODE);
        self.append_int(encoded, (len - MIN_NRUN_LEN) as i64);
        encoded.push(N_CODE);
    }

    /// Encode a match
    fn encode_match(&self, ref_pos: u32, len: Option<u32>, pred_pos: u32, encoded: &mut Vec<u8>) {
        let dif_pos = (ref_pos as i32) - (pred_pos as i32);
        self.append_int(encoded, dif_pos as i64);

        if let Some(match_len) = len {
            encoded.push(b',');
            self.append_int(encoded, (match_len - self.min_match_len) as i64);
        }

        encoded.push(b'.');
    }

    /// Append integer as ASCII decimal
    fn append_int(&self, text: &mut Vec<u8>, mut x: i64) {
        if x == 0 {
            text.push(b'0');
            return;
        }

        if x < 0 {
            text.push(b'-');
            x = -x;
        }

        let mut digits = Vec::new();
        while x > 0 {
            digits.push(b'0' + (x % 10) as u8);
            x /= 10;
        }

        // Reverse to get correct order
        digits.reverse();
        text.extend(digits);
    }

    /// Find best match in reference for the given position
    fn find_best_match(
        &self,
        hash: u64,
        target: &[u8],
        text_pos: usize,
        max_len: usize,
        no_prev_literals: usize,
    ) -> Option<(u32, u32, u32)> {
        // Returns (ref_pos, len_bck, len_fwd)

        let positions = self.ht.get(&hash)?;

        let mut best_ref_pos = 0;
        let mut best_len_bck = 0;
        let mut best_len_fwd = 0;
        let mut min_to_update = self.min_match_len as usize;

        for &pos in positions.iter().take(MAX_NO_TRIES) {
            let h_pos = (pos as usize) * HASHING_STEP;

            // Bounds check
            if h_pos >= self.reference.len() {
                continue;
            }

            let ref_ptr = &self.reference[h_pos..];
            let text_ptr = &target[text_pos..];

            // Forward match
            let f_len = Self::matching_length(text_ptr, ref_ptr, max_len);

            if f_len >= self.key_len as usize {
                // Backward match
                let mut b_len = 0;
                let max_back = no_prev_literals.min(h_pos).min(text_pos);
                while b_len < max_back {
                    if target[text_pos - b_len - 1] != self.reference[h_pos - b_len - 1] {
                        break;
                    }
                    b_len += 1;
                }

                if b_len + f_len > min_to_update {
                    best_len_bck = b_len as u32;
                    best_len_fwd = f_len as u32;
                    best_ref_pos = h_pos as u32;
                    min_to_update = b_len + f_len;
                }
            }
        }

        if (best_len_bck + best_len_fwd) as usize >= self.min_match_len as usize {
            Some((best_ref_pos, best_len_bck, best_len_fwd))
        } else {
            None
        }
    }

    /// Count matching length between two sequences
    fn matching_length(s1: &[u8], s2: &[u8], max_len: usize) -> usize {
        let mut len = 0;
        let max = max_len.min(s1.len()).min(s2.len());
        while len < max && s1[len] == s2[len] {
            len += 1;
        }
        len
    }

    /// Encode target sequence relative to reference
    pub fn encode(&mut self, target: &Contig) -> Vec<u8> {
        let mut encoded = Vec::new();

        // Optimization: if target equals reference, return empty
        if target.len() == self.reference.len() - (self.key_len as usize)
            && target
                .iter()
                .zip(self.reference.iter())
                .all(|(a, b)| a == b)
        {
            return encoded;
        }

        let text_size = target.len();
        let mut i = 0;
        let mut pred_pos = 0u32;
        let mut no_prev_literals = 0usize;
        let mut x_prev: Option<u64> = None;

        while i + (self.key_len as usize) < text_size {
            // Get k-mer code
            let x = if let Some(prev) = x_prev {
                if no_prev_literals > 0 {
                    self.get_code_skip1(prev, &target[i..])
                } else {
                    self.get_code(&target[i..])
                }
            } else {
                self.get_code(&target[i..])
            };

            x_prev = x;

            if x.is_none() {
                // Check for N-run
                let nrun_len = self.get_nrun_len(&target[i..], text_size - i);

                if nrun_len >= MIN_NRUN_LEN {
                    self.encode_nrun(nrun_len, &mut encoded);
                    i += nrun_len as usize;
                    no_prev_literals = 0;
                } else {
                    // Single literal
                    self.encode_literal(target[i], &mut encoded);
                    i += 1;
                    pred_pos += 1;
                    no_prev_literals += 1;
                }
                continue;
            }

            // Try to find match
            let hash = MurMur64Hash::hash(x.unwrap());
            let max_len = text_size - i;

            if let Some((match_pos, len_bck, len_fwd)) =
                self.find_best_match(hash, target, i, max_len, no_prev_literals)
            {
                // Handle backward extension
                if len_bck > 0 {
                    for _ in 0..len_bck {
                        encoded.pop();
                    }
                    i -= len_bck as usize;
                    pred_pos -= len_bck;
                }

                // Check if this is a match to end of sequence
                let total_len = len_bck + len_fwd;
                let len_to_encode = if i + (total_len as usize) == text_size
                    && (match_pos as usize) + (total_len as usize)
                        == self.reference.len() - (self.key_len as usize)
                {
                    None // Match to end
                } else {
                    Some(total_len)
                };

                self.encode_match(match_pos - len_bck, len_to_encode, pred_pos, &mut encoded);

                pred_pos = match_pos - len_bck + total_len;
                i += total_len as usize;
                no_prev_literals = 0;
            } else {
                // No match, encode literal
                self.encode_literal(target[i], &mut encoded);
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
            }
        }

        // Encode remaining bases as literals
        while i < text_size {
            self.encode_literal(target[i], &mut encoded);
            i += 1;
        }

        encoded
    }

    /// Decode encoded sequence using reference
    pub fn decode(&self, encoded: &[u8]) -> Vec<u8> {
        let mut decoded = Vec::new();
        let mut pred_pos = 0usize;
        let mut i = 0;

        while i < encoded.len() {
            if self.is_literal(encoded[i]) {
                let c = self.decode_literal(encoded[i]);
                let actual_c = if c == b'!' {
                    self.reference[pred_pos]
                } else {
                    c
                };
                decoded.push(actual_c);
                pred_pos += 1;
                i += 1;
            } else if encoded[i] == N_RUN_STARTER_CODE {
                let (len, consumed) = self.decode_nrun(&encoded[i..]);
                decoded.resize(decoded.len() + len as usize, N_CODE);
                i += consumed;
            } else {
                // It's a match
                let (ref_pos, len, consumed) = self.decode_match(&encoded[i..], pred_pos);
                let actual_len = if len == u32::MAX {
                    self.reference.len() - ref_pos
                } else {
                    len as usize
                };
                decoded.extend_from_slice(&self.reference[ref_pos..ref_pos + actual_len]);
                pred_pos = ref_pos + actual_len;
                i += consumed;
            }
        }

        decoded
    }

    /// Check if byte is a literal
    fn is_literal(&self, c: u8) -> bool {
        (b'A'..=b'A' + 20).contains(&c) || c == b'!'
    }

    /// Decode a literal
    fn decode_literal(&self, c: u8) -> u8 {
        if c == b'!' {
            b'!'
        } else {
            c - b'A'
        }
    }

    /// Decode an N-run, returns (length, bytes_consumed)
    fn decode_nrun(&self, data: &[u8]) -> (u32, usize) {
        let mut i = 1; // Skip starter code
        let (raw_len, len_bytes) = self.read_int(&data[i..]);
        i += len_bytes;
        i += 1; // Skip N_CODE suffix
        ((raw_len as u32) + MIN_NRUN_LEN, i)
    }

    /// Decode a match, returns (ref_pos, length, bytes_consumed)
    fn decode_match(&self, data: &[u8], pred_pos: usize) -> (usize, u32, usize) {
        let mut i = 0;
        let (raw_pos, pos_bytes) = self.read_int(&data[i..]);
        i += pos_bytes;

        let ref_pos = ((pred_pos as i64) + raw_pos) as usize;

        let len = if data[i] == b',' {
            i += 1; // Skip comma
            let (raw_len, len_bytes) = self.read_int(&data[i..]);
            i += len_bytes;
            i += 1; // Skip period
            (raw_len as u32) + self.min_match_len
        } else {
            i += 1; // Skip period
            u32::MAX // Sentinel for "to end of sequence"
        };

        (ref_pos, len, i)
    }

    /// Read ASCII decimal integer, returns (value, bytes_consumed)
    fn read_int(&self, data: &[u8]) -> (i64, usize) {
        let mut i = 0;
        let mut is_neg = false;

        if data[i] == b'-' {
            is_neg = true;
            i += 1;
        }

        let mut x = 0i64;
        while i < data.len() && data[i] >= b'0' && data[i] <= b'9' {
            x = x * 10 + ((data[i] - b'0') as i64);
            i += 1;
        }

        if is_neg {
            x = -x;
        }

        (x, i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_literal() {
        let reference = vec![0, 0, 0, 1, 1, 1];
        let target = vec![0, 1, 2, 3];

        let mut lz = LZDiff::new(18);
        lz.prepare(&reference);

        let encoded = lz.encode(&target);
        let decoded = lz.decode(&encoded);

        assert_eq!(target, decoded);
    }

    #[test]
    fn test_identical_sequences() {
        let reference = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let target = reference.clone();

        let mut lz = LZDiff::new(18);
        lz.prepare(&reference);

        let encoded = lz.encode(&target);
        // Should be empty (optimization)
        assert_eq!(encoded.len(), 0);

        // Special handling for empty encoding
        let decoded = if encoded.is_empty() && target.len() == reference.len() {
            reference.clone()
        } else {
            lz.decode(&encoded)
        };

        assert_eq!(target, decoded);
    }
}
