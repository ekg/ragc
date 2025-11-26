// LZ Diff Encoding
// Encodes a target sequence as differences from a reference sequence

#![allow(clippy::same_item_push)]

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
    reference_len: usize,       // Original length before padding
    // Legacy map used by encode/decode paths (OK if not bit-identical)
    ht: HashMap<u64, Vec<u32>>, // Hash table: kmer_hash -> list of positions (legacy)
    // Linear-probing table for exact coding-cost matching with C++
    ht_lp: Vec<u32>,            // stores (i / HASHING_STEP) or u32::MAX for empty
    ht_mask: u64,
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
            reference_len: 0,
            ht: HashMap::new(),
            ht_lp: Vec::new(),
            ht_mask: 0,
            min_match_len,
            key_len,
            key_mask,
        }
    }

    /// Prepare the encoder with a reference sequence
    pub fn prepare(&mut self, reference: &Contig) {
        self.reference = reference.clone();
        self.reference_len = reference.len(); // Store original length before padding
                                              // Add padding for key_len
        self.reference
            .resize(self.reference.len() + self.key_len as usize, 31);

        // Pre-allocate hash table capacity based on reference length
        // to avoid repeated rehashing (major source of page faults!)
        // Each entry in ht corresponds to one k-mer position (every HASHING_STEP bases)
        let expected_entries = (self.reference.len() / HASHING_STEP) + 1;
        if self.ht.capacity() < expected_entries {
            self.ht = HashMap::with_capacity(expected_entries);
        }

        self.build_index();
        self.build_index_lp();
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
                // Pre-allocate capacity to avoid repeated 0→4→8→16 growths
                self.ht
                    .entry(hash)
                    .or_insert_with(|| Vec::with_capacity(4))
                    .push((i / HASHING_STEP) as u32);
            }
            i += HASHING_STEP;
        }
    }

    /// Build linear-probing hash table (exactly like C++ CLZDiffBase::make_index32)
    fn build_index_lp(&mut self) {
        // Count valid k-mer positions as in C++ prepare_index (sparse mode)
        let mut ht_size: u64 = 0;
        let mut no_prev_valid: u32 = 0;
        let mut cnt_mod: u32 = 0;
        let key_len_mod: u32 = self.key_len % (HASHING_STEP as u32);
        for &c in &self.reference {
            if c < 4 { no_prev_valid += 1; } else { no_prev_valid = 0; }
            cnt_mod += 1;
            if cnt_mod == HASHING_STEP as u32 { cnt_mod = 0; }
            if cnt_mod == key_len_mod && no_prev_valid >= self.key_len {
                ht_size += 1;
            }
        }

        // Adjust size by load factor (0.7) and round to power of two then double
        let mut ht_size = (ht_size as f64 / 0.7) as u64;
        if ht_size == 0 { ht_size = 1; }
        while (ht_size & (ht_size - 1)) != 0 { ht_size &= ht_size - 1; }
        ht_size <<= 1;
        if ht_size < 8 { ht_size = 8; }

        self.ht_mask = ht_size - 1;
        self.ht_lp.clear();
        self.ht_lp.resize(ht_size as usize, u32::MAX);

        // Insert positions with linear probing (sparse step)
        let ref_len = self.reference.len();
        let mut i = 0usize;
        while i + (self.key_len as usize) < ref_len {
            if let Some(code) = self.get_code(&self.reference[i..]) {
                let base = (MurMur64Hash::hash(code) & self.ht_mask) as usize;
                for j in 0..MAX_NO_TRIES {
                    let idx = (base + j) & (self.ht_mask as usize);
                    if self.ht_lp[idx] == u32::MAX {
                        self.ht_lp[idx] = (i / HASHING_STEP) as u32;
                        break;
                    }
                }
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

        // Write digits directly to output (in reverse), then reverse just that portion
        let start_pos = text.len();
        while x > 0 {
            text.push(b'0' + (x % 10) as u8);
            x /= 10;
        }

        // Reverse just the digits we added
        text[start_pos..].reverse();
    }

    /// Find best match using legacy map (used by encode)
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

    /// Find best match using linear-probing table (exactly like C++ for cost vectors)
    fn find_best_match_lp(
        &self,
        hash: u64,
        target: &[u8],
        text_pos: usize,
        max_len: usize,
        no_prev_literals: usize,
    ) -> Option<(u32, u32, u32)> {
        if self.ht_lp.is_empty() { return None; }

        let mut best_ref_pos = 0u32;
        let mut best_len_bck = 0u32;
        let mut best_len_fwd = 0u32;
        let mut min_to_update = self.min_match_len as usize;

        let mut ht_pos = (hash & self.ht_mask) as usize;

        for j in 0..MAX_NO_TRIES {
            let idx = (ht_pos + j) & (self.ht_mask as usize);
            let slot = self.ht_lp[idx];
            if slot == u32::MAX { break; }

            let h_pos = (slot as usize) * HASHING_STEP;
            if h_pos >= self.reference.len() { continue; }

            let ref_ptr = &self.reference[h_pos..];
            let text_ptr = &target[text_pos..];
            let f_len = Self::matching_length(text_ptr, ref_ptr, max_len);
            if f_len >= self.key_len as usize {
                let mut b_len = 0usize;
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
        // Pre-allocate capacity to avoid repeated reallocations
        // Typical LZ compression achieves 2-4:1, so estimate capacity as target_len / 2
        let mut encoded = Vec::with_capacity(target.len() / 2);

        // Optimization: if target equals reference, return empty
        if target.len() == self.reference_len
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
                self.find_best_match_lp(hash, target, i, max_len, no_prev_literals)
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
                    && (match_pos as usize) + (total_len as usize) == self.reference_len
                {
                    None // Match to end
                } else {
                    Some(total_len)
                };

                let adjusted_match_pos = match_pos - len_bck;

                // C++ AGC optimization (lz_diff.cpp lines 769-779): when match_pos == pred_pos,
                // convert preceding literals that match the reference to '!' for better compression.
                // IMPORTANT: This must be done BEFORE encode_match, so the last bytes in buffer are literals.
                // The '!' character is decoded by looking up reference[pred_pos].
                if adjusted_match_pos == pred_pos {
                    let e_size = encoded.len();
                    // C++: for (uint32_t i = 1; i < e_size && i < match_pos; ++i)
                    let max_scan = e_size.min(adjusted_match_pos as usize);
                    for scan_i in 1..max_scan {
                        let enc_idx = e_size - scan_i;
                        let c = encoded[enc_idx];
                        // Stop if not a literal (A-Z range)
                        if c < b'A' || c > b'Z' {
                            break;
                        }
                        // Check if literal matches reference at corresponding position
                        let base = c - b'A';
                        let ref_idx = adjusted_match_pos as usize - scan_i;
                        if base == self.reference[ref_idx] {
                            encoded[enc_idx] = b'!';
                        }
                    }
                }

                self.encode_match(adjusted_match_pos, len_to_encode, pred_pos, &mut encoded);

                pred_pos = adjusted_match_pos + total_len;
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
                    // Match to end: use original reference length (before padding)
                    self.reference_len - ref_pos
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

    /// Get coding cost vector for target sequence
    /// This computes the per-position cost of encoding the target against the reference
    /// Returns a vector where v_costs[i] is the cost of encoding position i
    /// If prefix_costs=true, match cost is placed at start of match; otherwise at end
    pub fn get_coding_cost_vector(&self, target: &Contig, prefix_costs: bool) -> Vec<u32> {
        let mut v_costs = Vec::with_capacity(target.len());

        if self.reference.is_empty() {
            return v_costs;
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
                    let tc = self.coding_cost_nrun(nrun_len);
                    if prefix_costs {
                        v_costs.push(tc);
                        for _ in 1..nrun_len {
                            v_costs.push(0);
                        }
                    } else {
                        for _ in 1..nrun_len {
                            v_costs.push(0);
                        }
                        v_costs.push(tc);
                    }
                    i += nrun_len as usize;
                    no_prev_literals = 0;
                } else {
                    // Single literal: cost is 1
                    v_costs.push(1);
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
                        v_costs.pop();
                    }
                    i -= len_bck as usize;
                    pred_pos -= len_bck;
                }

                let total_len = len_bck + len_fwd;
                let tc = self.coding_cost_match(match_pos - len_bck, total_len, pred_pos);

                if prefix_costs {
                    v_costs.push(tc);
                    for _ in 1..total_len {
                        v_costs.push(0);
                    }
                } else {
                    for _ in 1..total_len {
                        v_costs.push(0);
                    }
                    v_costs.push(tc);
                }

                pred_pos = match_pos - len_bck + total_len;
                i += total_len as usize;
                no_prev_literals = 0;
            } else {
                // No match, literal cost is 1
                v_costs.push(1);
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
            }
        }

        // Remaining bases are literals
        while i < text_size {
            v_costs.push(1);
            i += 1;
        }

        v_costs
    }

    /// Compute decimal digit length like C++ int_len()
    fn int_len(x: u32) -> u32 {
        if x < 10 { 1 }
        else if x < 100 { 2 }
        else if x < 1_000 { 3 }
        else if x < 10_000 { 4 }
        else if x < 100_000 { 5 }
        else if x < 1_000_000 { 6 }
        else if x < 10_000_000 { 7 }
        else if x < 100_000_000 { 8 }
        else if x < 1_000_000_000 { 9 }
        else { 10 }
    }

    /// Compute coding cost for N-run (matches C++ coding_cost_Nrun)
    fn coding_cost_nrun(&self, len: u32) -> u32 {
        let delta = len - MIN_NRUN_LEN;
        // starter + decimal digits + suffix
        1 + Self::int_len(delta) + 1
    }

    /// Compute coding cost for match (matches C++ coding_cost_match)
    fn coding_cost_match(&self, match_pos: u32, len: u32, pred_pos: u32) -> u32 {
        let dif_pos = (match_pos as i32) - (pred_pos as i32);
        let pos_digits = if dif_pos >= 0 {
            Self::int_len(dif_pos as u32)
        } else {
            Self::int_len((-dif_pos) as u32) + 1 // sign
        };

        let delta = len - self.min_match_len;
        let len_digits = Self::int_len(delta);

        pos_digits + len_digits + 2 // pos + ',' + len + '.'
    }

    /// Compute uint_len like C++ CLZDiff_V2::uint_len (caps at 8 digits)
    fn uint_len_v2(x: u32) -> u32 {
        if x < 10 { 1 }
        else if x < 100 { 2 }
        else if x < 1_000 { 3 }
        else if x < 10_000 { 4 }
        else if x < 100_000 { 5 }
        else if x < 1_000_000 { 6 }
        else if x < 10_000_000 { 7 }
        else { 8 }
    }

    /// Compute int_len like C++ CLZDiff_V2::int_len
    fn int_len_v2(x: i32) -> u32 {
        if x >= 0 {
            Self::uint_len_v2(x as u32)
        } else {
            1 + Self::uint_len_v2((-x) as u32)
        }
    }

    /// Compute cost_match like C++ CLZDiff_V2::cost_match
    /// Note: len == u32::MAX means "match to end of sequence" (no length encoding)
    fn cost_match_v2(&self, ref_pos: u32, len: u32, pred_pos: u32) -> u32 {
        let dif_pos = (ref_pos as i32) - (pred_pos as i32);
        let mut r = Self::int_len_v2(dif_pos);

        if len != u32::MAX {
            r += 1 + Self::uint_len_v2(len - self.min_match_len);
        }

        r + 1 // +1 for '.' terminator
    }

    /// Compute cost_Nrun like C++ CLZDiff_V2::cost_Nrun
    fn cost_nrun_v2(len: u32) -> u32 {
        2 + Self::uint_len_v2(len - MIN_NRUN_LEN)
    }

    /// Estimate encoding cost without actually encoding (matches C++ CLZDiff_V2::Estimate)
    /// This is faster than full encode and used for terminator grouping decisions.
    ///
    /// # Arguments
    /// * `target` - Target sequence to estimate encoding cost for
    /// * `bound` - Early termination bound (return early if cost exceeds this)
    ///
    /// # Returns
    /// Estimated encoding cost in bytes
    pub fn estimate(&self, target: &Contig, bound: u32) -> u32 {
        if self.ht_lp.is_empty() {
            return target.len() as u32; // No index, cost is all literals
        }

        let text_size = target.len() as u32;

        // Quick check for equal sequences
        if text_size == self.reference_len as u32 {
            if target.iter().zip(self.reference.iter()).all(|(a, b)| a == b) {
                return 0; // Equal sequences
            }
        }

        let mut est_cost = 0u32;
        let mut i = 0u32;
        let mut pred_pos = 0u32;
        let mut no_prev_literals = 0u32;
        let mut x_prev: Option<u64> = None;
        let text_ptr = target.as_slice();

        while (i + self.key_len) < text_size {
            // Early termination
            if est_cost > bound {
                return est_cost;
            }

            // Get k-mer code
            let x = if x_prev.is_some() && no_prev_literals > 0 {
                self.get_code_skip1(x_prev.unwrap(), &text_ptr[i as usize..])
            } else {
                self.get_code(&text_ptr[i as usize..])
            };
            x_prev = x;

            if x.is_none() {
                // Check for N-run
                let nrun_len = self.get_nrun_len(&text_ptr[i as usize..], (text_size - i) as usize);

                if nrun_len >= MIN_NRUN_LEN {
                    est_cost += Self::cost_nrun_v2(nrun_len);
                    i += nrun_len;
                    no_prev_literals = 0;
                } else {
                    // Single literal
                    est_cost += 1;
                    i += 1;
                    pred_pos += 1;
                    no_prev_literals += 1;
                }
                continue;
            }

            // Look up k-mer in linear-probing hash table
            let hash = MurMur64Hash::hash(x.unwrap());
            let max_len = (text_size - i) as usize;

            if let Some((match_pos, len_bck, len_fwd)) =
                self.find_best_match_lp(hash, target, i as usize, max_len, no_prev_literals as usize)
            {
                let total_len = len_bck + len_fwd;
                // CRITICAL: C++ AGC's Estimate uses match_pos directly (NOT adjusted by len_bck)
                // This differs from the actual encode which does adjust for backward extension

                // Check if this is a match to end of sequence
                // C++ AGC: i + len_bck + len_fwd == text_size && match_pos + len_bck + len_fwd == reference.size() - key_len
                let is_end_match = (i + total_len) == text_size
                    && (match_pos + total_len) as usize == self.reference_len;

                if is_end_match {
                    est_cost += self.cost_match_v2(match_pos, u32::MAX, pred_pos);
                } else {
                    est_cost += self.cost_match_v2(match_pos, total_len, pred_pos);
                }

                // C++ AGC: pred_pos = match_pos + len_bck + len_fwd (NOT adjusted)
                pred_pos = match_pos + total_len;
                i += total_len;
                no_prev_literals = 0;
            } else {
                // No match, literal cost is 1
                est_cost += 1;
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
            }
        }

        // Remaining bases are literals
        est_cost += text_size - i;

        est_cost
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
