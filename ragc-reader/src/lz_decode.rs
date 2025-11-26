// LZ Diff Decoding (reader-only)
// Decodes a target sequence from differences against a reference

use ragc_core::Contig;

const N_CODE: u8 = 4;
const N_RUN_STARTER_CODE: u8 = 30;
const MIN_NRUN_LEN: u32 = 4;

/// LZ Diff decoder
pub struct LZDecoder {
    reference: Vec<u8>,
    reference_len: usize,
    min_match_len: u32,
}

impl LZDecoder {
    /// Create a new LZ decoder with the given minimum match length
    pub fn new(min_match_len: u32) -> Self {
        LZDecoder {
            reference: Vec::new(),
            reference_len: 0,
            min_match_len,
        }
    }

    /// Prepare the decoder with a reference sequence
    pub fn prepare(&mut self, reference: &Contig) {
        self.reference = reference.clone();
        self.reference_len = reference.len();
        // Add padding for key_len (matches encoder behavior)
        let key_len = self.min_match_len - 3; // HASHING_STEP=4, so key_len = min_match_len - 4 + 1
        self.reference.resize(self.reference.len() + key_len as usize, 31);
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

    fn is_literal(&self, c: u8) -> bool {
        (b'A'..=b'A' + 20).contains(&c) || c == b'!'
    }

    fn decode_literal(&self, c: u8) -> u8 {
        if c == b'!' { b'!' } else { c - b'A' }
    }

    fn decode_nrun(&self, data: &[u8]) -> (u32, usize) {
        let mut i = 1;
        let (raw_len, len_bytes) = self.read_int(&data[i..]);
        i += len_bytes;
        i += 1;
        ((raw_len as u32) + MIN_NRUN_LEN, i)
    }

    fn decode_match(&self, data: &[u8], pred_pos: usize) -> (usize, u32, usize) {
        let mut i = 0;
        let (raw_pos, pos_bytes) = self.read_int(&data[i..]);
        i += pos_bytes;

        let ref_pos = ((pred_pos as i64) + raw_pos) as usize;

        let len = if data[i] == b',' {
            i += 1;
            let (raw_len, len_bytes) = self.read_int(&data[i..]);
            i += len_bytes;
            i += 1;
            (raw_len as u32) + self.min_match_len
        } else {
            i += 1;
            u32::MAX
        };

        (ref_pos, len, i)
    }

    fn read_int(&self, data: &[u8]) -> (i64, usize) {
        let mut i = 0;
        let mut is_neg = false;

        if !data.is_empty() && data[i] == b'-' {
            is_neg = true;
            i += 1;
        }

        let mut x = 0i64;
        while i < data.len() && data[i] >= b'0' && data[i] <= b'9' {
            x = x * 10 + ((data[i] - b'0') as i64);
            i += 1;
        }

        if is_neg { x = -x; }
        (x, i)
    }
}
