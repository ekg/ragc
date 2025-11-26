// Contig preprocessing functions matching C++ AGC genome_io logic

/// Convert table matching C++ AGC's cnv_num[128]
/// Maps ASCII characters to numeric DNA codes
const CNV_NUM: [u8; 128] = [
    // 0-15: Standard DNA codes
    b'A', b'C', b'G', b'T', b'N', b'R', b'Y', b'S', b'W', b'K', b'M', b'B', b'D', b'H', b'V', b'U',
    // 16-63: Spaces
    b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
    b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
    b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
    // 64-79: ASCII uppercase letters
    b' ',   0,  11,   1,  12,  30,  30,   2,  13,  30,  30,   9,  30,  10,   4,  30,
    // 80-95: ASCII uppercase letters continued
     30,  30,   5,   7,   3,  15,  14,   8,  30,   6,  30,  30,  30,  30,  30,  30,
    // 96-111: ASCII lowercase letters
    b' ',   0,  11,   1,  12,  30,  30,   2,  13,  30,  30,   9,  30,  10,   4,  30,
    // 112-127: ASCII lowercase letters continued
     30,  30,   5,   7,   3,  15,  14,   8,  30,   6,  30,  30,  30,  30,  30,  30
];

/// Preprocess raw contig by converting ASCII characters to numeric codes
///
/// Matches C++ AGC's preprocess_raw_contig() logic:
/// - Filters bytes with high bits set (>= 64, i.e., letter characters)
/// - Converts them using cnv_num lookup table
/// - Removes all other characters (like spaces, newlines)
/// - Uses loop unrolling for performance (processes 4 bytes at a time)
///
/// # Arguments
/// * `contig` - Mutable vector to process in place
///
/// # Behavior
/// - Modifies contig in place
/// - Resizes to contain only converted characters
pub fn preprocess_raw_contig(contig: &mut Vec<u8>) {
    let len = contig.len();
    let mut in_pos = 0usize;
    let mut out_pos = 0usize;

    // Handle remainder (len % 4) using Duff's device pattern
    match len % 4 {
        3 => {
            let c = contig[in_pos];
            in_pos += 1;
            if c >> 6 != 0 {  // c >= 64
                contig[out_pos] = CNV_NUM[c as usize];
                out_pos += 1;
            }
            // Fall through to case 2
            let c = contig[in_pos];
            in_pos += 1;
            if c >> 6 != 0 {
                contig[out_pos] = CNV_NUM[c as usize];
                out_pos += 1;
            }
            // Fall through to case 1
            let c = contig[in_pos];
            in_pos += 1;
            if c >> 6 != 0 {
                contig[out_pos] = CNV_NUM[c as usize];
                out_pos += 1;
            }
        },
        2 => {
            let c = contig[in_pos];
            in_pos += 1;
            if c >> 6 != 0 {
                contig[out_pos] = CNV_NUM[c as usize];
                out_pos += 1;
            }
            // Fall through to case 1
            let c = contig[in_pos];
            in_pos += 1;
            if c >> 6 != 0 {
                contig[out_pos] = CNV_NUM[c as usize];
                out_pos += 1;
            }
        },
        1 => {
            let c = contig[in_pos];
            in_pos += 1;
            if c >> 6 != 0 {
                contig[out_pos] = CNV_NUM[c as usize];
                out_pos += 1;
            }
        },
        _ => {} // len % 4 == 0, nothing to do
    }

    // Process remaining bytes 4 at a time (loop unrolling)
    while in_pos < len {
        let c = contig[in_pos];
        in_pos += 1;
        if c >> 6 != 0 {
            contig[out_pos] = CNV_NUM[c as usize];
            out_pos += 1;
        }

        let c = contig[in_pos];
        in_pos += 1;
        if c >> 6 != 0 {
            contig[out_pos] = CNV_NUM[c as usize];
            out_pos += 1;
        }

        let c = contig[in_pos];
        in_pos += 1;
        if c >> 6 != 0 {
            contig[out_pos] = CNV_NUM[c as usize];
            out_pos += 1;
        }

        let c = contig[in_pos];
        in_pos += 1;
        if c >> 6 != 0 {
            contig[out_pos] = CNV_NUM[c as usize];
            out_pos += 1;
        }
    }

    contig.truncate(out_pos);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_preprocess_raw_contig_acgt() {
        let mut contig = b"ACGT".to_vec();
        preprocess_raw_contig(&mut contig);
        assert_eq!(contig, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_preprocess_raw_contig_with_spaces() {
        let mut contig = b"A C G T".to_vec();
        preprocess_raw_contig(&mut contig);
        // Spaces (ASCII 32, < 64) should be filtered out
        assert_eq!(contig, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_preprocess_raw_contig_mixed() {
        let mut contig = b"ACGTN\nATGC".to_vec();
        preprocess_raw_contig(&mut contig);
        // Newline (ASCII 10, < 64) should be filtered out
        assert_eq!(contig, vec![0, 1, 2, 3, 4, 0, 3, 2, 1]);
    }
}
