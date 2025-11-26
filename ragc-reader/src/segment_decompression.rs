// Segment Decompression
// ZSTD decompression with tuple unpacking

use anyhow::Result;
use ragc_core::Contig;

/// Unpack tuples back to bytes (matches C++ AGC tuples2bytes)
fn tuples_to_bytes(tuples: &[u8]) -> Vec<u8> {
    if tuples.is_empty() {
        return Vec::new();
    }

    let marker = tuples[tuples.len() - 1];
    let no_bytes = marker >> 4;
    let trailing_bytes = marker & 0xf;

    if no_bytes == 1 {
        return tuples[..tuples.len() - 1].to_vec();
    }

    let output_size = (tuples.len() - 2) * (no_bytes as usize) + (trailing_bytes as usize);
    let mut result = vec![0u8; output_size];

    match no_bytes {
        2 => unpack_tuples::<2, 16>(&tuples[..tuples.len() - 1], &mut result, output_size),
        3 => unpack_tuples::<3, 6>(&tuples[..tuples.len() - 1], &mut result, output_size),
        4 => unpack_tuples::<4, 4>(&tuples[..tuples.len() - 1], &mut result, output_size),
        _ => panic!("Invalid no_bytes: {no_bytes}"),
    }

    result
}

fn unpack_tuples<const N: usize, const MAX: u8>(
    tuples: &[u8],
    output: &mut [u8],
    output_size: usize,
) {
    let mut i = 0;
    let mut j = 0;

    while j + N <= output_size {
        let mut c = tuples[i] as u32;
        for k in (0..N).rev() {
            output[j + k] = (c % (MAX as u32)) as u8;
            c /= MAX as u32;
        }
        i += 1;
        j += N;
    }

    let n = output_size % N;
    if n > 0 {
        let mut c = tuples[i] as u32;
        for k in (0..n).rev() {
            output[j + k] = (c % (MAX as u32)) as u8;
            c /= MAX as u32;
        }
    }
}

/// Decompress a segment based on marker byte
/// - Marker 0: plain ZSTD decompression
/// - Marker 1+: ZSTD + tuple unpacking
pub fn decompress_segment_with_marker(compressed: &[u8], marker: u8) -> Result<Contig> {
    if compressed.is_empty() {
        return Ok(Vec::new());
    }

    if marker == 0 {
        zstd::decode_all(compressed)
            .map_err(|e| anyhow::anyhow!("ZSTD decompression failed: {e}"))
    } else {
        let tuples = zstd::decode_all(compressed)
            .map_err(|e| anyhow::anyhow!("ZSTD decompression failed: {e}"))?;
        Ok(tuples_to_bytes(&tuples))
    }
}
