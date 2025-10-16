// Variable-length integer encoding
// Rust equivalent of write/read template methods in src/common/archive.h

use std::io::{self, Write, Read};

/// Write a value with variable-length encoding
/// Format: [num_bytes: u8][value bytes in big-endian]
/// Returns number of bytes written
pub fn write_varint<W: Write>(writer: &mut W, value: u64) -> io::Result<usize> {
    // Count number of bytes needed
    let mut no_bytes = 0u8;
    let mut tmp = value;
    while tmp > 0 {
        no_bytes += 1;
        tmp >>= 8;
    }

    // Special case: value is 0
    if no_bytes == 0 {
        writer.write_all(&[0])?;
        return Ok(1);
    }

    // Write number of bytes
    writer.write_all(&[no_bytes])?;

    // Write bytes in big-endian order (most significant first)
    for i in (0..no_bytes).rev() {
        let byte = ((value >> (i * 8)) & 0xff) as u8;
        writer.write_all(&[byte])?;
    }

    Ok((no_bytes + 1) as usize)
}

/// Read a value with variable-length encoding
/// Returns (value, bytes_read)
pub fn read_varint<R: Read>(reader: &mut R) -> io::Result<(u64, usize)> {
    // Read number of bytes
    let mut no_bytes_buf = [0u8; 1];
    reader.read_exact(&mut no_bytes_buf)?;
    let no_bytes = no_bytes_buf[0];

    // Special case: value is 0
    if no_bytes == 0 {
        return Ok((0, 1));
    }

    // Read bytes and assemble value
    let mut value = 0u64;
    for _ in 0..no_bytes {
        let mut byte_buf = [0u8; 1];
        reader.read_exact(&mut byte_buf)?;
        value <<= 8;
        value += byte_buf[0] as u64;
    }

    Ok((value, (no_bytes + 1) as usize))
}

/// Write a fixed 8-byte unsigned integer
pub fn write_fixed_u64<W: Write>(writer: &mut W, value: u64) -> io::Result<usize> {
    let bytes = value.to_le_bytes();
    writer.write_all(&bytes)?;
    Ok(8)
}

/// Read a fixed 8-byte unsigned integer
pub fn read_fixed_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

/// Write a varint to a byte vector (convenience wrapper)
pub fn encode_varint(value: u64) -> Vec<u8> {
    let mut buf = Vec::new();
    write_varint(&mut buf, value).expect("Writing to Vec should not fail");
    buf
}

/// Read a varint from a byte slice (convenience wrapper)
pub fn decode_varint(bytes: &[u8]) -> io::Result<(u64, usize)> {
    let mut cursor = std::io::Cursor::new(bytes);
    read_varint(&mut cursor)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_varint_roundtrip() {
        let test_values = vec![
            0u64,
            1,
            127,
            128,
            255,
            256,
            65535,
            65536,
            0xFFFFFFFF,
            0x1234567890ABCDEF,
            u64::MAX,
        ];

        for value in test_values {
            let encoded = encode_varint(value);
            let (decoded, _) = decode_varint(&encoded).unwrap();
            assert_eq!(value, decoded, "Roundtrip failed for value {}", value);
        }
    }

    #[test]
    fn test_varint_encoding_lengths() {
        assert_eq!(encode_varint(0).len(), 1); // [0]
        assert_eq!(encode_varint(1).len(), 2); // [1, 0x01]
        assert_eq!(encode_varint(255).len(), 2); // [1, 0xff]
        assert_eq!(encode_varint(256).len(), 3); // [2, 0x01, 0x00]
        assert_eq!(encode_varint(65535).len(), 3); // [2, 0xff, 0xff]
        assert_eq!(encode_varint(65536).len(), 4); // [3, 0x01, 0x00, 0x00]
    }

    #[test]
    fn test_fixed_u64_roundtrip() {
        let test_values = vec![0u64, 1, 42, 0xDEADBEEF, u64::MAX];

        for value in test_values {
            let mut buf = Vec::new();
            write_fixed_u64(&mut buf, value).unwrap();
            assert_eq!(buf.len(), 8);

            let decoded = read_fixed_u64(&mut std::io::Cursor::new(&buf)).unwrap();
            assert_eq!(value, decoded);
        }
    }
}
