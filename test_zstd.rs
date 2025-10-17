use zstd::stream::{encode_all, decode_all};

fn main() {
    // Test data - 70032 bytes (70031 + separator)
    let mut data = vec![0u8; 70031];
    for (i, byte) in data.iter_mut().enumerate() {
        *byte = (i % 256) as u8;
    }
    data.push(0xFF); // separator
    
    println!("Original data len: {}", data.len());
    
    // Compress
    let compressed = encode_all(&data[..], 3).expect("Compression failed");
    println!("Compressed len: {}", compressed.len());
    
    // Add marker byte (like our code does)
    let mut compressed_with_marker = compressed;
    compressed_with_marker.push(0);
    println!("Compressed with marker len: {}", compressed_with_marker.len());
    
    // Decompress (removing marker byte first)
    let to_decompress = &compressed_with_marker[..compressed_with_marker.len()-1];
    let decompressed = decode_all(to_decompress).expect("Decompression failed");
    println!("Decompressed len: {}", decompressed.len());
    
    // Verify
    if decompressed == data {
        println!("SUCCESS: Data matches!");
    } else {
        println!("FAIL: Data does not match!");
    }
}
