// Integration test for variable-length integer encoding
// Should produce output identical to C++ test_varint

use ragc_common::{write_varint, read_varint, write_fixed_u64, read_fixed_u64};
use std::io::Cursor;

fn print_bytes(bytes: &[u8]) {
    for (i, byte) in bytes.iter().enumerate() {
        if i > 0 {
            print!(" ");
        }
        print!("{:02x}", byte);
    }
}

fn main() {
    // Test 1: Varint encoding for various values
    println!("# Test 1: Varint encoding");

    let test_values: Vec<u64> = vec![
        0,
        1,
        127,
        128,
        255,
        256,
        65535,
        65536,
        0xFFFFFFFF,
        0x1234567890ABCDEF,
        0xFFFFFFFFFFFFFFFF,
    ];

    for val in &test_values {
        let mut buf = Vec::new();
        let bytes_written = write_varint(&mut buf, *val).unwrap();
        print!("{:x}\t{}\t", val, bytes_written);
        print_bytes(&buf);
        println!();
    }

    // Test 2: Varint roundtrip (encode + decode)
    println!("\n# Test 2: Varint roundtrip");
    for val in &test_values {
        let mut buf = Vec::new();
        write_varint(&mut buf, *val).unwrap();

        let mut cursor = Cursor::new(&buf);
        let (decoded, _bytes_read) = read_varint(&mut cursor).unwrap();

        println!(
            "{:x}\t{:x}\t{}",
            val,
            decoded,
            if *val == decoded { "OK" } else { "FAIL" }
        );
    }

    // Test 3: Fixed u64 encoding
    println!("\n# Test 3: Fixed u64 encoding");

    let fixed_test_values: Vec<u64> = vec![
        0,
        1,
        42,
        0xDEADBEEF,
        0xFFFFFFFFFFFFFFFF,
    ];

    for val in &fixed_test_values {
        let mut buf = Vec::new();
        write_fixed_u64(&mut buf, *val).unwrap();
        print!("{:x}\t8\t", val);
        print_bytes(&buf);
        println!();
    }

    // Test 4: Fixed u64 roundtrip
    println!("\n# Test 4: Fixed u64 roundtrip");
    for val in &fixed_test_values {
        let mut buf = Vec::new();
        write_fixed_u64(&mut buf, *val).unwrap();

        let mut cursor = Cursor::new(&buf);
        let decoded = read_fixed_u64(&mut cursor).unwrap();

        println!(
            "{:x}\t{:x}\t{}",
            val,
            decoded,
            if *val == decoded { "OK" } else { "FAIL" }
        );
    }
}
