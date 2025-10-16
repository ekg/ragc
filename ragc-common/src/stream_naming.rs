// Stream naming utilities for AGC archives
// Matches C++ utils.cpp stream naming conventions

/// Convert integer to custom base64 encoding (matches C++ int_to_base64)
/// Uses digits: 0-9, A-Z, a-z, _, #
pub fn int_to_base64(mut n: u32) -> String {
    const DIGITS: &[u8; 64] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_#";

    let mut result = String::new();

    loop {
        result.push(DIGITS[(n & 0x3f) as usize] as char);
        n /= 64;
        if n == 0 {
            break;
        }
    }

    result
}

/// Get stream prefix based on archive version
pub fn stream_prefix(archive_version: u32) -> &'static str {
    if archive_version < 3000 {
        "seg-"
    } else {
        "x"
    }
}

/// Get base stream name
pub fn stream_base(archive_version: u32, n: u32) -> String {
    if archive_version < 3000 {
        format!("seg-{}", n)
    } else {
        format!("x{}", int_to_base64(n))
    }
}

/// Get reference stream name
pub fn stream_ref_name(archive_version: u32, n: u32) -> String {
    if archive_version < 3000 {
        format!("seg-{}-ref", n)
    } else {
        format!("x{}r", int_to_base64(n))
    }
}

/// Get delta stream name
pub fn stream_delta_name(archive_version: u32, n: u32) -> String {
    if archive_version < 3000 {
        format!("seg-{}-delta", n)
    } else {
        format!("x{}d", int_to_base64(n))
    }
}

/// Get reference extension
pub fn stream_ref_ext(archive_version: u32) -> &'static str {
    if archive_version < 3000 {
        "-ref"
    } else {
        "r"
    }
}

/// Get delta extension
pub fn stream_delta_ext(archive_version: u32) -> &'static str {
    if archive_version < 3000 {
        "-delta"
    } else {
        "d"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int_to_base64() {
        // Test basic values
        assert_eq!(int_to_base64(0), "0");
        assert_eq!(int_to_base64(1), "1");
        assert_eq!(int_to_base64(9), "9");
        assert_eq!(int_to_base64(10), "A");
        assert_eq!(int_to_base64(35), "Z");
        assert_eq!(int_to_base64(36), "a");
        assert_eq!(int_to_base64(61), "z");
        assert_eq!(int_to_base64(62), "_");
        assert_eq!(int_to_base64(63), "#");

        // Test multi-digit values
        assert_eq!(int_to_base64(64), "01");  // 64 = 1*64 + 0
        assert_eq!(int_to_base64(65), "11");  // 65 = 1*64 + 1
    }

    #[test]
    fn test_stream_naming_v3() {
        let version = 3000;

        assert_eq!(stream_prefix(version), "x");
        assert_eq!(stream_base(version, 0), "x0");
        assert_eq!(stream_base(version, 1), "x1");
        assert_eq!(stream_base(version, 10), "xA");

        assert_eq!(stream_ref_name(version, 0), "x0r");
        assert_eq!(stream_delta_name(version, 0), "x0d");

        assert_eq!(stream_ref_ext(version), "r");
        assert_eq!(stream_delta_ext(version), "d");
    }

    #[test]
    fn test_stream_naming_legacy() {
        let version = 2000;

        assert_eq!(stream_prefix(version), "seg-");
        assert_eq!(stream_base(version, 0), "seg-0");
        assert_eq!(stream_base(version, 1), "seg-1");

        assert_eq!(stream_ref_name(version, 0), "seg-0-ref");
        assert_eq!(stream_delta_name(version, 0), "seg-0-delta");

        assert_eq!(stream_ref_ext(version), "-ref");
        assert_eq!(stream_delta_ext(version), "-delta");
    }
}
