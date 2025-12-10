//! Cached environment variable lookups for debug flags
//! These are checked once at startup and cached, avoiding ~30% CPU overhead from getenv calls

use std::sync::OnceLock;

// Each flag has a unique static cache and a public accessor function

static DEBUG_LZ_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_lz() -> bool {
    *DEBUG_LZ_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_LZ").is_ok())
}

static DEBUG_REF_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_ref() -> bool {
    *DEBUG_REF_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_REF").is_ok())
}

static DEBUG_LZ_DECODE_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_lz_decode() -> bool {
    *DEBUG_LZ_DECODE_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_LZ_DECODE").is_ok())
}

static DEBUG_LZ_DECODE_FULL_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_lz_decode_full() -> bool {
    *DEBUG_LZ_DECODE_FULL_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_LZ_DECODE_FULL").is_ok())
}

static DEBUG_REF_WRITE_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_ref_write() -> bool {
    *DEBUG_REF_WRITE_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_REF_WRITE").is_ok())
}

static DEBUG_COMPRESSION_SIZES_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_compression_sizes() -> bool {
    *DEBUG_COMPRESSION_SIZES_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_COMPRESSION_SIZES").is_ok())
}

static DEBUG_IS_DIR_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_is_dir() -> bool {
    *DEBUG_IS_DIR_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_IS_DIR").is_ok())
}

static DEBUG_SPLIT_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_split() -> bool {
    *DEBUG_SPLIT_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_SPLIT").is_ok())
}

static DEBUG_SPLIT_FIND_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_split_find() -> bool {
    *DEBUG_SPLIT_FIND_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_SPLIT_FIND").is_ok())
}

static DEBUG_SPLIT_REF_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_split_ref() -> bool {
    *DEBUG_SPLIT_REF_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_SPLIT_REF").is_ok())
}

static DEBUG_SPLIT_MAP_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_split_map() -> bool {
    *DEBUG_SPLIT_MAP_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_SPLIT_MAP").is_ok())
}

static DEBUG_TERM_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_term() -> bool {
    *DEBUG_TERM_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_TERM").is_ok())
}

static DEBUG_FALLBACK2_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_fallback2() -> bool {
    *DEBUG_FALLBACK2_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_FALLBACK2").is_ok())
}

static DEBUG_SEGMENT_COVERAGE_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_segment_coverage() -> bool {
    *DEBUG_SEGMENT_COVERAGE_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_SEGMENT_COVERAGE").is_ok())
}

static DEBUG_ENDPOS_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_endpos() -> bool {
    *DEBUG_ENDPOS_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_ENDPOS").is_ok())
}

static DEBUG_OVERLAP_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_overlap() -> bool {
    *DEBUG_OVERLAP_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_OVERLAP").is_ok())
}

static TRACE_ALL_SPLITS_CACHE: OnceLock<bool> = OnceLock::new();
pub fn trace_all_splits() -> bool {
    *TRACE_ALL_SPLITS_CACHE.get_or_init(|| std::env::var("RAGC_TRACE_ALL_SPLITS").is_ok())
}

static TRACE_GROUPS_CACHE: OnceLock<bool> = OnceLock::new();
pub fn trace_groups() -> bool {
    *TRACE_GROUPS_CACHE.get_or_init(|| std::env::var("RAGC_TRACE_GROUPS").is_ok())
}

static TRACE_GROUP_CACHE: OnceLock<bool> = OnceLock::new();
pub fn trace_group() -> bool {
    *TRACE_GROUP_CACHE.get_or_init(|| std::env::var("RAGC_TRACE_GROUP").is_ok())
}

static TEST_LZ_ENCODING_CACHE: OnceLock<bool> = OnceLock::new();
pub fn test_lz_encoding() -> bool {
    *TEST_LZ_ENCODING_CACHE.get_or_init(|| std::env::var("RAGC_TEST_LZ_ENCODING").is_ok())
}

static ASSERT_VERBOSE_CACHE: OnceLock<bool> = OnceLock::new();
pub fn assert_verbose() -> bool {
    *ASSERT_VERBOSE_CACHE.get_or_init(|| std::env::var("RAGC_ASSERT_VERBOSE").is_ok())
}

static DEBUG_DECODE_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_decode() -> bool {
    *DEBUG_DECODE_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_DECODE").is_ok())
}

static DEBUG_RECONSTRUCT_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_reconstruct() -> bool {
    *DEBUG_RECONSTRUCT_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_RECONSTRUCT").is_ok())
}

// Boolean env vars with specific value checks
static DEBUG_LZ_ENABLED_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_lz_enabled() -> bool {
    *DEBUG_LZ_ENABLED_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_LZ").unwrap_or_default() == "1")
}

static SPLIT_ALL_CACHE: OnceLock<bool> = OnceLock::new();
pub fn split_all() -> bool {
    *SPLIT_ALL_CACHE.get_or_init(|| std::env::var("RAGC_SPLIT_ALL").unwrap_or_default() == "1")
}

static SPLIT_CREATE_GROUPS_CACHE: OnceLock<bool> = OnceLock::new();
pub fn split_create_groups() -> bool {
    *SPLIT_CREATE_GROUPS_CACHE.get_or_init(|| std::env::var("RAGC_SPLIT_CREATE_GROUPS").unwrap_or_default() == "1")
}

static GROUP_LOG_CACHE: OnceLock<bool> = OnceLock::new();
pub fn group_log() -> bool {
    *GROUP_LOG_CACHE.get_or_init(|| std::env::var("RAGC_GROUP_LOG").unwrap_or_default() == "1")
}

static DEBUG_FALLBACK2_ENABLED_CACHE: OnceLock<bool> = OnceLock::new();
pub fn debug_fallback2_enabled() -> bool {
    *DEBUG_FALLBACK2_ENABLED_CACHE.get_or_init(|| std::env::var("RAGC_DEBUG_FALLBACK2").unwrap_or_default() == "1")
}

// Optional: cache the assert path (less hot but still called)
static ASSERT_CPP_ARCHIVE_CACHE: OnceLock<Option<String>> = OnceLock::new();
pub fn assert_cpp_archive() -> Option<&'static str> {
    ASSERT_CPP_ARCHIVE_CACHE.get_or_init(|| std::env::var("RAGC_ASSERT_CPP_ARCHIVE").ok())
        .as_deref()
}
