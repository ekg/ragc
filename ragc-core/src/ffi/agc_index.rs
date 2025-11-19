// Rust FFI bindings to C++ AGC index writing (Archive + Collection_V3)

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::ptr;

// C++ bridge functions (from agc_index.cpp)
extern "C" {
    fn agc_archive_create_writer() -> *mut std::ffi::c_void;
    fn agc_archive_open(archive_ptr: *mut std::ffi::c_void, path: *const c_char) -> bool;
    fn agc_archive_register_stream(archive_ptr: *mut std::ffi::c_void, name: *const c_char) -> i32;
    fn agc_archive_add_part(
        archive_ptr: *mut std::ffi::c_void,
        stream_id: i32,
        data: *const u8,
        data_len: usize,
        metadata: u64,
    ) -> bool;
    fn agc_archive_close(archive_ptr: *mut std::ffi::c_void) -> bool;
    fn agc_archive_destroy(archive_ptr: *mut std::ffi::c_void);

    fn agc_collection_create() -> *mut std::ffi::c_void;
    fn agc_collection_set_archives(
        coll_ptr: *mut std::ffi::c_void,
        out_archive_ptr: *mut std::ffi::c_void,
        num_threads: u32,
        batch_size: usize,
        segment_size: u32,
        kmer_length: u32,
    ) -> bool;
    fn agc_collection_complete_serialization(coll_ptr: *mut std::ffi::c_void);
    fn agc_collection_store_contig_batch(
        coll_ptr: *mut std::ffi::c_void,
        id_from: u32,
        id_to: u32,
    );
    fn agc_collection_get_no_samples(coll_ptr: *mut std::ffi::c_void) -> usize;
    fn agc_collection_destroy(coll_ptr: *mut std::ffi::c_void);

    fn agc_collection_register_sample_contig(
        coll_ptr: *mut std::ffi::c_void,
        sample_name: *const c_char,
        contig_name: *const c_char,
    ) -> bool;

    fn agc_collection_add_segments_placed(
        coll_ptr: *mut std::ffi::c_void,
        placements: *const CSegmentPlacement,
        count: usize,
    );
}

#[repr(C)]
struct CSegmentPlacement {
    sample_name: *const c_char,
    contig_name: *const c_char,
    seg_part_no: u32,
    group_id: i32,
    in_group_id: i32,
    is_rev_comp: bool,
    data_size: u32,
}

/// Rust-friendly segment placement structure
pub struct SegmentPlacement {
    pub sample_name: String,
    pub contig_name: String,
    pub seg_part_no: u32,
    pub group_id: u32,
    pub in_group_id: u32,
    pub is_rev_comp: bool,
    pub data_size: u32,
}

/// Safe wrapper around C++ AGC Archive
pub struct CppArchive {
    ptr: *mut std::ffi::c_void,
}

impl CppArchive {
    pub fn new_writer() -> Self {
        unsafe {
            let ptr = agc_archive_create_writer();
            assert!(!ptr.is_null(), "Failed to create C++ archive");
            Self { ptr }
        }
    }

    pub fn open(&mut self, path: &str) -> anyhow::Result<()> {
        let c_path = CString::new(path)?;
        unsafe {
            if agc_archive_open(self.ptr, c_path.as_ptr()) {
                Ok(())
            } else {
                Err(anyhow::anyhow!("Failed to open archive: {}", path))
            }
        }
    }

    pub fn register_stream(&mut self, name: &str) -> anyhow::Result<i32> {
        let c_name = CString::new(name)?;
        unsafe {
            let stream_id = agc_archive_register_stream(self.ptr, c_name.as_ptr());
            if stream_id >= 0 {
                Ok(stream_id)
            } else {
                Err(anyhow::anyhow!("Failed to register stream: {}", name))
            }
        }
    }

    pub fn add_part(&mut self, stream_id: i32, data: &[u8], metadata: u64) -> anyhow::Result<()> {
        unsafe {
            if agc_archive_add_part(self.ptr, stream_id, data.as_ptr(), data.len(), metadata) {
                Ok(())
            } else {
                Err(anyhow::anyhow!("Failed to add part to stream {}", stream_id))
            }
        }
    }

    pub fn close(&mut self) -> anyhow::Result<()> {
        unsafe {
            if agc_archive_close(self.ptr) {
                Ok(())
            } else {
                Err(anyhow::anyhow!("Failed to close archive"))
            }
        }
    }

    pub fn as_ptr(&mut self) -> *mut std::ffi::c_void {
        self.ptr
    }
}

impl Drop for CppArchive {
    fn drop(&mut self) {
        unsafe {
            agc_archive_destroy(self.ptr);
        }
    }
}

// Archive must be Send + Sync for multi-threading
unsafe impl Send for CppArchive {}
unsafe impl Sync for CppArchive {}

/// Safe wrapper around C++ AGC Collection_V3
pub struct CppCollection {
    ptr: *mut std::ffi::c_void,
}

impl CppCollection {
    pub fn new() -> Self {
        unsafe {
            let ptr = agc_collection_create();
            assert!(!ptr.is_null(), "Failed to create C++ collection");
            Self { ptr }
        }
    }

    pub fn set_archives(
        &mut self,
        archive: &mut CppArchive,
        num_threads: u32,
        batch_size: usize,
        segment_size: u32,
        kmer_length: u32,
    ) -> anyhow::Result<()> {
        unsafe {
            if agc_collection_set_archives(
                self.ptr,
                archive.as_ptr(),
                num_threads,
                batch_size,
                segment_size,
                kmer_length,
            ) {
                Ok(())
            } else {
                Err(anyhow::anyhow!("Failed to set archives"))
            }
        }
    }

    pub fn register_sample_contig(&mut self, sample_name: &str, contig_name: &str) -> anyhow::Result<()> {
        let c_sample_name = CString::new(sample_name)?;
        let c_contig_name = CString::new(contig_name)?;

        unsafe {
            if agc_collection_register_sample_contig(
                self.ptr,
                c_sample_name.as_ptr(),
                c_contig_name.as_ptr(),
            ) {
                Ok(())
            } else {
                Err(anyhow::anyhow!(
                    "Failed to register sample '{}' contig '{}'",
                    sample_name,
                    contig_name
                ))
            }
        }
    }

    pub fn add_segments_placed(&mut self, placements: &[SegmentPlacement]) -> anyhow::Result<()> {
        if placements.is_empty() {
            return Ok(());
        }

        // Convert Rust strings to C strings (must keep them alive during FFI call)
        let c_sample_names: Vec<CString> = placements
            .iter()
            .map(|p| CString::new(p.sample_name.as_str()).unwrap())
            .collect();

        let c_contig_names: Vec<CString> = placements
            .iter()
            .map(|p| CString::new(p.contig_name.as_str()).unwrap())
            .collect();

        // Build C-compatible array
        let c_placements: Vec<CSegmentPlacement> = placements
            .iter()
            .enumerate()
            .map(|(i, p)| CSegmentPlacement {
                sample_name: c_sample_names[i].as_ptr(),
                contig_name: c_contig_names[i].as_ptr(),
                seg_part_no: p.seg_part_no,
                group_id: p.group_id as i32,
                in_group_id: p.in_group_id as i32,
                is_rev_comp: p.is_rev_comp,
                data_size: p.data_size,
            })
            .collect();

        unsafe {
            agc_collection_add_segments_placed(self.ptr, c_placements.as_ptr(), c_placements.len());
        }

        Ok(())
    }

    pub fn complete_serialization(&mut self) {
        unsafe {
            agc_collection_complete_serialization(self.ptr);
        }
    }

    pub fn store_contig_batch(&mut self, id_from: u32, id_to: u32) {
        unsafe {
            agc_collection_store_contig_batch(self.ptr, id_from, id_to);
        }
    }

    pub fn get_no_samples(&self) -> usize {
        unsafe {
            agc_collection_get_no_samples(self.ptr)
        }
    }
}

impl Drop for CppCollection {
    fn drop(&mut self) {
        unsafe {
            agc_collection_destroy(self.ptr);
        }
    }
}

// Collection must be Send + Sync for multi-threading
unsafe impl Send for CppCollection {}
unsafe impl Sync for CppCollection {}
