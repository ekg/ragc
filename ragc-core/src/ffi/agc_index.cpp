// FFI bridge to C++ AGC index writing (Archive + Collection_V3)
// Provides C-compatible interface for Rust

#include <memory>
#include <vector>
#include <string>
#include <cstring>

// C++ AGC headers (from /home/erik/ragc/agc/src/common/)
#include "archive.h"
#include "collection_v3.h"
#include "defs.h"

extern "C" {

// ========== Archive Functions ==========

void* agc_archive_create_writer() {
    auto* archive = new std::shared_ptr<CArchive>(std::make_shared<CArchive>(false));
    return archive;
}

bool agc_archive_open(void* archive_ptr, const char* path) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    return (*archive)->Open(std::string(path));
}

int agc_archive_register_stream(void* archive_ptr, const char* name) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    return (*archive)->RegisterStream(std::string(name));
}

bool agc_archive_add_part(
    void* archive_ptr,
    int stream_id,
    const uint8_t* data,
    size_t data_len,
    uint64_t metadata
) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    std::vector<uint8_t> data_vec(data, data + data_len);
    return (*archive)->AddPartBuffered(stream_id, data_vec, metadata);
}

bool agc_archive_close(void* archive_ptr) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    return (*archive)->Close();
}

void agc_archive_destroy(void* archive_ptr) {
    auto* archive = static_cast<std::shared_ptr<CArchive>*>(archive_ptr);
    delete archive;
}

// ========== Collection_V3 Functions ==========

void* agc_collection_create() {
    auto* coll = new std::shared_ptr<CCollection_V3>(std::make_shared<CCollection_V3>());
    return coll;
}

bool agc_collection_set_archives(
    void* coll_ptr,
    void* out_archive_ptr,
    uint32_t num_threads,
    size_t batch_size,
    uint32_t segment_size,
    uint32_t kmer_length
) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);
    auto* out_archive = static_cast<std::shared_ptr<CArchive>*>(out_archive_ptr);

    // in_archive is nullptr for writing
    return (*coll)->set_archives(nullptr, *out_archive, num_threads, batch_size, segment_size, kmer_length);
}

void agc_collection_complete_serialization(void* coll_ptr) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);
    (*coll)->complete_serialization();
}

void agc_collection_store_contig_batch(void* coll_ptr, uint32_t id_from, uint32_t id_to) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);
    (*coll)->store_contig_batch(id_from, id_to);
}

size_t agc_collection_get_no_samples(void* coll_ptr) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);
    return (*coll)->get_no_samples();
}

void agc_collection_destroy(void* coll_ptr) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);
    delete coll;
}

// ========== Batch Segment Placement ==========

struct CSegmentPlacement {
    const char* sample_name;
    const char* contig_name;
    uint32_t seg_part_no;
    int32_t group_id;
    int32_t in_group_id;
    bool is_rev_comp;
    uint32_t data_size;
};

bool agc_collection_register_sample_contig(
    void* coll_ptr,
    const char* sample_name,
    const char* contig_name
) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);
    return (*coll)->register_sample_contig(
        std::string(sample_name),
        std::string(contig_name)
    );
}

void agc_collection_add_segments_placed(
    void* coll_ptr,
    const CSegmentPlacement* placements,
    size_t count
) {
    auto* coll = static_cast<std::shared_ptr<CCollection_V3>*>(coll_ptr);

    // Convert C array to C++ vector
    std::vector<segments_to_place_t> segments;
    segments.reserve(count);

    for (size_t i = 0; i < count; ++i) {
        const auto& p = placements[i];
        segments.emplace_back(
            std::string(p.sample_name),
            std::string(p.contig_name),
            p.seg_part_no,
            p.group_id,
            p.in_group_id,
            p.is_rev_comp,
            p.data_size
        );
    }

    (*coll)->add_segments_placed(segments);
}

} // extern "C"
