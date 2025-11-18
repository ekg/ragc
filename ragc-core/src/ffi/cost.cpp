// Minimal C++ cost engine exposed via C ABI for Rust FFI
#include <cstdint>
#include <vector>
#include <algorithm>

namespace {
static inline uint64_t murmur64(uint64_t h) {
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33; return h;
}

struct CostEngine {
    static constexpr uint32_t HASHING_STEP = 4;
    static constexpr uint32_t MIN_NRUN_LEN = 4;
    static constexpr uint8_t N_CODE = 4;
    static constexpr uint8_t INVALID = 31;
    static constexpr uint32_t MAX_NO_TRIES = 64;

    uint32_t min_match_len;
    uint32_t key_len;
    uint64_t key_mask;
    std::vector<uint8_t> reference; // padded
    std::vector<uint32_t> ht; // LP table
    uint64_t ht_mask{};

    explicit CostEngine(uint32_t mml): min_match_len(mml) {
        key_len = min_match_len - HASHING_STEP + 1u;
        key_mask = ~0ull >> (64 - 2 * key_len);
    }

    void prepare(const uint8_t* ref, size_t ref_len) {
        reference.assign(ref, ref + ref_len);
        reference.resize(reference.size() + key_len, INVALID);
        build_index();
    }

    void build_index() {
        uint64_t sz = 0; uint32_t no_prev_valid = 0; uint32_t cnt_mod = 0; uint32_t key_len_mod = key_len % HASHING_STEP;
        for (auto c: reference) { if (c<4) ++no_prev_valid; else no_prev_valid = 0; if (++cnt_mod == HASHING_STEP) cnt_mod = 0; if (cnt_mod==key_len_mod && no_prev_valid>=key_len) ++sz; }
        uint64_t ht_size = (uint64_t)(sz / 0.7);
        while (ht_size & (ht_size - 1)) ht_size &= ht_size - 1; ht_size <<= 1; if (ht_size < 8) ht_size = 8;
        ht_mask = ht_size - 1;
        ht.assign(ht_size, UINT32_MAX);
        for (size_t i=0; i + key_len < reference.size(); i += HASHING_STEP) {
            uint64_t code = get_code(&reference[i]); if (code == ~0ull) continue;
            uint64_t base = murmur64(code) & ht_mask;
            for (uint32_t j=0;j<MAX_NO_TRIES;++j) { size_t idx=(base+j)&ht_mask; if (ht[idx]==UINT32_MAX) { ht[idx]=(uint32_t)(i/HASHING_STEP); break; } }
        }
    }

    inline uint64_t get_code(const uint8_t* s) const {
        uint64_t x=0; uint32_t i = key_len % 4;
        switch(i){case 3: if(s[0]>3) return ~0ull; x=(x<<2)+s[0]; ++s; [[fallthrough]]; case 2: if(s[0]>3) return ~0ull; x=(x<<2)+s[0]; ++s; [[fallthrough]]; case 1: if(s[0]>3) return ~0ull; x=(x<<2)+s[0]; ++s;}
        for(; i<key_len;){ if(*s>3) return ~0ull; x=(x<<2)+*s++; ++i; if(*s>3) return ~0ull; x=(x<<2)+*s++; ++i; if(*s>3) return ~0ull; x=(x<<2)+*s++; ++i; if(*s>3) return ~0ull; x=(x<<2)+*s++; ++i; }
        return x;
    }

    static inline uint32_t int_len(uint32_t x){ if(x<10) return 1; if(x<100) return 2; if(x<1000) return 3; if(x<10000) return 4; if(x<100000) return 5; if(x<1000000) return 6; if(x<10000000) return 7; if(x<100000000) return 8; if(x<1000000000) return 9; return 10; }
    inline uint32_t cost_match(uint32_t ref_pos, uint32_t len, uint32_t pred_pos) const { int dif=(int)ref_pos-(int)pred_pos; uint32_t r = dif>=0? int_len((uint32_t)dif) : int_len((uint32_t)(-dif))+1; r += int_len(len - min_match_len) + 2; return r; }
    inline uint32_t cost_nrun(uint32_t len) const { return 2 + int_len(len - MIN_NRUN_LEN); }

    std::vector<uint32_t> get_costs(const uint8_t* text, size_t text_len, bool prefix) const {
        std::vector<uint32_t> v; v.reserve(text_len);
        uint32_t i=0,pred=0,no_prev=0; uint64_t x_prev=~0ull;
        while (i + key_len < text_len) {
            uint64_t x = (x_prev!=~0ull && no_prev>0) ? get_code_skip1(x_prev, &text[i]) : get_code(&text[i]);
            x_prev = x;
            if (x == ~0ull) {
                uint32_t nrun = get_nrun_len(&text[i], (uint32_t)(text_len - i));
                if (nrun >= MIN_NRUN_LEN) { auto tc = cost_nrun(nrun); if(prefix){ v.emplace_back(tc); v.insert(v.end(), nrun-1, 0);} else { v.insert(v.end(), nrun-1, 0); v.emplace_back(tc);} i+=nrun; no_prev=0; continue; }
                v.emplace_back(1); ++i; ++pred; ++no_prev; continue;
            }
            uint64_t base = murmur64(x) & ht_mask; uint32_t ref_pos=0,b=0,f=0; if (!find_best(base, &text[i], (uint32_t)(text_len - i), no_prev, ref_pos, b, f)) { v.emplace_back(1); ++i; ++pred; ++no_prev; continue; }
            if (b) { for(uint32_t k=0;k<b;++k) v.pop_back(); pred -= b; i -= b; }
            uint32_t tc = cost_match(ref_pos - b, b+f, pred);
            if (prefix) { v.emplace_back(tc); v.insert(v.end(), b+f-1, 0); } else { v.insert(v.end(), b+f-1, 0); v.emplace_back(tc); }
            pred = ref_pos - b + b + f; i += b + f; no_prev = 0;
        }
        for (; i<text_len; ++i) v.emplace_back(1);
        return v;
    }

    inline uint64_t get_code_skip1(uint64_t code, const uint8_t* s) const { s += key_len - 1; if (*s>3) return ~0ull; code = ((code<<2) & key_mask) + *s; return code; }
    inline uint32_t get_nrun_len(const uint8_t* s, uint32_t max_len) const { if (s[0]!=N_CODE || s[1]!=N_CODE || s[2]!=N_CODE) return 0; uint32_t len=3; while (len<max_len && s[len]==N_CODE) ++len; return len; }
    inline bool find_best(uint64_t base, const uint8_t* s, uint32_t max_len, uint32_t no_prev, uint32_t& ref_pos, uint32_t& len_b, uint32_t& len_f) const {
        len_b=0; len_f=0; uint32_t min_upd=min_match_len; const uint8_t* ref_ptr = reference.data();
        for(uint32_t j=0;j<MAX_NO_TRIES;++j){ size_t idx=(base+j)&ht_mask; if (ht[idx]==UINT32_MAX) break; uint32_t h_pos = ht[idx]*HASHING_STEP; const uint8_t* p = ref_ptr + h_pos; uint32_t f=0; uint32_t maxf = std::min<uint32_t>(max_len, (uint32_t)(reference.size() - h_pos)); while (f<maxf && s[f]==p[f]) ++f; if (f>=key_len){ uint32_t b=0; uint32_t maxb = std::min<uint32_t>(no_prev, h_pos); while (b<maxb && s[-(int)b-1]==p[-(int)b-1]) ++b; if (b+f>min_upd){ len_b=b; len_f=f; ref_pos=h_pos; min_upd=b+f; } } }
        return len_b + len_f >= min_match_len;
    }
};
}

extern "C" size_t agc_cost_vector(int prefix, const uint8_t* ref, size_t ref_len, const uint8_t* text, size_t text_len, uint32_t min_match_len, uint32_t* out_costs) {
    CostEngine eng(min_match_len);
    eng.prepare(ref, ref_len);
    auto v = eng.get_costs(text, text_len, prefix != 0);
    if (out_costs) {
        for (size_t i=0;i<v.size();++i) out_costs[i] = v[i];
    }
    return v.size();
}

static inline void reverse_complement(const uint8_t* in, size_t n, std::vector<uint8_t>& out) {
    out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        uint8_t b = in[i];
        if (b < 4) out[n - 1 - i] = (uint8_t)(3 - b);
        else out[n - 1 - i] = b;
    }
}

// Compute best split position (index into text) matching C++ folding/orientation
// front_lt_mid and mid_lt_back represent (front_kmer < middle_kmer) and (middle_kmer < back_kmer)
extern "C" int agc_best_split(
    const uint8_t* left_ref, size_t left_len,
    const uint8_t* right_ref, size_t right_len,
    const uint8_t* text, size_t text_len,
    uint32_t min_match_len, uint32_t k,
    int front_lt_mid, int mid_lt_back,
    uint32_t* out_best_pos,
    uint32_t* out_seg2_start
) {
    if (!left_ref || !right_ref || !text || !out_best_pos || !out_seg2_start) return 0;

    CostEngine left_eng(min_match_len); left_eng.prepare(left_ref, left_len);
    CostEngine right_eng(min_match_len); right_eng.prepare(right_ref, right_len);

    std::vector<uint8_t> text_rc; reverse_complement(text, text_len, text_rc);

    // Left costs
    std::vector<uint32_t> v1;
    if (front_lt_mid) {
        v1 = left_eng.get_costs(text, text_len, /*prefix*/true);
        // forward cumulative
        uint32_t acc=0; for (auto& c : v1) { acc += c; c = acc; }
    } else {
        v1 = left_eng.get_costs(text_rc.data(), text_len, /*prefix*/false);
        // reverse costs for rc
        std::reverse(v1.begin(), v1.end());
        // forward cumulative
        uint32_t acc=0; for (auto& c : v1) { acc += c; c = acc; }
    }

    // Right costs (folded to cumulative aligned with text indices)
    std::vector<uint32_t> v2;
    if (mid_lt_back) {
        v2 = right_eng.get_costs(text, text_len, /*prefix*/false);
        // reverse cumulative (right -> left)
        uint32_t acc=0; for (auto it = v2.rbegin(); it != v2.rend(); ++it) { acc += *it; *it = acc; }
    } else {
        v2 = right_eng.get_costs(text_rc.data(), text_len, /*prefix*/true);
        // forward cumulative then reverse
        uint32_t acc=0; for (auto& c : v2) { acc += c; c = acc; }
        std::reverse(v2.begin(), v2.end());
    }

    if (v1.size() != v2.size() || v1.empty()) return 0;

    // Argmin (prefer leftmost minimal)
    uint32_t best_pos = 0; uint64_t best_sum = (uint64_t)UINT32_MAX + 1;
    for (size_t i=0;i<v1.size();++i) {
        uint64_t s = (uint64_t)v1[i] + (uint64_t)v2[i];
        if (s < best_sum) { best_sum = s; best_pos = (uint32_t)i; }
    }

    // Degenerate position rules
    if (best_pos < k + 1u) best_pos = 0;
    if ((size_t)best_pos + k + 1u > v1.size()) best_pos = (uint32_t)v1.size();

    *out_best_pos = best_pos;
    // Compute byte start for second segment with k-byte overlap: [seg2_start .. end], left ends at seg2_start + k
    // C++ uses k-byte overlap centered around the k-mer; for odd k, midpoint rounds up
    // Use half = (k+1)/2 so that seg2_start advances by one compared to plain k/2 for odd k.
    uint32_t half = (k + 1u) >> 1; // integer ceil(k/2)
    uint32_t seg2_start = 0;
    if (best_pos > half) seg2_start = best_pos - half;
    *out_seg2_start = seg2_start;
    return 1;
}
