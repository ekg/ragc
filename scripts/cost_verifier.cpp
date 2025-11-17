// Minimal C++ coding cost verifier mirroring C++ AGC CLZDiffBase::GetCodingCostVector
#include <bits/stdc++.h>
using namespace std;

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
    vector<uint8_t> reference; // padded with INVALID
    vector<uint32_t> ht; // linear probing table (pos/HASHING_STEP), u32::max => empty
    uint64_t ht_mask{};

    explicit CostEngine(uint32_t mml): min_match_len(mml) {
        key_len = min_match_len - HASHING_STEP + 1u;
        key_mask = ~0ull >> (64 - 2 * key_len);
    }

    void prepare(const vector<uint8_t>& ref) {
        reference = ref;
        reference.resize(reference.size() + key_len, INVALID);
        build_index();
    }

    void build_index() {
        // Count valid positions for sparse mode
        uint64_t sz = 0; uint32_t no_prev_valid = 0; uint32_t cnt_mod = 0; uint32_t key_len_mod = key_len % HASHING_STEP;
        for (auto c: reference) {
            if (c < 4) ++no_prev_valid; else no_prev_valid = 0;
            if (++cnt_mod == HASHING_STEP) cnt_mod = 0;
            if (cnt_mod == key_len_mod && no_prev_valid >= key_len) ++sz;
        }
        uint64_t ht_size = (uint64_t)(sz / 0.7);
        while (ht_size & (ht_size - 1)) ht_size &= ht_size - 1; ht_size <<= 1; if (ht_size < 8) ht_size = 8;
        ht_mask = ht_size - 1;
        ht.assign(ht_size, numeric_limits<uint32_t>::max());

        // Insert with linear probing
        size_t ref_len = reference.size();
        for (size_t i = 0; i + key_len < ref_len; i += HASHING_STEP) {
            uint64_t code = get_code(&reference[i]); if (code == ~0ull) continue;
            uint64_t base = murmur64(code) & ht_mask;
            for (uint32_t j = 0; j < MAX_NO_TRIES; ++j) {
                size_t idx = (base + j) & ht_mask;
                if (ht[idx] == numeric_limits<uint32_t>::max()) { ht[idx] = (uint32_t)(i / HASHING_STEP); break; }
            }
        }
    }

    uint64_t get_code(const uint8_t* s) const {
        uint64_t x = 0; uint32_t i = key_len % 4;
        switch (i) {
            case 3: if (s[0] > 3) return ~0ull; x = (x<<2) + s[0]; ++s; [[fallthrough]];
            case 2: if (s[0] > 3) return ~0ull; x = (x<<2) + s[0]; ++s; [[fallthrough]];
            case 1: if (s[0] > 3) return ~0ull; x = (x<<2) + s[0]; ++s;
        }
        for (; i < key_len; ) {
            if (*s > 3) return ~0ull; x = (x<<2) + *s++; ++i;
            if (*s > 3) return ~0ull; x = (x<<2) + *s++; ++i;
            if (*s > 3) return ~0ull; x = (x<<2) + *s++; ++i;
            if (*s > 3) return ~0ull; x = (x<<2) + *s++; ++i;
        }
        return x;
    }

    uint64_t get_code_skip1(uint64_t code, const uint8_t* s) const {
        s += key_len - 1; if (*s > 3) return ~0ull; code = ((code<<2) & key_mask) + *s; return code;
    }

    uint32_t get_Nrun_len(const uint8_t* s, uint32_t max_len) const {
        if (s[0] != N_CODE || s[1] != N_CODE || s[2] != N_CODE) return 0;
        uint32_t len = 3; while (len < max_len && s[len] == N_CODE) ++len; return len;
    }

    static inline uint32_t int_len(uint32_t x) {
        if (x < 10) return 1; if (x < 100) return 2; if (x < 1000) return 3; if (x < 10000) return 4;
        if (x < 100000) return 5; if (x < 1000000) return 6; if (x < 10000000) return 7; if (x < 100000000) return 8;
        if (x < 1000000000) return 9; return 10;
    }

    uint32_t coding_cost_match(uint32_t ref_pos, uint32_t len, uint32_t pred_pos) const {
        int dif = (int)ref_pos - (int)pred_pos; uint32_t r = dif >= 0 ? int_len((uint32_t)dif) : int_len((uint32_t)(-dif)) + 1;
        r += int_len(len - min_match_len) + 2; return r;
    }
    uint32_t coding_cost_Nrun(uint32_t len) const { return 2 + int_len(len - MIN_NRUN_LEN); }

    bool find_best_match(uint64_t base, const uint8_t* s, uint32_t max_len, uint32_t no_prev_literals,
                         uint32_t& ref_pos, uint32_t& len_bck, uint32_t& len_fwd) const {
        len_bck = len_fwd = 0; uint32_t min_to_update = min_match_len; const uint8_t* ref_ptr = reference.data();
        for (uint32_t j = 0; j < MAX_NO_TRIES; ++j) {
            size_t idx = (base + j) & ht_mask; if (ht[idx] == numeric_limits<uint32_t>::max()) break;
            uint32_t h_pos = ht[idx] * HASHING_STEP; const uint8_t* p = ref_ptr + h_pos;
            // forward
            uint32_t f_len = 0; uint32_t maxf = min<uint32_t>((uint32_t)max_len, (uint32_t)(reference.size() - h_pos));
            while (f_len < maxf && s[f_len] == p[f_len]) ++f_len;
            if (f_len >= key_len) {
                uint32_t b_len = 0; uint32_t maxb = min(no_prev_literals, h_pos);
                while (b_len < maxb && s[-(int)b_len - 1] == p[-(int)b_len - 1]) ++b_len;
                if (b_len + f_len > min_to_update) { len_bck = b_len; len_fwd = f_len; ref_pos = h_pos; min_to_update = b_len + f_len; }
            }
        }
        return len_bck + len_fwd >= min_match_len;
    }

    vector<uint32_t> get_costs(const vector<uint8_t>& text, bool prefix_costs) const {
        vector<uint32_t> v; v.reserve(text.size());
        uint32_t i = 0, pred_pos = 0, no_prev_literals = 0; uint64_t x_prev = ~0ull;
        while (i + key_len < text.size()) {
            uint64_t x = (x_prev != ~0ull && no_prev_literals > 0) ? get_code_skip1(x_prev, &text[i]) : get_code(&text[i]);
            x_prev = x;
            if (x == ~0ull) {
                uint32_t nrun = get_Nrun_len(&text[i], (uint32_t)(text.size() - i));
                if (nrun >= MIN_NRUN_LEN) {
                    auto tc = coding_cost_Nrun(nrun);
                    if (prefix_costs) { v.emplace_back(tc); v.insert(v.end(), nrun - 1, 0); }
                    else { v.insert(v.end(), nrun - 1, 0); v.emplace_back(tc); }
                    i += nrun; no_prev_literals = 0; continue;
                } else { v.emplace_back(1); ++i; ++pred_pos; ++no_prev_literals; continue; }
            }
            uint64_t base = murmur64(x) & ht_mask; uint32_t ref_pos=0, len_bck=0, len_fwd=0;
            if (!find_best_match(base, &text[i], (uint32_t)(text.size() - i), no_prev_literals, ref_pos, len_bck, len_fwd)) {
                v.emplace_back(1); ++i; ++pred_pos; ++no_prev_literals; continue;
            } else {
                if (len_bck) { for (uint32_t k=0;k<len_bck;++k) v.pop_back(); pred_pos -= len_bck; i -= len_bck; }
                uint32_t tc = coding_cost_match(ref_pos - len_bck, len_bck + len_fwd, pred_pos);
                if (prefix_costs) { v.emplace_back(tc); v.insert(v.end(), len_bck + len_fwd - 1, 0); }
                else { v.insert(v.end(), len_bck + len_fwd - 1, 0); v.emplace_back(tc); }
                pred_pos = ref_pos - len_bck + len_bck + len_fwd; i += len_bck + len_fwd; no_prev_literals = 0;
            }
        }
        for (; i < text.size(); ++i) v.emplace_back(1);
        return v;
    }
};

static vector<uint8_t> read_all(const string& path) {
    ifstream f(path, ios::binary); vector<uint8_t> v((istreambuf_iterator<char>(f)), {}); return v;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: cost_verifier <prefix:0|1> <ref.bin> <text.bin>\n";
        return 1;
    }
    bool prefix = string(argv[1]) == "1";
    auto ref = read_all(argv[2]); auto txt = read_all(argv[3]);
    // Default min_match_len=20 to match CLI default; override via env if needed
    uint32_t mml = 20; if (const char* e = getenv("RAGC_MML")) mml = (uint32_t)stoi(e);
    CostEngine eng(mml); eng.prepare(ref);
    auto costs = eng.get_costs(txt, prefix);
    // Print summary
    cout << "LEN " << costs.size() << "\n";
    size_t best = 0; uint64_t best_sum = UINT64_MAX; // Note: split uses cumulative sums afterwards
    // For comparison purposes, just print first/last few entries
    cout << "HEAD "; for (size_t i=0; i<min<size_t>(costs.size(), 20); ++i) cout << costs[i] << (i+1<20? ',':'\n');
    cout << "TAIL "; for (size_t i=costs.size() - min<size_t>(costs.size(), 20); i<costs.size(); ++i) cout << costs[i] << (i+1<costs.size()? ',':'\n');
    // Not computing argmin of combined costs here; ragc combines left/right cumulative costs.
    return 0;
}

