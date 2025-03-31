// clang++ -O3 -std=c++17 -march=native standalone.cc
// GCC also compiles, but the performance is slightly worse
// Requires an x86 CPU with at least AVX2

// On Linux, make sure madvise is enabled to use hugepages
// (https://en.algorithmica.org/hpc/cpu-cache/paging/#changing-page-size)

// #pragma GCC optimize("O3")
#pragma once
#pragma GCC target("avx2,bmi")

#include <bits/stdc++.h>
#include <x86intrin.h>
#include <sys/mman.h>

template <typename K>
class SP_Tree {
    static const K B = 16, INF = std::numeric_limits<K>::max();
    K H, S;
    K *btree;
public:
    K N; // <- change these

    SP_Tree() = default;

    void build(std::vector<K>& q) {
        N = q.size();
        H = height(N);
        S = offset(N);
        prepare(q.data());
    }

    constexpr K blocks(K n) {
        return (n + B - 1) / B;
    }

    constexpr K prev_keys(K n) {
        return (blocks(n) + B) / (B + 1) * B;
    }

    constexpr K height(K n) {
        return (n <= B ? 1 : height(prev_keys(n)) + 1);
    }

    constexpr K offset(K h) {
        K k = 0, n = N;
        while (h--) {
            k += blocks(n) * B;
            n = prev_keys(n);
        }
        return k;
    }

    void permute(K *node) {
        const __m256i perm_mask = _mm256_set_epi32(3, 2, 1, 0, 7, 6, 5, 4);
        __m256i* middle = (__m256i*) (node + 4);
        __m256i x = _mm256_loadu_si256(middle);
        x = _mm256_permutevar8x32_epi32(x, perm_mask);
        _mm256_storeu_si256(middle, x);
    }

    void prepare(K *a) {
        const K P = 1 << 21, T = (4 * S + P - 1) / P * P;
        btree = (K*) std::aligned_alloc(P, T);
        #ifdef __linux__
        madvise(btree, T, MADV_HUGEPAGE);
        #endif
        // Windows:
        // btree = (K*) VirtualAlloc(NULL, T, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);

        for (K i = N; i < S; i++)
            btree[i] = INF;

        memcpy(btree, a, 4 * N);

        for (K h = 1; h < H; h++) {
            for (K i = 0; i < offset(h + 1) - offset(h); i++) {
                K k = i / B,
                    j = i - k * B;
                k = k * (B + 1) + j + 1;
                for (K l = 0; l < h - 1; l++)
                    k *= (B + 1);
                btree[offset(h) + i] = (k * B < N ? btree[k * B] : INF);
            }
        }

        for (K i = offset(1); i < S; i += B)
            permute(btree + i);
    }

    unsigned direct_rank(__m256i x, K* y) {
        __m256i a = _mm256_load_si256((__m256i*) y);
        __m256i b = _mm256_load_si256((__m256i*) (y + 8));

        __m256i ca = _mm256_cmpgt_epi32(a, x);
        __m256i cb = _mm256_cmpgt_epi32(b, x);

        K mb = _mm256_movemask_ps((__m256) cb);
        K ma = _mm256_movemask_ps((__m256) ca);

        unsigned mask = (1 << 16);
        mask |= mb << 8;
        mask |= ma;

        return __tzcnt_u32(mask);
    }

    unsigned permuted_rank(__m256i x, K* y) {
        __m256i a = _mm256_load_si256((__m256i*) y);
        __m256i b = _mm256_load_si256((__m256i*) (y + 8));

        __m256i ca = _mm256_cmpgt_epi32(a, x);
        __m256i cb = _mm256_cmpgt_epi32(b, x);

        __m256i c = _mm256_packs_epi32(ca, cb);
        unsigned mask = _mm256_movemask_epi8(c);

        return __tzcnt_u32(mask);
    }

    K lower_bound(K _x) {
        unsigned k = 0;
        __m256i x = _mm256_set1_epi32(_x - 1);
        for (K h = H - 1; h > 0; h--) {
            unsigned i = permuted_rank(x, btree + offset(h) + k);
            k = k * (B + 1) + (i << 3);
        }
        unsigned i = direct_rank(x, btree + k);
        return btree[k + i];
    }

    K lower_bound_index(K _x) {
        unsigned k = 0;
        __m256i x = _mm256_set1_epi32(_x - 1);
        for (K h = H - 1; h > 0; h--) {
            unsigned i = permuted_rank(x, btree + offset(h) + k);
            k = k * (B + 1) + (i << 3);
        }
        unsigned i = direct_rank(x, btree + k);
        return k + i;
    }

    K max_not_greater_index(K _x) {
        unsigned k = 0;
        __m256i x = _mm256_set1_epi32(_x - 1);
        for (K h = H - 1; h > 0; h--) {
            unsigned i = permuted_rank(x, btree + offset(h) + k);
            k = k * (B + 1) + (i << 3);
        }
        unsigned i = direct_rank(x, btree + k);
        K idx = k + i;
        if (idx >= N)
            idx = N - 1;
        else if (btree[idx] >_x && idx > 0)
            idx--;
        return idx;
    }

    // K baseline(K x) {
    //     return *std::lower_bound(a.begin(), a.end(), x);
    // }

    double timeit(std::string function_name, std::vector<K> q) {
        clock_t start = clock();

        K checksum1;

        if (function_name == "baseline") {
            for (K i = 0; i < q.size(); i++) {
                // checksum1 = baseline(q[i]);
            }
        } else if (function_name == "lower_bound") {
            for (K i = 0; i < q.size(); i++) {
                checksum1 = lower_bound(q[i]);
            }
        }
        double seconds = double(clock() - start) / CLOCKS_PER_SEC;
        printf("Checksum: %d\n", checksum1);

        return 1e9 * seconds / q.size();
    }

    void save_tree(const std::string& filename) const {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file for writing.");
        }

        size_t btree_size = S;
        file.write(reinterpret_cast<const char*>(&btree_size), sizeof(btree_size));
        file.write(reinterpret_cast<const char*>(btree), btree_size * sizeof(K));

        // Write other member variables
        file.write(reinterpret_cast<const char*>(&N), sizeof(N));
        file.write(reinterpret_cast<const char*>(&H), sizeof(H));
        file.write(reinterpret_cast<const char*>(&S), sizeof(S));

        if (!file.good()) {
            throw std::runtime_error("Error during serialization.");
        }
    }

    void  load_tree(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading.");
    }

    size_t btree_size;
    file.read(reinterpret_cast<char*>(&btree_size), sizeof(btree_size));
    const K P = 1 << 21;
    btree = (K*) std::aligned_alloc(P, btree_size * sizeof(K));
    file.read(reinterpret_cast<char*>(btree), btree_size * sizeof(K));

    // Read other member variables
    file.read(reinterpret_cast<char*>(&N), sizeof(N));
    file.read(reinterpret_cast<char*>(&H), sizeof(H));
    file.read(reinterpret_cast<char*>(&S), sizeof(S));

    if (!file.good()) {
        throw std::runtime_error("Error during deserialization.");
    }
}

};