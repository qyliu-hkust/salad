#include <immintrin.h>
#include <iostream>
#include <vp2intersect.h>
#include <vp2union.hpp>
// #include <hash_table.hpp>
#include <limits.h>
#include <optPFD_encode.hpp>
//
using namespace std;
//
//const int align_val = 64;
//const int key_nums = 16;
//
const int align_val = 32;
template <typename T>
T* aligned_new(uint64_t num_elements) {
    void* ptr = std::aligned_alloc(align_val, num_elements * sizeof(T));
    if (!ptr) throw std::bad_alloc();
    return static_cast<T*>(ptr);
}

uint32_t simsimd_intersect_u32_serial(uint32_t const *a, uint32_t const *b,uint32_t a_length, uint32_t b_length,uint32_t* intersect) {
    uint32_t* const intersect_start = intersect;

    const uint32_t* a_end = a + a_length;
    const uint32_t* b_end = b + b_length;

    while (a < a_end && b < b_end) {
        if (*a < *b) {
            ++a;
        } else if (*a > *b) {
            ++b;
        } else {
            // Found a match
            *intersect++ = *a;
            ++a;
            ++b;
        }
    }
    return intersect - intersect_start;
}

void simsimd_intersect_u32_ice(
    uint32_t const* shorter, uint32_t const* longer,
    uint32_t shorter_length, uint32_t longer_length,
    uint32_t& results) {

    uint32_t intersection_count = 0;
    uint32_t shorter_idx = 0, longer_idx = 0;
    uint32_t longer_load_size;
    __mmask16 longer_mask;

    while (shorter_idx < shorter_length && longer_idx < longer_length) {
        // Load `shorter_member` and broadcast it to shorter vector, load `longer_members_vec` from memory.
        uint32_t longer_remaining = longer_length - longer_idx;
        uint32_t shorter_member = shorter[shorter_idx];
        __m512i shorter_member_vec = _mm512_set1_epi32(*(int*)&shorter_member);
        __m512i longer_members_vec;
        if (longer_remaining < 16) {
            longer_load_size = longer_remaining;
            longer_mask = (__mmask16)_bzhi_u32(0xFFFF, longer_remaining);
        } else {
            longer_load_size = 16;
            longer_mask = 0xFFFF;
        }
        longer_members_vec = _mm512_maskz_loadu_epi32(longer_mask, (__m512i const*)(longer + longer_idx));

        // Compare `shorter_member` with each element in `longer_members_vec`,
        // and jump to the position of the match. There can be only one match at most!
        __mmask16 equal_mask = _mm512_mask_cmpeq_epu32_mask(longer_mask, shorter_member_vec, longer_members_vec);
        uint32_t equal_count = equal_mask != 0;
        intersection_count += equal_count;

        // When comparing a scalar against a sorted array, we can find three types of elements:
        // - entries that scalar is greater than,
        // - entries that scalar is equal to,
        // - entries that scalar is less than,
        // ... in that order! Any of them can be an empty set.
        __mmask16 greater_mask = _mm512_mask_cmplt_epu32_mask(longer_mask, longer_members_vec, shorter_member_vec);
        uint32_t greater_count = _mm_popcnt_u32(greater_mask);
        uint32_t smaller_exists = longer_load_size > greater_count - equal_count;

        // Advance the first array:
        // - to the next element, if a match was found,
        // - to the next element, if the current element is smaller than any elements in the second array.
        shorter_idx += equal_count | smaller_exists;
        // Advance the second array:
        // - to the next element after match, if a match was found,
        // - to the first element that is greater than the current element in the first array, if no match was found.
        longer_idx += greater_count + equal_count;

        // At any given cycle, take one entry from shorter array and compare it with multiple from the longer array.
        // For that, we need to swap the arrays if necessary.
        if ((shorter_length - shorter_idx) > (longer_length - longer_idx)) {
            uint32_t const* temp_array = shorter;
            shorter = longer, longer = temp_array;
            uint32_t temp_length = shorter_length;
            shorter_length = longer_length, longer_length = temp_length;
            uint32_t temp_idx = shorter_idx;
            shorter_idx = longer_idx, longer_idx = temp_idx;
        }
    }
    results = intersection_count;
}

uint32_t simsimd_intersect_u32(uint32_t const* a, uint32_t const* b, uint32_t a_length, uint32_t b_length, uint32_t *intersect) {
    uint32_t result = 0;
    uint32_t const* const a_end = a + a_length;
    uint32_t const* const b_end = b + b_length;
    uint32_t c = 0;
    union vec_t {
        __m512i zmm;
        uint32_t u32[16];
    } a_vec, b_vec;

    while (a + 16 < a_end && b + 16 < b_end) {
        a_vec.zmm = _mm512_loadu_si512((__m512i const*)a);
        b_vec.zmm = _mm512_loadu_si512((__m512i const*)b);

        // Intersecting registers with `_mm512_2intersect_epi16_mask` involves a lot of shuffling
        // and comparisons, so we want to avoid it if the slices don't overlap at all
        uint32_t a_min;
        uint32_t a_max = a_vec.u32[15];
        uint32_t b_min = b_vec.u32[0];
        uint32_t b_max = b_vec.u32[15];

        // If the slices don't overlap, advance the appropriate pointer
        while (a_max < b_min && a + 32 < a_end) {
            a += 16;
            a_vec.zmm = _mm512_loadu_si512((__m512i const*)a);
            a_max = a_vec.u32[15];
        }
        a_min = a_vec.u32[0];
        while (b_max < a_min && b + 32 < b_end) {
            b += 16;
            b_vec.zmm = _mm512_loadu_si512((__m512i const*)b);
            b_max = b_vec.u32[15];
        }
        b_min = b_vec.u32[0];

        // Now we are likely to have some overlap, so we can intersect the registers
        __mmask16 a_matches = _mm512_2intersect_epi32_mask(a_vec.zmm, b_vec.zmm);
        _mm512_mask_compressstoreu_epi32(intersect + c, a_matches, a_vec.zmm);
        // std::vector<uint32_t> intersection(16);
        // _mm512_mask_compressstoreu_epi32(intersection.data(), a_matches, a_vec);
        //
        // // 打印交集结果
        // for (int i = 0; i < 16; ++i) {
        //     std::cout << result[i] << " ";
        // }
        // std::cout << std::endl;


        c += _mm_popcnt_u32(a_matches); // The `_popcnt32` symbol isn't recognized by MSVC
        // cerr << _mm_popcnt_u32(a_matches) << endl;
        // c += __builtin_popcount(a_matches);

        // Determine the number of entries to skip in each array, by comparing
        // every element in the vector with the last (largest) element in the other array
        __m512i a_last_broadcasted = _mm512_set1_epi32(*(int const*)&a_max);
        __m512i b_last_broadcasted = _mm512_set1_epi32(*(int const*)&b_max);
        __mmask16 a_step_mask = _mm512_cmple_epu32_mask(a_vec.zmm, b_last_broadcasted);
        __mmask16 b_step_mask = _mm512_cmple_epu32_mask(b_vec.zmm, a_last_broadcasted);
        a += 16 - __lzcnt16((uint16_t)a_step_mask);
        b += 16 - __lzcnt16((uint16_t)b_step_mask);
    }

    // Handle the tail:
    c += simsimd_intersect_u32_serial(a, b, a_end - a, b_end - b, intersect + c);
    // result += c; // And merge it with the main body result
    return c;
}

alignas(align_val) int32_t *aa = aligned_new<int32_t>(32);
alignas(align_val) int32_t *ress = aligned_new<int32_t>(32);

int main() {
    // 生成两组随机数
    size_t n = 100000;
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n * 2];
    std::vector<uint32_t> intersect(n * 3);
    std::vector<uint32_t> intersect1(n * 3);
    for (int i = 0; i < n; i++) {
        a[i] = i;
        b[i] = i;
        b[i + n] = i + n;
    }

    sort(a, a + n);
    sort (b, b + n * 2);

    size_t result = 0;
    auto start = std::chrono::high_resolution_clock::now();
    // result = simsimd_intersect_u32_serial(a, b, n, n, intersect.data());
    // cerr << "Result 0: " << result << endl;
    result = union_u32_simd(a, b, n, n * 2, intersect.data());
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    cout << "Duration 1: " << duration.count() << " us" << endl;
    cerr << "Result 1: " << result << endl;
    //
    start = std::chrono::high_resolution_clock::now();
    result = union_scalar(a, b, n, n * 2, intersect.data());
    // result = union_u32_normal(a, b, n, n * 2, intersect.data());
    // result = simsimd_intersect_u32(a, b, n, n, intersect1.data());
    // simsimd_intersect_u32_ice(a, b, n, n, result);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    cout << "Duration 2: " << duration.count() << " us" << endl;
    cerr << "Result 2: " << result << endl;

    return 0;
}
