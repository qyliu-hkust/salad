#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"
#include "hash_table.hpp"
#include "perf_event.hpp"

namespace pgm_sequence {
    template <typename T>
    class HugePageAllocator {
        constexpr static size_t PAGE_SIZE = 2 * 1024 * 1024; // 2MB
    public:
        using value_type = T;

        // 构造函数和析构函数
        HugePageAllocator() = default;
        ~HugePageAllocator() = default;

        // 分配内存
        T* allocate(std::size_t n) {
            const size_t size = (n * sizeof(T) + PAGE_SIZE - 1) / PAGE_SIZE * PAGE_SIZE;
            // const size_t size = n * sizeof(T);
            void* ptr = mmap(nullptr, size,
                             PROT_WRITE | PROT_READ,
                             MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                             -1, 0);
            // cerr << size << " " << n * sizeof(T) / 4096 << " " << ptr << endl;
            if (ptr == MAP_FAILED) {
                std::cerr << "Failed to allocate huge pages: " << std::strerror(errno) << std::endl;
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        // 释放内存
        void deallocate(T* p, std::size_t n) noexcept {
            const size_t size = (n * sizeof(T) + PAGE_SIZE - 1) / PAGE_SIZE * PAGE_SIZE;
            munmap(p, size);
        }

        // 比较操作符（必须定义）
        bool operator == (const HugePageAllocator&) const { return true; }
        bool operator != (const HugePageAllocator&) const { return false; }
    };

    template <typename K>
    class pgm_enumerator{

        typedef int64_t Simd_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;

        public:

        struct segment{
            K first;
            Intercept_Value intercept; // 32 bits
            uint8_t slope_exponent;
            uint64_t slope_significand;
            Covered_Value covered; // covered
            inline segment(K first, Intercept_Value intercept, uint8_t slope_exponent, uint64_t slope_significand, Covered_Value covered) :
                    first(first), intercept(intercept), slope_exponent(slope_exponent), slope_significand(slope_significand), covered(covered) {}
        };

        // struct segment_query {
        //     K first;
        //     K value;
        //     Covered_Value covered; // covered
        //     inline segment_query(K first, K value, Covered_Value covered) :
        //             first(first), value(value), covered(covered) {}
        // };

        uint64_t n;

        std::vector<segment> segments;

        // std::vector<segment_query> segments_query;

        std::vector<Correction_Value> corrections_vector; // corrections for decode

        void load_copy(uint64_t data_size, std::vector<Correction_Value> corrections_vector) {
            this -> n = data_size;
            this -> corrections_vector = corrections_vector;
            this -> next_first_value = segments[1].intercept + corrections_vector[segments[1].first];
            this -> current_correction = segments[0].intercept;
            this -> current_value = segments[0].intercept + corrections_vector[0];
        }

        // used for Query Test
        K current_value = INT_MAX;
        K next_first_value = INT_MAX;
        Correction_Value current_correction = 0;
        Covered_Value current_pos = 0;
        Covered_Value current_segment = 0;
        std::vector<K, HugePageAllocator<K>> current_value_vector;
        // std::vector<K> current_value_vector;

        K docid() {
            return current_value;
        }

        void warm_up() {
            for (auto k = 0; k < 5; k++) {
                 for (auto i = 0;i < n; i++)
                     current_value_vector[i] = 0;
            }
        }

        void query_init(const std::string decode_type, const std::string query_type) {
            if (decode_type == "normal") {
                if (query_type == "intersection") {
                    current_pos = 0;
                    current_segment = 0;
                    total_skip = 0;
                    current_value_vector.resize(n);
                } else if (query_type == "union") {
                    current_pos = 0;
                    current_segment = 0;
                    current_value = INT_MAX - 1;
                    total_skip = 0;
                    current_value_vector.resize(n);
                }
            } else if (decode_type == "simd") {
                current_pos = 0;
                current_segment = 0;
                current_value_vector.resize(n);
                simd_init();
                vector<Correction_Value> ().swap(corrections_vector);
                vector<segment> ().swap(segments);
            }
        }

        void decode_query(const std::string decode_type) {
            if (decode_type == "normal") {
                normal_decode(current_value_vector.data());
            } else if (decode_type == "simd") {
                simd_decode_512i(current_value_vector.data());
            }
        }

        K next(const std::string &decode_type) {
            return current_value_vector[++current_pos];
        }

        long double total_skip = 0;

        std::vector<K> segment_decode(K idx) {
            std::vector<K> output;
            auto it = segments.begin() + idx;
            auto& s = *it;
            auto covered = s.covered;
            auto significand = s.slope_significand;
            auto exponent = s.slope_exponent;
            auto intercept = s.intercept;
            auto first = s.first;
            int32_t last_correction = intercept;
            for (Covered_Value j = first; j < first + covered; ++j) {
                last_correction = last_correction + corrections_vector[j];
                output.push_back(((j * significand) >> exponent) + last_correction);
            }
            return output;
        }

        K nextgeq(K posting_value) {
            if (current_value >= posting_value) {
                return current_value;
            }
            for (auto it = segments.begin() + current_segment; it < segments.end(); it++) {
                if (posting_value >= next_first_value) {
                    // if (current_pos <= it -> first)
                    //     total_skip += it ->covered;
                    // else
                    //     total_skip += it -> covered - (current_pos - it -> first);

                    auto it_next = std::next(std::next(it));
                    if (it_next < segments.end())
                        next_first_value = it_next -> intercept + corrections_vector[it_next -> first];
                    else
                        next_first_value = INT_MAX;
                    continue;
                }

                auto covered = it -> covered;
                auto first = it -> first;
                auto intercept = it -> intercept;
                auto slope_significand = it -> slope_significand;
                auto slope_exponent = it -> slope_exponent;

                if (current_segment != it - segments.begin()) {
                    current_pos = first;
                    current_segment = it - segments.begin();
                    current_correction = intercept;
                    auto it_next = std::next(it);
                    if (it_next < segments.end())
                        next_first_value = it_next -> intercept + corrections_vector[it_next -> first];
                    else
                        next_first_value = INT_MAX;
                }

                for (Covered_Value j = current_pos - first; j < covered; ++j) {
                    current_correction = current_correction + corrections_vector[j + first];
                    current_value = ((j * slope_significand) >> slope_exponent) + current_correction;
                    if (current_value >= posting_value) {
                        current_pos = j + first + 1;
                        return current_value;
                    }
                }
            }
            current_pos = n;
            return INT_MAX;
        }

        void normal_decode(K* output) {
            uint32_t pointer = 0;
            auto start = std::chrono::high_resolution_clock::now();
            for (auto it = segments.begin(); it != segments.end(); ++it) {
                auto intercept = it -> intercept;
                auto exponent = it -> slope_exponent;
                auto significand = it -> slope_significand;
                auto covered = it -> covered;
                int32_t last_correction = intercept;
                for (int32_t j = 0; j < covered; ++j) {
                    last_correction = last_correction + corrections_vector[pointer++];
                    *output++ = ((j * significand) >> exponent) + last_correction;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration = duration.count();
        }

        // used for SIMD
        constexpr static K key_nums = 8;
        constexpr static K align_val = 64; // for avx512
        std::vector<segment> segments_sort; // resorted segments
        std::vector<Simd_Value*> slope_significand_simd;
        std::vector<Simd_Value*> slope_exponent_simd;
        std::vector<Intercept_Value*> intercept_simd;
        std::vector<Correction_Value*> corrections_simd;
        // std::vector<Correction_Value*, HugePageAllocator<Correction_Value*>> corrections_simd;
        std::vector<Correction_Value> corrections_vector_residual;
        std::vector<Covered_Value*> first_simd;
        std::vector<Covered_Value*> covered_simd;
        std::vector<Covered_Value> cover_length;
        HashTable<Covered_Value, Covered_Value> decode_result_map;

        uint64_t total_calculated = 0;
        uint64_t total_calculated_add = 0;
        uint64_t conversion_time = 0;
        uint64_t total_duration = 0;

        uint64_t idx = 0;

        template <typename T>
        T* aligned_new(uint64_t num_elements) {
            void* ptr = std::aligned_alloc(align_val, num_elements * sizeof(T));
            if (!ptr) throw std::bad_alloc();
            return static_cast<T*>(ptr);
        }

        template <typename T>
        static void aligned_delete(T* ptr) {
            std::free(ptr);
        }

        constexpr static size_t HUGE_PAGE_SIZE = 2 * 1024 * 1024;

        template <typename T>
        T* aligned_new_huge(uint64_t num_elements) {
            const size_t size = (num_elements * sizeof(T) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE;
            void* ptr = mmap(nullptr, size,
                             PROT_WRITE | PROT_READ,
                             MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                             -1, 0);
            if (ptr == MAP_FAILED) {
                cerr << "Failed to allocate huge pages: " << std::strerror(errno) << " " << size << " " << num_elements << std::endl;
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        template <typename T>
        static void aligned_delete_huge(T* p, std::size_t num_elements) noexcept {
            const size_t size = (num_elements * sizeof(T) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE;
            munmap(p, size);
        }

        void free_memory(string simd_type = "simd") {
            if (simd_type == "simd") {
                for (auto i = 0; i < slope_significand_simd.size(); i++) {
                    aligned_delete(slope_significand_simd[i]);
                    aligned_delete(slope_exponent_simd[i]);
                    aligned_delete(intercept_simd[i]);
                    aligned_delete(corrections_simd[i]);
                    // aligned_delete_huge(corrections_simd[i], cover_length[i]);
                    aligned_delete(first_simd[i]);
                    aligned_delete(covered_simd[i]);
                }

                vector<Covered_Value> ().swap(cover_length);
                // vector<Correction_Value*, HugePageAllocator<Correction_Value*>> ().swap(corrections_simd);
                vector<Correction_Value*> ().swap(corrections_simd);
                vector<Correction_Value> ().swap(corrections_vector_residual);
                vector<Simd_Value*> ().swap(slope_significand_simd);
                vector<Simd_Value*> ().swap(slope_exponent_simd);
                vector<Intercept_Value*> ().swap(intercept_simd);
                vector<Covered_Value*> ().swap(first_simd);
                vector<Covered_Value*> ().swap(covered_simd);
                vector<segment> ().swap(segments_sort);
                // vector<K> ().swap(current_value_vector);
                vector<K, HugePageAllocator<K>> ().swap(current_value_vector);
            }
        }

        // our SIMD
        void simd_init() {
            segments_sort = segments;
            std::sort(segments_sort.begin(), segments_sort.end(), [](const segment &a, const segment &b) {return a.covered > b.covered;}); // sorted by covered
            Covered_Value min_cover = INT_MAX;
            Covered_Value max_min_covered = 0;
            idx = 0;

            for (auto it = segments_sort.begin(); it + key_nums < segments_sort.end(); it = it + key_nums) {
                alignas(align_val) Simd_Value *slope_significand_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Simd_Value *slope_exponent_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Intercept_Value *intercept_simd_tmp = aligned_new<Intercept_Value>(key_nums);
                alignas(align_val) Covered_Value *covered_tmp = aligned_new<Covered_Value>(key_nums);
                alignas(align_val) Covered_Value *first_tmp = aligned_new<Covered_Value>(key_nums);

                std::vector<segment> simd_segments(it, it + key_nums);
                if (simd_segments.back().covered < 32) {
                    break;
                }
                std::sort(simd_segments.begin(), simd_segments.end(), [](const segment &a, const segment &b) {return a.first < b.first;}); // part sorted
                idx += key_nums;

                for (auto i = 0, its = simd_segments.begin(); its != simd_segments.end(); ++its, ++i) {
                    auto covered = its -> covered;
                    slope_significand_simd_tmp[i] = static_cast<Simd_Value>(its -> slope_significand);
                    slope_exponent_simd_tmp[i] = static_cast<Simd_Value>(its -> slope_exponent);
                    intercept_simd_tmp[i] = static_cast<Intercept_Value>(its -> intercept);
                    covered_tmp[i] = static_cast<Covered_Value> (covered);
                    first_tmp[i] = static_cast<Covered_Value> (its -> first);
                    min_cover = min_cover < covered ? min_cover : covered;
                    // max_cover = max_cover > covered ? max_cover : covered;
                }

                slope_significand_simd.emplace_back(slope_significand_simd_tmp);
                slope_exponent_simd.emplace_back(slope_exponent_simd_tmp);
                intercept_simd.emplace_back(intercept_simd_tmp);
                first_simd.emplace_back(first_tmp);
                covered_simd.emplace_back(covered_tmp);
                max_min_covered = min_cover - min_cover % 2;
                cover_length.emplace_back(max_min_covered);
            }
            create_corrections();
            create_corrections_residual();
        }

        void create_corrections() {
             corrections_simd.resize(cover_length.size(), 0);
            for (int i = 0;i < cover_length.size(); i++) {
                Covered_Value cover_length_tmp = cover_length[i];
                Covered_Value *first_tmp = first_simd[i];
                alignas(align_val) Correction_Value *corrections_simd_tmp = aligned_new<Correction_Value>(cover_length_tmp * key_nums);
                // Correction_Value *corrections_simd_tmp = aligned_new_huge<Correction_Value>(cover_length_tmp * key_nums);
                for (Covered_Value j = 0; j < cover_length_tmp; j++) {
                    for (Covered_Value k = 0; k < key_nums; k++) {
                        corrections_simd_tmp[j * key_nums + k] = corrections_vector[first_tmp[k] + j];
                    }
                }
                corrections_simd[i] = corrections_simd_tmp;
            }
        }

        Covered_Value count_residual_size() {
            Covered_Value corrections_vector_residual_size = 0;
            for(int i = 0; i < cover_length.size(); i++) {
                Covered_Value *cover_tmp = covered_simd[i];
                Covered_Value cover_length_tmp = cover_length[i];
                for (Covered_Value k = 0; k < key_nums; k++) {
                    auto covered = cover_tmp[k];
                    if (cover_length_tmp < covered)
                        corrections_vector_residual_size += covered - cover_length_tmp;
                }
            }
            for(auto it = segments_sort.begin() + idx; it < segments_sort.end(); it++)
                corrections_vector_residual_size += it -> covered;
            return corrections_vector_residual_size;
        }

        void create_corrections_residual() {
             corrections_vector_residual.resize(count_residual_size());
             Covered_Value correct_pointers = 0;
             for(int i = 0; i < cover_length.size(); i++) {
                 Covered_Value *first_tmp = first_simd[i];
                 Covered_Value *cover_tmp = covered_simd[i];
                 Covered_Value cover_length_tmp = cover_length[i];
                 for (K k = 0; k < key_nums; k++) {
                     if (cover_length_tmp < cover_tmp[k]) {
                         for (Covered_Value pos = cover_length_tmp; pos < cover_tmp[k]; pos++){
                             corrections_vector_residual[correct_pointers++] = corrections_vector[first_tmp[k] + pos];
                         }
                     }
                 }
             }

             for(auto it = segments_sort.begin() + idx; it < segments_sort.end(); it++) {
                 for (Covered_Value pos = it -> first; pos < (it -> first + it -> covered); pos++){
                     corrections_vector_residual[correct_pointers++] = corrections_vector[pos];
                 }
             }

            if (count_residual_size() != correct_pointers) {
                cerr << "Error: create_corrections_residual" << count_residual_size() << " " << correct_pointers << endl;
            }
        }

        void simd_decode_512i(K* output){
            Correction_Value* correct_pointers = corrections_vector_residual.data();
            total_calculated = 0;
            total_calculated_add = 0;
            total_duration = 0;
            // std::array<K*, key_nums> first_pointer;
            __m512i rerange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            alignas(align_val) Correction_Value *last_correction32 = aligned_new<Correction_Value>(8);

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length.size(); i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them
                Covered_Value  *first_tmp = first_simd[i];
                K *p0 = output + first_tmp[0];
                K *p1 = output + first_tmp[1];
                K *p2 = output + first_tmp[2];
                K *p3 = output + first_tmp[3];
                K *p4 = output + first_tmp[4];
                K *p5 = output + first_tmp[5];
                K *p6 = output + first_tmp[6];
                K *p7 = output + first_tmp[7];

                // 0
                __m256i intercept_v = _mm256_load_epi32(intercept_simd[i]);
                Correction_Value *corrections_p = corrections_simd[i];
                __m256i corrections_v = _mm256_load_epi32(corrections_p);
                corrections_p += 8;
                intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                __m512i result_v = _mm512_castsi256_si512(intercept_v); // lower 256 bits

                const Simd_Value *slope_significand_p = slope_significand_simd[i];
                const Simd_Value *slope_exponent_p = slope_exponent_simd[i];
                const __m512i slope_significand_v_tmp1 = _mm512_load_epi64(slope_significand_p);
                __m512i slope_significand_v1 = slope_significand_v_tmp1;
                __m512i slope_exponent_v1 = _mm512_load_epi64(slope_exponent_p);

                corrections_v = _mm256_load_epi32(corrections_p);
                corrections_p += 8;
                intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                __m512i slope_correct_v = _mm512_srlv_epi64(slope_significand_v1, slope_exponent_v1);
                __m256i int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                int32_v = _mm256_add_epi32(int32_v, intercept_v);
                result_v = _mm512_inserti64x4(result_v, int32_v, 1); // upper 256 bits
                result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                __m128i t0 = _mm512_extracti32x4_epi32(result_v, 0);
                __m128i t1 = _mm512_extracti32x4_epi32(result_v, 1);
                __m128i t2 = _mm512_extracti32x4_epi32(result_v, 2);
                __m128i t3 = _mm512_extracti32x4_epi32(result_v, 3);

                *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // a0, a1
                *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // a0, a1
                *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // a0, a1
                *((int64_t*)(p4)) = _mm_extract_epi64(t2, 0);  // a0, a1
                *((int64_t*)(p5)) = _mm_extract_epi64(t2, 1);  // a0, a1
                *((int64_t*)(p6)) = _mm_extract_epi64(t3, 0);  // a0, a1
                *((int64_t*)(p7)) = _mm_extract_epi64(t3, 1);  // a0, a1

                p0 = p0 + 2;
                p1 = p1 + 2;
                p2 = p2 + 2;
                p3 = p3 + 2;
                p4 = p4 + 2;
                p5 = p5 + 2;
                p6 = p6 + 2;
                p7 = p7 + 2;

                const Covered_Value cover_length_tmp = cover_length[i];
                total_calculated += cover_length_tmp;

                for (Covered_Value j = 2; j < cover_length_tmp; j += 2) {
                    corrections_v = _mm256_load_epi32(corrections_p);
                    corrections_p += 8;
                    intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                    slope_significand_v1 = _mm512_add_epi64(slope_significand_v1, slope_significand_v_tmp1);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v1, slope_exponent_v1);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    int32_v = _mm256_add_epi32(int32_v, intercept_v);
                    result_v = _mm512_castsi256_si512(int32_v); // lower 256 bits

                    corrections_v = _mm256_load_epi32(corrections_p);
                    corrections_p += 8;
                    intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                    slope_significand_v1 = _mm512_add_epi64(slope_significand_v1, slope_significand_v_tmp1);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v1, slope_exponent_v1);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    int32_v = _mm256_add_epi32(int32_v, intercept_v);
                    result_v = _mm512_inserti64x4(result_v, int32_v, 1); // upper 256 bits

                    result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                    t0 = _mm512_extracti32x4_epi32(result_v, 0);
                    t1 = _mm512_extracti32x4_epi32(result_v, 1);
                    t2 = _mm512_extracti32x4_epi32(result_v, 2);
                    t3 = _mm512_extracti32x4_epi32(result_v, 3);

                    *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                    *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // a0, a1
                    *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // a0, a1
                    *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // a0, a1
                    *((int64_t*)(p4)) = _mm_extract_epi64(t2, 0);  // a0, a1
                    *((int64_t*)(p5)) = _mm_extract_epi64(t2, 1);  // a0, a1
                    *((int64_t*)(p6)) = _mm_extract_epi64(t3, 0);  // a0, a1
                    *((int64_t*)(p7)) = _mm_extract_epi64(t3, 1);  // a0, a1

                    p0 = p0 + 2;
                    p1 = p1 + 2;
                    p2 = p2 + 2;
                    p3 = p3 + 2;
                    p4 = p4 + 2;
                    p5 = p5 + 2;
                    p6 = p6 + 2;
                    p7 = p7 + 2;
                }

                _mm256_store_epi32(last_correction32, intercept_v);

                const Covered_Value *cover_tmp = covered_simd[i];
                for (int k = 0; k < 8; k++) {
                    const Covered_Value covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        const Simd_Value slope_significand = slope_significand_p[k];
                        const Simd_Value slope_exponent = slope_exponent_p[k];
                        int32_t last_correction = last_correction32[k];
                        K* p_tmp = output + first_tmp[k] + cover_length_tmp;
                        for (Correction_Value pos = cover_length_tmp; pos < covered; pos++) {
                            last_correction = last_correction + *correct_pointers++;
                            *p_tmp++ = ((slope_significand * pos) >> slope_exponent) + last_correction;
                        }
                    }
                }
            }


            const auto end_iter = segments_sort.end();
            for (auto it = segments_sort.begin() + idx; it < end_iter; ++it) {
                const auto& seg = *it;
                int32_t last_correction = seg.intercept;
                K* dst = output + seg.first;
                for (Covered_Value pos = 0; pos < seg.covered; ++pos) {
                    last_correction = last_correction + *correct_pointers++;
                    *dst++ = ((seg.slope_significand * pos) >> seg.slope_exponent) + last_correction;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration = duration.count();
            total_calculated = total_calculated * 8;
        }

        // for compress residual
        inline Correction_Value liner_spline(Correction_Value x1, Correction_Value x2) {
            return (x1 + x2) >> 1;
        }

        inline Correction_Value compare_max(Correction_Value x1, Correction_Value x2) {
            return x1 > x2 ? x1 : x2;
        }

        inline Correction_Value compare_min(Correction_Value x1, Correction_Value x2) {
            return x1 < x2 ? x1 : x2;
        }

        std::vector<Correction_Value> residuals_compress; // corrections for decode
        Correction_Value residual_max, residual_min, residual_5, residual_10, residual_20, residual_25, residual_50, residual_75, residual_80, residual_90, residual_95;
        uint64_t Epsilon_Data = 0;

        void spline_compress() {
            for (auto i = 0; i < this -> segments.size(); i++) {
                auto seg = this -> segments[i];
                for (auto j = seg.first + 1; j < seg.first + seg.covered - 1; j+=2) {
                    residuals_compress.push_back(liner_spline(this -> corrections_vector[j - 1], this -> corrections_vector[j + 1]) - corrections_vector[j]);
                }
            }
            std::sort(residuals_compress.begin(), residuals_compress.end());
            residual_min = residuals_compress[0];
            residual_25 = residuals_compress[residuals_compress.size() / 4 * 1];
            residual_50 = residuals_compress[residuals_compress.size() / 4 * 2];
            residual_75 = residuals_compress[residuals_compress.size() / 4 * 3];
            residual_max = residuals_compress[residuals_compress.size() - 1];
        }

        int over_data = 0;
        int over_num = 0;
        int last_max = 0;
        int max_distance = -10;
        int compress_bit = 3;
        uint64_t total_residual_bit_size_distance = 0;
        uint64_t total_residual_bit_size_flag = 0;

        void second_difference_compress() {
            over_data = Epsilon_Data >> compress_bit;
            for (auto i = 0; i < this -> segments.size() - 1; i++) {
                auto seg = this -> segments[i];
                residuals_compress.push_back(abs(corrections_vector[seg.first]));

                last_max = seg.first;
                for (int j = seg.first + 1; j < seg.first + seg.covered && j < n; j++) {
                    int diff = abs(corrections_vector[j] - corrections_vector[j - 1]);
                    residuals_compress.push_back(diff);
                    if (diff <= over_data) {
                        total_residual_bit_size_distance += BIT_WIDTH(over_data);
                    }
                    else {
                        over_num++;
                        if (j - last_max - 1 > max_distance) {
                            max_distance = j - last_max - 1;
                            last_max = j;
                        }
                    }
                    // residuals_compress.push_back(abs(corrections_vector[j]));
                }
            }
            total_residual_bit_size_flag = total_residual_bit_size_distance + over_num * (BIT_WIDTH(Epsilon_Data) + 1 + BIT_WIDTH(over_data));
            total_residual_bit_size_distance = total_residual_bit_size_distance + over_num * (BIT_WIDTH(Epsilon_Data) + 1 + BIT_WIDTH(max_distance));

            std::sort(residuals_compress.begin(), residuals_compress.end());
            residual_min = residuals_compress[0];
            residual_5 = residuals_compress[residuals_compress.size() / 20 * 1];
            residual_10 = residuals_compress[residuals_compress.size() / 10 * 1];
            residual_20 = residuals_compress[residuals_compress.size() / 5 * 1];
            residual_25 = residuals_compress[residuals_compress.size() / 4 * 1];
            residual_50 = residuals_compress[residuals_compress.size() / 4 * 2];
            residual_75 = residuals_compress[residuals_compress.size() / 4 * 3];
            residual_80 = residuals_compress[residuals_compress.size() / 5 * 4];
            residual_90 = residuals_compress[residuals_compress.size() / 10 * 9];
            residual_95 = residuals_compress[residuals_compress.size() / 20 * 19];
            residual_max = residuals_compress[residuals_compress.size() - 1];
        }
    };
}

