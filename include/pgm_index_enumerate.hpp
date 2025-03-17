#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"
#include "hash_table.hpp"

namespace pgm_sequence {
    template <typename K>
    class pgm_enumerator{

        typedef int64_t Simd_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;
        public:

        // pgm_enumerator() = default;

        struct segment{
            K first;
            Intercept_Value intercept; // 32 bits
            uint8_t slope_exponent;
            uint64_t slope_significand;
            Covered_Value covered; // covered
            inline segment(K first, Intercept_Value intercept, uint8_t slope_exponent, uint64_t slope_significand, Covered_Value covered) :
                    first(first), intercept(intercept), slope_exponent(slope_exponent), slope_significand(slope_significand), covered(covered) {}
        };

        struct segment_query {
            K first;
            K value;
            Covered_Value covered; // covered
            inline segment_query(K first, K value, Covered_Value covered) :
                    first(first), value(value), covered(covered) {}
        };

        uint64_t n;

        std::vector<segment> segments;

        std::vector<segment_query> segments_query;

        std::vector<Correction_Value> corrections_vector; // corrections for decode

        void load_copy(uint64_t data_size, std::vector<Correction_Value> corrections_vector) {
            this -> n = data_size;
            this -> corrections_vector = corrections_vector;
            this -> next_first_value = segments[1].intercept + corrections_vector[segments[1].first];
        }

        // used for Query Test
        K current_value = INT_MAX;
        K next_first_value = INT_MAX;
        Covered_Value current_pos = 0;
        Covered_Value current_segment = 0;
        std::vector<K> current_value_vector;

        K docid() {
            return current_value;
        }

        void warm_up() {
            for (auto k = 0; k < 5; k++) {
                // if (!corrections_vector.empty()) {
                //     for (auto i = 0; i < n; i++) {
                //         auto tmp = corrections_vector[i];
                //         asm volatile("" : : "r,m"(tmp) : "memory");
                //     }
                // }
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
                normal_decode_query();
            } else if (decode_type == "simd") {
                simd_decode_512i_query();
            }
        }

        K next(const std::string &decode_type) {
            return current_value_vector[++current_pos];
//             if (decode_type == "normal") {
//                 return current_value_vector[++current_pos];
//                 // next_normal();
//             } else if (decode_type == "simd") {
//                 return current_value_vector[++current_pos];
// //                next_simd();
//             } else if (decode_type == "simd_simple") {
//                 // return nextgeq_simd_simple(posting_value);
//             } else {
//                 throw std::invalid_argument("The decode type " + decode_type + " is not supported.");
//             }
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
            for (Covered_Value j = 0; j < covered; ++j) {
                output.push_back(((j * significand) >> exponent) + intercept + corrections_vector[j +first]);
            }
            return output;
        }

        int find_relative_segment_num(std::vector<K> s0) {
            K smin_value = s0[0];
            K smax_value = s0[s0.size() - 1];
            int relative_sement_num = 0;
            for (auto it = segments.begin() + current_segment; it < std::prev(segments.end()); it++) {
                K seg_min = it -> intercept + corrections_vector[it -> first];
                K seg_max = 0;
                if (it < std::prev(std::prev(segments.end())))
                    seg_max = std::next(it) -> intercept + corrections_vector[std::next(it) -> first];
                else
                    seg_max = ((((it -> covered) - 1) * it -> slope_significand) >> it -> slope_exponent) + it -> intercept + corrections_vector[it -> covered - 1 + it -> first];
                if (seg_max < smin_value || smax_value < seg_min) {
                    continue;
                } else {
                    relative_sement_num++;
                }
            }
            return relative_sement_num;
        }

        K nextgeq(K posting_value) {

            // for (auto i = current_pos; i < n; i++) {
            //     if (current_value_vector[i] >= posting_value) {
            //         current_pos = i;
            //         return current_value_vector[i];
            //     }
            // }
            //
            // current_pos = n;
            // return INT_MAX;

            for (auto it = segments.begin() + current_segment; it < std::prev(segments.end()); it++) {
                if (posting_value >= next_first_value) {
                    // if (current_pos <= it -> first)
                    //     total_skip += it ->covered;
                    // else
                    //     total_skip += it -> covered - (current_pos - it -> first);

                    auto it_next = std::next(std::next(it));
                    if (it_next < std::prev(segments.end()))
                        // next_first_value = current_value_vector[it_next -> first];
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
                    auto it_next = std::next(it);
                    if (it_next < std::prev(segments.end()))
                        next_first_value = it_next -> intercept + corrections_vector[it_next -> first];
                    else
                        next_first_value = INT_MAX;
                }

                for (Covered_Value j = current_pos - first; j < covered; ++j) {
                    K value = ((j * slope_significand) >> slope_exponent) + intercept + corrections_vector[j + first];
                    if (value >= posting_value) {
                        current_pos = j + first;
                        return value;
                    }
                }
            }
            current_pos = n;
            return INT_MAX;
        }

        void next_normal() {
            if (current_value == INT_MAX) {
                return;
            }
            else if (current_pos >= n) {
                current_value = INT_MAX;
                return;
            }

            auto it = segments.begin() + current_segment;
            while (std::next(it) -> first <= current_pos && current_segment < segments.size() - 1) {
                it = std::next(it);
                current_segment++;
            }
            current_value = (((current_pos - it -> first) * it -> slope_significand) >> it -> slope_exponent) + it -> intercept + corrections_vector[current_pos];
            current_pos++;
        }

        void normal_decode_query() {
            uint32_t pointer = 0;
            for (auto it = segments.begin(); it < std::prev(segments.end()); ++it) {
                auto intercept = it -> intercept;
                auto exponent = it -> slope_exponent;
                auto significand = it -> slope_significand;
                auto covered = it -> covered;
                for (Covered_Value j = 0; j < covered; ++j) {
                    current_value_vector[pointer] = ((j * significand) >> exponent) + intercept + corrections_vector[pointer];
                    pointer++;
                }
            }
        }

        std::vector<K> normal_decode(bool warm_up=true) {
            std::vector<K> output;
            output.resize(n);
            if (warm_up) {
                for (int warm_time = 0; warm_time < 5; warm_time++) {
                    for (int i = 0; i < n; i ++) {
                        output[i] = 0;
                    }
                }
            }
            uint32_t pointer = 0;
            auto start = std::chrono::high_resolution_clock::now();
            for (auto it = segments.begin(); it != std::prev(segments.end()); ++it) {
                auto first = it -> first;
                auto intercept = it -> intercept;
                auto exponent = it -> slope_exponent;
                auto significand = it -> slope_significand;
                auto covered = it -> covered;
                for (Covered_Value j = 0; j < covered; ++j) {
                    output[pointer] = ((j * significand) >> exponent) + intercept + corrections_vector[pointer];
                    pointer++;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration = duration.count();
            return output;
        }

        // used for SIMD
        constexpr static K key_nums = 8;
        constexpr static K align_val = 64; // 64 bytes for avx512
        std::vector<segment> segments_sort; // resorted segments
        std::vector<Simd_Value*> slope_significand_simd;
        std::vector<Simd_Value*> slope_exponent_simd;
        std::vector<Intercept_Value*> intercept_simd;
        std::vector<Correction_Value*> corrections_simd;
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
        void aligned_delete(T* ptr) {
            std::free(ptr);
        }

        void free_memory(string simd_type = "simd") {
            if (simd_type == "simd") {
                for (auto i = 0; i < slope_significand_simd.size(); i++) {
                    aligned_delete(slope_significand_simd[i]);
                    aligned_delete(slope_exponent_simd[i]);
                    aligned_delete(intercept_simd[i]);
                    aligned_delete(corrections_simd[i]);
                    aligned_delete(first_simd[i]);
                    aligned_delete(covered_simd[i]);
                }
                vector<Covered_Value> ().swap(cover_length);
                vector<Correction_Value*> ().swap(corrections_simd);
                vector<Correction_Value> ().swap(corrections_vector_residual);
                vector<Simd_Value*> ().swap(slope_significand_simd);
                vector<Simd_Value*> ().swap(slope_exponent_simd);
                vector<Intercept_Value*> ().swap(intercept_simd);
                vector<Covered_Value*> ().swap(first_simd);
                vector<Covered_Value*> ().swap(covered_simd);
                vector<segment> ().swap(segments_sort);
                vector<K> ().swap(current_value_vector);
            }
        }

        // simple SIMD
        void create_corrections_simple() { // correction used for simple simd
            Correction_Value key_nums_tmp = key_nums * 2;
            for (auto it = segments.begin(); it < std::prev(segments.end()); it++) {
                auto &s = *it;
                auto covered = s.covered;
                auto first = s.first;
                Covered_Value max_min_covered = 0;
                for (auto i = 0; i + key_nums_tmp < covered; i = i + key_nums_tmp) {
                    max_min_covered++;
                    alignas(align_val) Correction_Value *corrections_simd_tmp = aligned_new<Correction_Value>(key_nums_tmp);
                    for (Covered_Value j = 0; j < key_nums_tmp; j++) {
                        corrections_simd_tmp[j] = corrections_vector[j + i + first];
                    }
                    corrections_simd.emplace_back(corrections_simd_tmp);
                }
                cover_length.emplace_back(max_min_covered);
            }
        }

        void create_corrections_residual_simple() {
            // Covered_Value correct_residual_pointers = -1;
            Covered_Value covered_length_pointer = -1;
            Covered_Value key_nums_tmp = key_nums * 2;
            for (auto it = segments.begin(); it < std::prev(segments.end()); it++) {
                auto &s = *it;
                auto covered = s.covered;
                auto first = s.first;
                Covered_Value cover_length_tmp = cover_length[++covered_length_pointer];
                for (Covered_Value i = cover_length_tmp * key_nums_tmp; i < covered; i++) {
                    Correction_Value tmp = corrections_vector[i + first];
                    corrections_vector_residual.emplace_back(tmp);
                }
            }
        }

        std::vector<K> simd_decode_512i_simple(bool warm_up=true){
            std::vector<K> output(n);
            total_calculated = 0;
            total_calculated_add = 0;
            total_duration = 0;
            if (warm_up) {
                for (int warm_time = 0; warm_time < 5; warm_time++) {
                    for (int i = 0; i < n; i ++) {
                        output[i] = 0;
                    }
                }
            }
            Covered_Value pointers = -1;
            Covered_Value correct_pointers = -1;
            Covered_Value correct_residual_pointers = -1;
            Covered_Value covered_length_pointer = -1;
            Covered_Value key_nums_tmp = key_nums * 2;
            conversion_time = segments.size();
            alignas(align_val) Correction_Value *result_int32 = aligned_new<Correction_Value>(key_nums_tmp);
            auto start1 = std::chrono::high_resolution_clock::now();

            for (auto it = segments.begin(); it < std::prev(segments.end()); it++) {
                auto &s = *it;
                Covered_Value covered = s.covered;
                Simd_Value significand = s.slope_significand;
                Simd_Value exponent = s.slope_exponent;
                Intercept_Value intercept = s.intercept;

                Covered_Value covered_length_tmp = cover_length[++covered_length_pointer];

                if (covered_length_tmp > 0) {
                    __m512i slope_correct_v, result_v, contact_int32_v, corrections_v;
                    __m256i int32_v1, int32_v2;

                    __m512i slope_significand_v = _mm512_set_epi64(7 * significand, 6 * significand, 5 * significand, 4 * significand, 3 * significand, 2 * significand, 1 * significand, 0);
                    __m512i slope_significand_v_add = _mm512_set1_epi64(8 * significand);
                    __m512i slope_exponent_v = _mm512_set1_epi64(exponent);
                    __m512i intercept_v = _mm512_set1_epi32(intercept);

                    for (Covered_Value i = 0; i < covered_length_tmp; i++) {
                        Correction_Value *corrections_p = corrections_simd[++correct_pointers];
                        slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                        int32_v1 = _mm512_cvtepi64_epi32(slope_correct_v); // lower 8
                        slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_add);

                        slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                        int32_v2 = _mm512_cvtepi64_epi32(slope_correct_v); // upper 8
                        slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_add);

                        contact_int32_v = _mm512_inserti64x4(contact_int32_v, int32_v1, 0);
                        contact_int32_v = _mm512_inserti64x4(contact_int32_v, int32_v2, 1);

                        corrections_v = _mm512_load_epi32(corrections_p);
                        result_v = _mm512_add_epi32(intercept_v, corrections_v);
                        result_v = _mm512_add_epi32(result_v, contact_int32_v);

                        _mm512_store_epi32(result_int32, result_v);
                        for (Covered_Value k = 0; k < key_nums_tmp; k++)
                            output[++pointers] = result_int32[k];
                        total_calculated++;
                    }

                }

                for (Covered_Value j = covered_length_tmp * key_nums_tmp; j < covered; ++j){
                    output[++pointers] = ((j * significand) >> exponent) + intercept + corrections_vector_residual[++correct_residual_pointers];
                }
            }

            auto end1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
            total_duration = duration1.count();
            total_calculated_add = total_calculated * 2;
            total_calculated = total_calculated * key_nums_tmp;
            return output;
        }

        // our SIMD
        void simd_init(bool use_max = false) {
            segments_sort = segments;
            std::sort(segments_sort.begin(), segments_sort.end(), [](const segment &a, const segment &b) {return a.covered > b.covered;});
            Covered_Value max_cover = 0;
            Covered_Value min_cover = INT_MAX;
            Covered_Value max_min_covered = 0;
            idx = 0;

            for (typename std::vector<segment>::iterator it = segments_sort.begin(); it + key_nums < std::prev(segments_sort.end()); it = it + key_nums) {
                alignas(align_val) Simd_Value *slope_significand_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Simd_Value *slope_exponent_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Intercept_Value *intercept_simd_tmp = aligned_new<Intercept_Value>(key_nums * 2);
                alignas(align_val) Covered_Value *covered_tmp = aligned_new<Covered_Value>(key_nums);
                alignas(align_val) Covered_Value *first_tmp = aligned_new<Covered_Value>(key_nums);

                std::vector<segment> simd_segments(it, it + key_nums);
                std::sort(simd_segments.begin(), simd_segments.end(), [](const segment &a, const segment &b) {return a.first < b.first;}); // part sorted
                if (simd_segments.back().covered < 32) {
                    break;
                }
                idx += key_nums;

                for (auto i = 0, its = simd_segments.begin(); its != simd_segments.end(); ++its, ++i) {
                    auto &s = *its;
                    auto covered = s.covered;
                    slope_significand_simd_tmp[i] = static_cast<Simd_Value>(s.slope_significand);
                    slope_exponent_simd_tmp[i] = static_cast<Simd_Value>(s.slope_exponent);
                    intercept_simd_tmp[i] = static_cast<Intercept_Value>(s.intercept);
                    intercept_simd_tmp[i + key_nums] = static_cast<Intercept_Value>(s.intercept);
                    covered_tmp[i] = static_cast<Covered_Value> (covered);
                    first_tmp[i] = static_cast<Covered_Value> (s.first);
                    min_cover = min_cover < covered ? min_cover : covered;
                    max_cover = max_cover > covered ? max_cover : covered;
                }

                slope_significand_simd.emplace_back(slope_significand_simd_tmp);
                slope_exponent_simd.emplace_back(slope_exponent_simd_tmp);
                intercept_simd.emplace_back(intercept_simd_tmp);
                first_simd.emplace_back(first_tmp);
                covered_simd.emplace_back(covered_tmp);
                max_min_covered = use_max ? max_cover : min_cover;
                if (max_min_covered % 2 == 1) // we need to make it even
                    max_min_covered--;
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
                if (cover_length_tmp > 0)
                    for (Covered_Value k = 0; k < key_nums; k++) {
                        corrections_simd_tmp[k] = corrections_vector[first_tmp[k]];
                    }
                for (Covered_Value j = 1; j < cover_length_tmp; j++) {
                    for (Covered_Value k = 0; k < key_nums; k++) {
                        corrections_simd_tmp[j * key_nums + k] = corrections_vector[first_tmp[k] + j];
                    }
                }
                corrections_simd[i] = corrections_simd_tmp;
            }
        }

        Covered_Value count_residual_size() {
            Covered_Value corrections_vector_residual_size = 0;
            for(int i = 0; i < slope_significand_simd.size(); i++) {
                Covered_Value *cover_tmp = covered_simd[i];
                Covered_Value cover_length_tmp = cover_length[i];
                for (Covered_Value k = 0; k < key_nums; k++) {
                    auto covered = cover_tmp[k];
                    if (cover_length_tmp < covered)
                        for (Covered_Value pos = cover_length_tmp; pos < covered; pos++)
                            corrections_vector_residual_size++;
                }
            }
            for(auto it = segments_sort.begin() + idx; it < std::prev(segments_sort.end()); it++) {
                auto& s = *it;
                auto covered = s.covered;
                for (Covered_Value pos = 0; pos < covered; pos++)
                    corrections_vector_residual_size++;
            }
            return corrections_vector_residual_size;
        }

        void create_corrections_residual() {
             corrections_vector_residual.resize( count_residual_size());
             Covered_Value correct_pointers = -1;
             Covered_Value key_nums_tmp = key_nums;
             for(int i = 0; i < slope_significand_simd.size(); i++) {
                 Covered_Value *first_tmp = first_simd[i];
                 Covered_Value *cover_tmp = covered_simd[i];
                 Covered_Value cover_length_tmp = cover_length[i];
                 for (K k = 0; k < key_nums_tmp; k++) {
                     if (cover_length_tmp < cover_tmp[k]) {
                         for (Simd_Value pos = cover_length_tmp; pos < cover_tmp[k]; pos++){
                             corrections_vector_residual[++correct_pointers] = corrections_vector[first_tmp[k] + pos];
                         }
                     }
                 }
             }
             for(auto it = segments_sort.begin() + idx; it < std::prev(segments_sort.end()); it++) {
                 for (Simd_Value pos = 0; pos < it -> covered; pos++){
                     corrections_vector_residual[++correct_pointers] = corrections_vector[it -> first + pos];
                 }
             }
        }

        std::vector<K> simd_decode_512i(bool warm_up=true){
            Correction_Value correct_pointers = -1;
            // alignas(align_val) K *output = aligned_new<K>(n);
            std::vector<K> output(n);
            total_calculated = 0;
            total_calculated_add = 0;
            total_duration = 0;
            if (warm_up) {
                for (int warm_time = 0; warm_time < 5; warm_time++) {
                    for (int i = 0; i < n; i ++) {
                        output[i] = 0;
                    }
                }
            }


            std::vector<Covered_Value*> first_pointer(8);
            // K *output_pointer = output.data();

            alignas(align_val) Correction_Value *result_int32 = aligned_new<Correction_Value>(16);
            // __m512i rerange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < slope_significand_simd.size(); i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them

                Covered_Value cover_length_tmp = cover_length[i];
                Covered_Value *cover_tmp = covered_simd[i];
                // Covered_Value  *first_tmp = first_simd[i];

                for (Covered_Value k = 0; k < 8; k++)
                    first_pointer[k] = output.data() + first_simd[i][k];
                //     // first_pointer[k + 8] = current_value_vector.data() + first_tmp[k] + 1;
                // }

                __m512i slope_significand_v_tmp = _mm512_load_epi64(slope_significand_simd[i]);
                __m512i slope_significand_v = slope_significand_v_tmp;
                __m512i slope_exponent_v = _mm512_load_epi64(slope_exponent_simd[i]);
                __m512i slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);

                __m256i int32_v = _mm256_setzero_si256();
                __m512i result_v = _mm512_castsi256_si512(int32_v);
                int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                result_v = _mm512_inserti64x4(result_v, int32_v, 1);

                __m512i intercept_v = _mm512_load_epi32(intercept_simd[i]);
                Correction_Value *corrections_p = corrections_simd[i];
                __m512i corrections_v = _mm512_load_epi32(corrections_p);

                result_v = _mm512_add_epi32(result_v, intercept_v);
                result_v = _mm512_add_epi32(result_v, corrections_v);
                // result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                _mm512_store_epi32(result_int32, result_v); // save the data

                for (uint8_t k = 0; k < 8; k++) {
                    *first_pointer[k]++ = result_int32[k];
                    // *first_pointer[k]++ = result_int32[k + 8];
                }
                for (uint8_t k = 8, first_idx = 0; first_idx < 8; first_idx++) {
                    *first_pointer[first_idx]++ = result_int32[k++];
                }

                // for (uint8_t k = 0, first_idx = 0; first_idx < 8; first_idx++) {
                //     *first_pointer[first_idx]++ = result_int32[k++];
                //     *first_pointer[first_idx]++ = result_int32[k++];
                // }

                // _mm512_storeu_epi32(output_pointer, result_v); // save the data
                // output_pointer += 16;


                total_calculated++;
                for (Covered_Value j = 2; j < cover_length_tmp; j= j + 2) {
                    slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    result_v = _mm512_castsi256_si512(int32_v);

                    slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    result_v = _mm512_inserti64x4(result_v, int32_v, 1);

                    result_v = _mm512_add_epi32(result_v , intercept_v);

                    corrections_p += 16;
                    corrections_v = _mm512_load_epi32(corrections_p);
                    result_v = _mm512_add_epi32(result_v, corrections_v);
                    // result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                    // _mm512_storeu_epi32(output_pointer, result_v);
                    // output_pointer += 16;
                    _mm512_store_epi32(result_int32, result_v);

                    // for (uint8_t k = 0, first_idx = 0; first_idx < 8; first_idx++) {
                    //     *first_pointer[first_idx]++ = result_int32[k++];
                    //     *first_pointer[first_idx]++ = result_int32[k++];
                    // }
                    // _mm512_i32scatter_epi32

                    for (uint8_t k = 0; k < 8; k++) {
                        *first_pointer[k]++ = result_int32[k];
                        // *first_pointer[k]++ = result_int32[k + 8];
                    }
                    for (uint8_t k = 8, first_idx = 0; first_idx < 8; first_idx++) {
                        *first_pointer[first_idx]++ = result_int32[k++];
                    }
                    total_calculated++;
                }

                for (uint8_t k = 0; k < 8; k++) {
                    auto covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        auto slope_significand = slope_significand_simd[i][k];
                        auto slope_exponent = slope_exponent_simd[i][k];
                        auto intercept = intercept_simd[i][k];
                        for (Covered_Value pos = cover_length_tmp; pos < covered; pos++)

                            // current_value_vector[decode_result_map.find_value(++pointers)] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];

                        // current_value_vector[first_tmp[k] + pos] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                            *first_pointer[k]++ = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                            // *output_pointer++ = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                    }
                }
            }

            for(auto it = segments_sort.begin() + idx; it < std::prev(segments_sort.end()); it++) {
                auto covered = it -> covered;
                auto slope_significand = it -> slope_significand;
                auto slope_exponent = it -> slope_exponent;
                auto intercept = it -> intercept;
                auto first = it -> first;
                for (Covered_Value pos = 0; pos < covered; pos++)
                    output[first + pos] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                    // *output_pointer++ = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];

                // current_value_vector[it -> first + pos] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
            }
            // std::sort(output.begin(), output.end());
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration = duration.count();
            total_calculated = total_calculated * 16;
            // aligned_delete(output);
            // output.clear();
            return output;
        }

        // void simd_decode_512i_query(){
        //     Correction_Value correct_pointers = -1;
        //     alignas(64) Correction_Value *result_int32 = aligned_new<Correction_Value>(16);
        //     __m512i rerange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
        //     std::array<Covered_Value*, 8> first_pointer;
        //     for (int i = 0; i < slope_significand_simd.size(); i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them
        //
        //         const Covered_Value cover_length_tmp = cover_length[i];
        //         Covered_Value *const cover_tmp = covered_simd[i];
        //         Covered_Value  *const first_tmp = first_simd[i];
        //
        //         first_pointer[0] = current_value_vector.data() + first_tmp[0];
        //         first_pointer[1] = current_value_vector.data() + first_tmp[1];
        //         first_pointer[2] = current_value_vector.data() + first_tmp[2];
        //         first_pointer[3] = current_value_vector.data() + first_tmp[3];
        //         first_pointer[4] = current_value_vector.data() + first_tmp[4];
        //         first_pointer[5] = current_value_vector.data() + first_tmp[5];
        //         first_pointer[6] = current_value_vector.data() + first_tmp[6];
        //         first_pointer[7] = current_value_vector.data() + first_tmp[7];
        //
        //         const __m512i slope_significand_v_tmp = _mm512_load_epi64(slope_significand_simd[i]);
        //         __m512i slope_significand_v = slope_significand_v_tmp;
        //         const __m512i slope_exponent_v = _mm512_load_epi64(slope_exponent_simd[i]);
        //         __m512i slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
        //
        //         __m256i int32_v = _mm256_setzero_si256();
        //         __m512i result_v = _mm512_castsi256_si512(int32_v);
        //         int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
        //         result_v = _mm512_inserti64x4(result_v, int32_v, 1);
        //
        //         const __m512i intercept_v = _mm512_load_epi32(intercept_simd[i]);
        //         Correction_Value *corrections_p = corrections_simd[i];
        //         __m512i corrections_v = _mm512_load_epi32(corrections_p);
        //
        //         result_v = _mm512_add_epi32(result_v, intercept_v);
        //         result_v = _mm512_add_epi32(result_v, corrections_v);
        //
        //         result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);
        //
        //         _mm512_store_epi32(result_int32, result_v); // save the data
        //         // for (uint8_t k = 0; k < 8; k++) {
        //         //     *first_pointer[k]++ = result_int32[k];
        //         // }
        //         // for (uint8_t k = 8, first_idx = 0; first_idx < 8; first_idx++) {
        //         //     *first_pointer[first_idx]++ = result_int32[k++];
        //         // }
        //         *first_pointer[0]++ = result_int32[0];
        //         *first_pointer[0]++ = result_int32[1];
        //         *first_pointer[1]++ = result_int32[2];
        //         *first_pointer[1]++ = result_int32[3];
        //         *first_pointer[2]++ = result_int32[4];
        //         *first_pointer[2]++ = result_int32[5];
        //         *first_pointer[3]++ = result_int32[6];
        //         *first_pointer[3]++ = result_int32[7];
        //         *first_pointer[4]++ = result_int32[8];
        //         *first_pointer[4]++ = result_int32[9];
        //         *first_pointer[5]++ = result_int32[10];
        //         *first_pointer[5]++ = result_int32[11];
        //         *first_pointer[6]++ = result_int32[12];
        //         *first_pointer[6]++ = result_int32[13];
        //         *first_pointer[7]++ = result_int32[14];
        //         *first_pointer[7]++ = result_int32[15];
        //
        //         // *first_pointer[0]++ = result_int32[0];
        //         // *first_pointer[0]++ = result_int32[8];
        //         // *first_pointer[1]++ = result_int32[1];
        //         // *first_pointer[1]++ = result_int32[9];
        //         // *first_pointer[2]++ = result_int32[2];
        //         // *first_pointer[2]++ = result_int32[10];
        //         // *first_pointer[3]++ = result_int32[3];
        //         // *first_pointer[3]++ = result_int32[11];
        //         // *first_pointer[4]++ = result_int32[4];
        //         // *first_pointer[4]++ = result_int32[12];
        //         // *first_pointer[5]++ = result_int32[5];
        //         // *first_pointer[5]++ = result_int32[13];
        //         // *first_pointer[6]++ = result_int32[6];
        //         // *first_pointer[6]++ = result_int32[14];
        //         // *first_pointer[7]++ = result_int32[7];
        //         // *first_pointer[7]++ = result_int32[15];
        //
        //
        //         // for (uint8_t k = 0, first_idx = 0; first_idx < 8; first_idx++) {
        //         //     *first_pointer[first_idx]++ = result_int32[k++];
        //         //     *first_pointer[first_idx]++ = result_int32[k++];
        //         // }
        //
        //         for (Covered_Value j = 2; j < cover_length_tmp; j= j + 2) {
        //             slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
        //             slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
        //             int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
        //             result_v = _mm512_castsi256_si512(int32_v);
        //
        //             slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
        //             slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
        //             int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
        //             result_v = _mm512_inserti64x4(result_v, int32_v, 1);
        //
        //             result_v = _mm512_add_epi32(result_v , intercept_v);
        //
        //             corrections_p += 16;
        //             corrections_v = _mm512_load_epi32(corrections_p);
        //             result_v = _mm512_add_epi32(result_v, corrections_v);
        //
        //             result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);
        //
        //             _mm512_store_epi32(result_int32, result_v);
        //             *first_pointer[0]++ = result_int32[0];
        //             *first_pointer[0]++ = result_int32[1];
        //             *first_pointer[1]++ = result_int32[2];
        //             *first_pointer[1]++ = result_int32[3];
        //             *first_pointer[2]++ = result_int32[4];
        //             *first_pointer[2]++ = result_int32[5];
        //             *first_pointer[3]++ = result_int32[6];
        //             *first_pointer[3]++ = result_int32[7];
        //             *first_pointer[4]++ = result_int32[8];
        //             *first_pointer[4]++ = result_int32[9];
        //             *first_pointer[5]++ = result_int32[10];
        //             *first_pointer[5]++ = result_int32[11];
        //             *first_pointer[6]++ = result_int32[12];
        //             *first_pointer[6]++ = result_int32[13];
        //             *first_pointer[7]++ = result_int32[14];
        //             *first_pointer[7]++ = result_int32[15];
        //
        //             // *first_pointer[0]++ = result_int32[0];
        //             // *first_pointer[0]++ = result_int32[8];
        //             // *first_pointer[1]++ = result_int32[1];
        //             // *first_pointer[1]++ = result_int32[9];
        //             // *first_pointer[2]++ = result_int32[2];
        //             // *first_pointer[2]++ = result_int32[10];
        //             // *first_pointer[3]++ = result_int32[3];
        //             // *first_pointer[3]++ = result_int32[11];
        //             // *first_pointer[4]++ = result_int32[4];
        //             // *first_pointer[4]++ = result_int32[12];
        //             // *first_pointer[5]++ = result_int32[5];
        //             // *first_pointer[5]++ = result_int32[13];
        //             // *first_pointer[6]++ = result_int32[6];
        //             // *first_pointer[6]++ = result_int32[14];
        //             // *first_pointer[7]++ = result_int32[7];
        //             // *first_pointer[7]++ = result_int32[15];
        //
        //         }
        //
        //         for (Covered_Value k = 0; k < 8; k++) {
        //             const Covered_Value covered = cover_tmp[k];
        //             if (cover_length_tmp < covered) {
        //                 const Simd_Value slope_significand = slope_significand_simd[i][k];
        //                 const Simd_Value slope_exponent = slope_exponent_simd[i][k];
        //                 const Intercept_Value intercept = intercept_simd[i][k];
        //                 for (Covered_Value pos = cover_length_tmp; pos < covered; pos++)
        //                     *first_pointer[k]++ = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
        //                         // current_value_vector[first_tmp[k] + pos] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
        //             }
        //         }
        //     }
        //
        //     const auto end_iter = std::prev(segments_sort.end());
        //     for (auto it = segments_sort.begin() + idx; it < end_iter; ++it) {
        //         const auto& seg = *it;
        //         K* dst = current_value_vector.data() + seg.first;
        //         for (Covered_Value pos = 0; pos < seg.covered; ++pos) {
        //             dst[pos] = ((seg.slope_significand * pos) >> seg.slope_exponent) + seg.intercept + corrections_vector_residual[++correct_pointers];
        //         }
        //     }
        // }

        void simd_decode_512i_query(){
            Correction_Value correct_pointers = -1;
            alignas(64) Correction_Value *result_int32 = aligned_new<Correction_Value>(16);
            std::array<K*, 8> first_pointer;
            __m512i rerange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            for (int i = 0; i < slope_significand_simd.size(); i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them

                const Covered_Value cover_length_tmp = cover_length[i];
                Covered_Value *const cover_tmp = covered_simd[i];
                Covered_Value  *const first_tmp = first_simd[i];

                first_pointer[0] = current_value_vector.data() + first_tmp[0];
                first_pointer[1] = current_value_vector.data() + first_tmp[1];
                first_pointer[2] = current_value_vector.data() + first_tmp[2];
                first_pointer[3] = current_value_vector.data() + first_tmp[3];
                first_pointer[4] = current_value_vector.data() + first_tmp[4];
                first_pointer[5] = current_value_vector.data() + first_tmp[5];
                first_pointer[6] = current_value_vector.data() + first_tmp[6];
                first_pointer[7] = current_value_vector.data() + first_tmp[7];

                const __m512i slope_significand_v_tmp = _mm512_load_epi64(slope_significand_simd[i]);
                __m512i slope_significand_v = slope_significand_v_tmp;
                const __m512i slope_exponent_v = _mm512_load_epi64(slope_exponent_simd[i]);
                __m512i slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);

                __m256i int32_v = _mm256_setzero_si256();
                __m512i result_v = _mm512_castsi256_si512(int32_v);
                int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                result_v = _mm512_inserti64x4(result_v, int32_v, 1);

                const __m512i intercept_v = _mm512_load_epi32(intercept_simd[i]);
                Correction_Value *corrections_p = corrections_simd[i];
                __m512i corrections_v = _mm512_load_epi32(corrections_p);

                result_v = _mm512_add_epi32(result_v, intercept_v);
                result_v = _mm512_add_epi32(result_v, corrections_v);

                result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                __m128i t0 = _mm512_extracti32x4_epi32(result_v, 0);
                __m128i t1 = _mm512_extracti32x4_epi32(result_v, 1);
                __m128i t2 = _mm512_extracti32x4_epi32(result_v, 2);
                __m128i t3 = _mm512_extracti32x4_epi32(result_v, 3);

                // 对每个 128 位块进行 shuffle 和存储
                __m128i t0_shuffled = _mm_shuffle_epi32(t0, _MM_SHUFFLE(3, 2, 3, 2));
                __m128i t1_shuffled = _mm_shuffle_epi32(t1, _MM_SHUFFLE(3, 2, 3, 2));
                __m128i t2_shuffled = _mm_shuffle_epi32(t2, _MM_SHUFFLE(3, 2, 3, 2));
                __m128i t3_shuffled = _mm_shuffle_epi32(t3, _MM_SHUFFLE(3, 2, 3, 2));

                // 合并存储操作
                _mm_storel_epi64((__m128i*)first_pointer[0], t0);
                _mm_storel_epi64((__m128i*)first_pointer[1], t0_shuffled);
                _mm_storel_epi64((__m128i*)first_pointer[2], t1);
                _mm_storel_epi64((__m128i*)first_pointer[3], t1_shuffled);
                _mm_storel_epi64((__m128i*)first_pointer[4], t2);
                _mm_storel_epi64((__m128i*)first_pointer[5], t2_shuffled);
                _mm_storel_epi64((__m128i*)first_pointer[6], t3);
                _mm_storel_epi64((__m128i*)first_pointer[7], t3_shuffled);


                // _mm_storel_epi64((__m128i*)first_pointer[0], t0);
                // _mm_storel_epi64((__m128i*)first_pointer[1], _mm_shuffle_epi32(t0, _MM_SHUFFLE(3, 2, 3, 2)));
                // _mm_storel_epi64((__m128i*)first_pointer[2], t1);
                // _mm_storel_epi64((__m128i*)first_pointer[3], _mm_shuffle_epi32(t1, _MM_SHUFFLE(3, 2, 3, 2)));
                // _mm_storel_epi64((__m128i*)first_pointer[4], t2);
                // _mm_storel_epi64((__m128i*)first_pointer[5], _mm_shuffle_epi32(t2, _MM_SHUFFLE(3, 2, 3, 2)));
                // _mm_storel_epi64((__m128i*)first_pointer[6], t3);
                // _mm_storel_epi64((__m128i*)first_pointer[7], _mm_shuffle_epi32(t3, _MM_SHUFFLE(3, 2, 3, 2)));


                first_pointer[0] += 2;
                first_pointer[1] += 2;
                first_pointer[2] += 2;
                first_pointer[3] += 2;
                first_pointer[4] += 2;
                first_pointer[5] += 2;
                first_pointer[6] += 2;
                first_pointer[7] += 2;

                for (Covered_Value j = 2; j < cover_length_tmp; j= j + 2) {
                    slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    result_v = _mm512_castsi256_si512(int32_v);

                    slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    result_v = _mm512_inserti64x4(result_v, int32_v, 1);

                    result_v = _mm512_add_epi32(result_v , intercept_v);

                    corrections_p += 16;
                    corrections_v = _mm512_load_epi32(corrections_p);
                    result_v = _mm512_add_epi32(result_v, corrections_v);

                    result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                    t0 = _mm512_extracti32x4_epi32(result_v, 0);
                    t1 = _mm512_extracti32x4_epi32(result_v, 1);
                    t2 = _mm512_extracti32x4_epi32(result_v, 2);
                    t3 = _mm512_extracti32x4_epi32(result_v, 3);

                    t0_shuffled = _mm_shuffle_epi32(t0, _MM_SHUFFLE(3, 2, 3, 2));
                    t1_shuffled = _mm_shuffle_epi32(t1, _MM_SHUFFLE(3, 2, 3, 2));
                    t2_shuffled = _mm_shuffle_epi32(t2, _MM_SHUFFLE(3, 2, 3, 2));
                    t3_shuffled = _mm_shuffle_epi32(t3, _MM_SHUFFLE(3, 2, 3, 2));

                    // 合并存储操作
                    _mm_storel_epi64((__m128i*)first_pointer[0], t0);
                    _mm_storel_epi64((__m128i*)first_pointer[1], t0_shuffled);
                    _mm_storel_epi64((__m128i*)first_pointer[2], t1);
                    _mm_storel_epi64((__m128i*)first_pointer[3], t1_shuffled);
                    _mm_storel_epi64((__m128i*)first_pointer[4], t2);
                    _mm_storel_epi64((__m128i*)first_pointer[5], t2_shuffled);
                    _mm_storel_epi64((__m128i*)first_pointer[6], t3);
                    _mm_storel_epi64((__m128i*)first_pointer[7], t3_shuffled);

                    // _mm_storel_epi64((__m128i*)first_pointer[0], t0);
                    // _mm_storel_epi64((__m128i*)first_pointer[1], _mm_shuffle_epi32(t0, _MM_SHUFFLE(3, 2, 3, 2)));
                    // _mm_storel_epi64((__m128i*)first_pointer[2], t1);
                    // _mm_storel_epi64((__m128i*)first_pointer[3], _mm_shuffle_epi32(t1, _MM_SHUFFLE(3, 2, 3, 2)));
                    // _mm_storel_epi64((__m128i*)first_pointer[4], t2);
                    // _mm_storel_epi64((__m128i*)first_pointer[5], _mm_shuffle_epi32(t2, _MM_SHUFFLE(3, 2, 3, 2)));
                    // _mm_storel_epi64((__m128i*)first_pointer[6], t3);
                    // _mm_storel_epi64((__m128i*)first_pointer[7], _mm_shuffle_epi32(t3, _MM_SHUFFLE(3, 2, 3, 2)));
                    first_pointer[0] += 2;
                    first_pointer[1] += 2;
                    first_pointer[2] += 2;
                    first_pointer[3] += 2;
                    first_pointer[4] += 2;
                    first_pointer[5] += 2;
                    first_pointer[6] += 2;
                    first_pointer[7] += 2;
                }

                for (Covered_Value k = 0; k < 8; k++) {
                    const Covered_Value covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        const Simd_Value slope_significand = slope_significand_simd[i][k];
                        const Simd_Value slope_exponent = slope_exponent_simd[i][k];
                        const Intercept_Value intercept = intercept_simd[i][k];
                        // Covered_Value* dst = first_pointer[k];
                        for (Covered_Value pos = cover_length_tmp; pos < covered; pos++)
                            *first_pointer[k]++ = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                                // *dst++ = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                    }
                }
            }

            const auto end_iter = std::prev(segments_sort.end());
            for (auto it = segments_sort.begin() + idx; it < end_iter; ++it) {
                const auto& seg = *it;
                K* dst = current_value_vector.data() + seg.first;
                for (Covered_Value pos = 0; pos < seg.covered; ++pos) {
                    dst[pos] = ((seg.slope_significand * pos) >> seg.slope_exponent) + seg.intercept + corrections_vector_residual[++correct_pointers];
                }
            }
        }

    };
}