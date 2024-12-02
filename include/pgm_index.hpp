// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once
#include <immintrin.h>

#include <climits>
#include <iostream>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <cmath>
#include "sdsl/bits.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "tools.hpp"
#include "piecewise_linear_model.hpp"

namespace pgm
{
    #define BIT_CEIL(x) ((x) < 2 ? 1u : 1u << (64u - __builtin_clzll((x) - 1)))
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

    template <typename K, uint64_t Epsilon = 64, uint64_t EpsilonRecursive = 0, typename Floating = float>
    // template <typename K, uint64_t Epsilon = 64, uint64_t EpsilonRecursive = 0, typename Floating = double>
    class PGMIndex {
    public:

        static_assert(Epsilon > 0);

        struct Segment;

        typedef int64_t Simd_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;

        uint64_t n;                           ///< The number of elements this index was built on.
        K first_pos;                     ///< The smallest element.
        std::vector<Segment> segments;      ///< The segments composing the index.
        std::vector<uint64_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

        sdsl::bit_vector signs; // signs for saving
        sdsl::rrr_vector<> signs_rrr; // signs for decode
        sdsl::int_vector<64> corrections; // corrections for saving
        std::vector<Correction_Value> corrections_vector; // corrections for decode

        /// Sentinel value to avoid bounds checking.
        static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity() : std::numeric_limits<K>::max();

        static constexpr bool auto_bpc = false;
        static constexpr uint64_t cache_line_bits = 64 * CHAR_BIT;
        static constexpr uint64_t extraction_density = auto_bpc ? 1 : BIT_CEIL(4 * cache_line_bits / BIT_WIDTH(Epsilon));

        using position_type = typename std::conditional_t<sizeof(K) <= 4, uint32_t, uint64_t>;
        using larger_signed_key_type = typename std::conditional_t<sizeof(K) <= 4, uint64_t, __int128>;
        using canonical_segment = typename internal::OptimalPiecewiseLinearModel<position_type, K>::CanonicalSegment;

        uint64_t errorPointCount;

        PGMIndex(): n(0), errorPointCount(0) {}

        explicit PGMIndex(const std::vector<K>& data) : PGMIndex(data.begin(), data.end()) {}

        template <typename RandomIt>
        PGMIndex(RandomIt begin, RandomIt end):
            n(std::distance(begin, end)),
            first_pos(n ? *begin : K(0)),
            segments(),
            levels_offsets(),
            errorPointCount(0),
            corrections((n) * BIT_WIDTH(Epsilon) / 64 + 5), // set 5 to avoid insufficient memory
            signs(n),
            signs_rrr() {
            build(begin, end, Epsilon, EpsilonRecursive, segments, levels_offsets, errorPointCount, corrections,  signs, signs_rrr);
        }

        template <typename RandomIt>
        static void build(RandomIt begin, RandomIt end,
                          uint64_t epsilon,
                          uint64_t epsilon_recursive, // we don't use this parameter to build the recursive PGM-index
                          std::vector<Segment>& segments,
                          std::vector<uint64_t>& levels_offsets,
                          uint64_t& errorPointCount,
                          sdsl::int_vector<64>& corrections,
                          sdsl::bit_vector& signs,
                          sdsl::rrr_vector<>& signs_rrr
                          ) {

            auto n = (uint64_t) std::distance(begin, end);
            if (n == 0)
                return;

            levels_offsets.push_back(0);
            segments.reserve(n / (epsilon * epsilon));

            if (*std::prev(--end) == sentinel)
                throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

            uint64_t corrections_offset = 0;

            std::vector<canonical_segment> canonical_segments;

            canonical_segments.reserve(epsilon > 0 ? n / (epsilon * epsilon) : n / 8);

            auto in_fun = [begin](auto i) { return std::pair<position_type, K>(i, begin[i]); };
            auto out_fun = [&canonical_segments](auto cs) { canonical_segments.push_back(cs); };
            auto n_segments = internal::make_segmentation_par(n, epsilon, in_fun, out_fun);
            auto last_n = n_segments;

            for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it)
            {
                auto i = it->get_first_x();
                auto j = std::next(it) != canonical_segments.end() ? std::next(it)->get_first_x() : n;
                uint8_t bpc = BIT_WIDTH(epsilon);
                segments.emplace_back(*it, n, i, j, begin, errorPointCount, corrections.data(), signs, bpc, corrections_offset); // build the segment
                corrections_offset += bpc * (j - i);
            }
            // segments contains every first pos , interception and slope
            segments.emplace_back(n);
            levels_offsets.emplace_back(segments.size());
            // signs_rrr = sdsl::rrr_vector<>(signs); // rrr is not efficient
        }

        uint64_t segments_count() const { return segments.empty() ? 0 : levels_offsets[1] - 1; }

        uint64_t height() const { return levels_offsets.size() - 1; }

        std::vector<Segment> get_segments() const { return segments; }

        auto segment_for_pos(const K& pos) const {
            return std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), pos));
        }

        uint64_t indexSegments_bytes() const {
            return segments.size() * sizeof(Segment) + levels_offsets.size() * sizeof(uint64_t);
        }

        uint64_t segment_size_in_bytes() const {
//            return segments_count() * (sizeof(Segment));
            return segments_count() * (sizeof(Covered_Value) + sizeof(Intercept_Value) + sizeof(uint64_t) + sizeof(uint8_t)); // first + Intercept + slope_significant + slope_exponent
        }

        uint64_t corrections_size_in_bytes() const {
            return corrections.bit_size() / CHAR_BIT;
        }

        uint64_t signs_size_in_bytes() const {
            return signs.bit_size() / CHAR_BIT;
        }

        uint64_t signs_rrr_size_in_bytes() const {
            return size_in_bytes(signs_rrr);
        }

        // used for SIMD
        constexpr static K key_nums = 8;
        constexpr static K align_val = 64; // 64 bytes for avx512
        std::vector<Simd_Value> signs_Floating; //correction is already Floating
        std::vector<Segment> segments_sort; // resorted segments
        std::vector<Simd_Value*> slope_simd;
        std::vector<Simd_Value*> slope_significand_simd;
        std::vector<Simd_Value*> slope_exponent_simd;
        std::vector<Intercept_Value*> intercept_simd;
        std::vector<Correction_Value*> corrections_simd;
        std::vector<Correction_Value> corrections_vector_residual;
        std::vector<Covered_Value*> first_simd;
        std::vector<Covered_Value*> covered_simd;
        std::vector<Covered_Value> cover_length;

        Simd_Value correction_pos = 0;
        uint64_t duration = 0;
        uint64_t total_calculated = 0;
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

        void free_memory() {
            for (auto i = 0; i < slope_significand_simd.size(); i++) {
                aligned_delete(slope_significand_simd[i]);
                aligned_delete(slope_exponent_simd[i]);
                aligned_delete(intercept_simd[i]);
                aligned_delete(corrections_simd[i]);
                aligned_delete(first_simd[i]);
                aligned_delete(covered_simd[i]);
            }
            vector<Covered_Value> ().swap(cover_length);
            vector<Correction_Value> ().swap(corrections_vector_residual);
            vector<Simd_Value*> ().swap(slope_significand_simd);
            vector<Simd_Value*> ().swap(slope_exponent_simd);
            vector<Intercept_Value*> ().swap(intercept_simd);
            vector<Correction_Value*> ().swap(corrections_simd);
            vector<Covered_Value*> ().swap(first_simd);
            vector<Covered_Value*> ().swap(covered_simd);
            vector<Segment> ().swap(segments_sort);
        }

        std::vector<Covered_Value> result_map_init() {
            std::vector<Covered_Value> result_map(n);
            Simd_Value pointers = -1;
            Simd_Value key_nums_tmp = key_nums;
            for(int i = 0; i < slope_significand_simd.size(); i++) {
                Covered_Value *first_tmp = first_simd[i];
                Covered_Value *cover_tmp = covered_simd[i];
                Covered_Value cover_length_tmp = cover_length[i];

                for (Covered_Value k = 0; k < key_nums_tmp; k++)
                    result_map[first_tmp[k]] = ++pointers;

                for (Covered_Value j = 1; j < cover_length_tmp; j++)
                    for (K k = 0; k < key_nums_tmp; k++)
                        result_map[first_tmp[k] + j] = ++pointers;

                for (K k = 0; k < key_nums_tmp; k++) {
                    auto first = first_tmp[k];
                    auto covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        for (Covered_Value pos = cover_length_tmp; pos < covered; pos++)
                            result_map[first + pos] = ++pointers;
                    }
                }
            }
            for(auto it = segments_sort.begin() + idx; it < std::prev(segments_sort.end()); it++) {
                auto& s = *it;
                auto covered = s.covered;
                auto first = s.first;
                for (Covered_Value pos = 0; pos < covered; pos++)
                    result_map[first + pos] = ++pointers;
            }
            return result_map;
        }

        void normal_init(){
            corrections_vector.resize(n, 0);
            for(auto &s : segments) {
                auto covered = s.covered;
                auto first = s.first;
                for (Covered_Value j = 0; j < covered; ++j)
                    corrections_vector[j + first] = get_correction(corrections.data(), n, j + first, signs); // change to rrr_vector
                    // corrections_vector[j + first] = get_correction(corrections.data(), n, j + first, signs_rrr); // change to rrr_vector
            }
        }

        void simd_init(bool use_max = true) {
            segments_sort = segments;
            sort(segments_sort.begin(), segments_sort.end());
            Covered_Value max_cover = 0;
            Covered_Value min_cover = INT_MAX;
            Covered_Value max_min_covered = 0;
            idx = 0;

            for (typename std::vector<Segment>::iterator it = segments_sort.begin(); it + key_nums < std::prev(segments_sort.end()); it = it + key_nums) {
                alignas(align_val) Simd_Value *slope_significand_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Simd_Value *slope_exponent_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Intercept_Value *intercept_simd_tmp = aligned_new<Intercept_Value>(key_nums * 2);
                alignas(align_val) Covered_Value *covered_tmp = aligned_new<Covered_Value>(key_nums);
                alignas(align_val) Covered_Value *first_tmp = aligned_new<Covered_Value>(key_nums);

                std::vector<Segment> simd_segments(it, it + key_nums);
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

        void create_corrections() {
             corrections_simd.resize(cover_length.size(), 0);
            for (int i = 0;i < cover_length.size(); i++) {
                Covered_Value cover_length_tmp = cover_length[i];
                Covered_Value *first_tmp = first_simd[i];
                // Simd_Value *cover_tmp = covered_simd[i];
                alignas(align_val) Correction_Value *corrections_simd_tmp = aligned_new<Correction_Value>(cover_length_tmp * key_nums);
                if (cover_length_tmp > 0)
                    for (Covered_Value k = 0; k < key_nums; k++) {
                        corrections_simd_tmp[k] = get_correction(corrections.data(), n, first_tmp[k], signs);
                        // corrections_simd_tmp[k] = get_correction(corrections.data(), n, first_tmp[k], signs_rrr);
                    }
                for (Covered_Value j = 1; j < cover_length_tmp; j++) {
                    for (Covered_Value k = 0; k < key_nums; k++) {
                        corrections_simd_tmp[j * key_nums + k] = get_correction(corrections.data(), n, j + first_tmp[k], signs);
                        // corrections_simd_tmp[j * key_nums + k] = get_correction(corrections.data(), n, j + first_tmp[k], signs_rrr);
                    }
                }
                corrections_simd[i] = corrections_simd_tmp;
            }
            create_corrections_residual();
        }

        void create_corrections_residual() {
             corrections_vector_residual.resize(count_residual_size());
             Covered_Value correct_pointers = -1;
             Covered_Value key_nums_tmp = key_nums;
             for(int i = 0; i < slope_significand_simd.size(); i++) {
                 Covered_Value *first_tmp = first_simd[i];
                 Covered_Value *cover_tmp = covered_simd[i];
                 Covered_Value cover_length_tmp = cover_length[i];
                 for (K k = 0; k < key_nums_tmp; k++) {
                     if (cover_length_tmp < cover_tmp[k]) {
                         for (Simd_Value pos = cover_length_tmp; pos < cover_tmp[k]; pos++){
                             corrections_vector_residual[++correct_pointers] = (get_correction(corrections.data(), n, first_tmp[k] + pos, signs));
                             // corrections_vector_residual[++correct_pointers] = (get_correction(corrections.data(), n, first_tmp[k] + pos, signs_rrr));
                         }
                     }
                 }
             }
             for(auto it = segments_sort.begin() + idx; it < std::prev(segments_sort.end()); it++) {
                 auto& s = *it;
                 for (Simd_Value pos = 0; pos < s.covered; pos++){
                     corrections_vector_residual[++correct_pointers] = get_correction(corrections.data(), n, s.first + pos, signs);
                     // corrections_vector_residual[++correct_pointers] = (get_correction(corrections.data(), n, s.first + pos, signs_rrr));
                 }
             }
         }

        void create_corrections_simple() {
            Covered_Value key_nums_tmp = key_nums * 2;
            for (auto it = segments.begin(); it < std::prev(segments.end()); it++) {
                auto &s = *it;
                auto covered = s.covered;
                auto first = s.first;
                Covered_Value max_min_covered = 0;
                for (auto i = 0; i + key_nums_tmp < covered; i = i + key_nums_tmp) {
                    max_min_covered++;
                    alignas(align_val) Correction_Value *corrections_simd_tmp = aligned_new<Intercept_Value>(key_nums_tmp);
                    for (Covered_Value j = 0; j < key_nums_tmp; j++) {
                        corrections_simd_tmp[j] = get_correction(corrections.data(), n, j + i + first, signs);
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
                    Correction_Value tmp = get_correction(corrections.data(), n, i + first, signs);
                    corrections_vector_residual.emplace_back(tmp);
                }
            }
        }

        std::vector<K> normal_decode() {
            std::vector<K> out;
            out.resize(n);
            auto start1 = std::chrono::high_resolution_clock::now();
            for (auto it = segments.begin(); it != std::prev(segments.end()); ++it) {
                auto& s = *it;
                auto covered = s.covered;
                auto significand = s.slope_significand;
                auto exponent = s.slope_exponent;
                auto intercept = s.intercept;
                auto first = s.first;
                for (Covered_Value j = 0; j < covered; ++j){
                    out[j + first] = ((j * significand) >> exponent) + intercept + corrections_vector[j +first];
                }
            }
            auto end1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
            duration = duration1.count();
            return out;
        }

        std::vector<K> simd_decode_512i(){
            std::vector<K> output;
            output.resize(n);
            Covered_Value pointers = -1;
            Correction_Value correct_pointers = -1;
            Covered_Value key_nums_tmp = key_nums * 2;
            alignas(align_val) Correction_Value *result_int32 = aligned_new<Correction_Value>(key_nums_tmp);
            auto start1 = std::chrono::high_resolution_clock::now();

            for(int i = 0; i < slope_significand_simd.size(); i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them
                Covered_Value cover_length_tmp = cover_length[i];
                if (cover_length_tmp == 0)
                    continue;
                Covered_Value *cover_tmp = covered_simd[i];

                __m512i slope_correct_v, result_v, contact_int32_v;
                __m256i int32_v2, int32_v1;

                Intercept_Value *intercept_p = intercept_simd[i];
                __m512i intercept_v = _mm512_load_epi32(intercept_p);
                Correction_Value *corrections_p = corrections_simd[i];
                __m512i corrections_v = _mm512_load_epi32(corrections_p);
                result_v = _mm512_add_epi32(intercept_v, corrections_v);

                Simd_Value *slope_significand_p = slope_significand_simd[i];
                __m512i slope_significand_v_tmp = _mm512_load_epi64(slope_significand_p);
                __m512i slope_significand_v = _mm512_load_epi64(slope_significand_p);

                Simd_Value *slope_exponent_p = slope_exponent_simd[i];
                __m512i slope_exponent_v = _mm512_load_epi64(slope_exponent_p);

                slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                int32_v2 = _mm512_cvtepi64_epi32(slope_correct_v);
                int32_v1 = _mm256_setzero_si256();
                contact_int32_v = _mm512_inserti64x4(contact_int32_v, int32_v1, 0);
                contact_int32_v = _mm512_inserti64x4(contact_int32_v, int32_v2, 1);
                result_v = _mm512_add_epi32(result_v, contact_int32_v);
                _mm512_store_epi32(result_int32, result_v); // save the data
                for (Covered_Value k = 0; k < key_nums_tmp; k++)
                    output[++pointers] = result_int32[k];

                for (Covered_Value j = 2; j < cover_length_tmp; j= j + 2) {
                    slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                    int32_v1 = _mm512_cvtepi64_epi32(slope_correct_v);
                    contact_int32_v = _mm512_inserti64x4(contact_int32_v, int32_v1, 0);

                    slope_significand_v = _mm512_add_epi64(slope_significand_v, slope_significand_v_tmp);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v, slope_exponent_v);
                    int32_v2 = _mm512_cvtepi64_epi32(slope_correct_v);

                    contact_int32_v = _mm512_inserti64x4(contact_int32_v, int32_v2, 1);
                    result_v = _mm512_add_epi32(intercept_v, contact_int32_v);

                    corrections_p += key_nums_tmp;
                    corrections_v = _mm512_load_epi32(corrections_p);
                    result_v = _mm512_add_epi32(result_v, corrections_v);
                    _mm512_store_epi32(result_int32, result_v);
                     for (Covered_Value k = 0; k < key_nums_tmp; k++)
                         output[++pointers] = result_int32[k];
                }

                for (Covered_Value k = 0; k < key_nums; k++) {
                    auto covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        auto slope_significand = slope_significand_p[k];
                        auto slope_exponent = slope_exponent_p[k];
                        auto intercept = intercept_p[k];
                        for (Covered_Value pos = cover_length_tmp; pos < covered; pos++)
                            output[++pointers] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
                    }
                }
            }
            for(auto it = segments_sort.begin() + idx; it < std::prev(segments_sort.end()); it++) {
                auto& s = *it;
                auto covered = s.covered;
                auto slope_significand = s.slope_significand;
                auto slope_exponent = s.slope_exponent;
                auto intercept = s.intercept;
                for (Covered_Value pos = 0; pos < covered; pos++)
                 output[++pointers] = ((slope_significand * pos) >> slope_exponent) + intercept + corrections_vector_residual[++correct_pointers];
             }
            auto end1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
            duration = duration1.count();
            return output;
        }

        std::vector<K> simd_decode_512i_simple(){
            std::vector<K> output;
            output.resize(n);
            Covered_Value pointers = -1;
            Covered_Value correct_pointers = -1;
            Covered_Value correct_residual_pointers = -1;
            Covered_Value covered_length_pointer = -1;
            Covered_Value key_nums_tmp = key_nums * 2;
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
                    Correction_Value *corrections_p = corrections_simd[++correct_pointers];
                    __m512i slope_correct_v, result_v, contact_int32_v, corrections_v;
                    __m256i int32_v1, int32_v2;

                    __m512i slope_significand_v = _mm512_set_epi64(7 * significand, 6 * significand, 5 * significand, 4 * significand, 3 * significand, 2 * significand, 1 * significand, 0);
                    __m512i slope_significand_v_add = _mm512_set1_epi64(8 * significand);
                    __m512i slope_exponent_v = _mm512_set1_epi64(exponent);
                    __m512i intercept_v = _mm512_set1_epi32(intercept);

                    for (Covered_Value i = 0; i < covered_length_tmp; i++) {
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
                        corrections_p += key_nums_tmp;
                    }

                }

                for (Covered_Value j = covered_length_tmp * key_nums_tmp; j < covered; ++j){
                    output[++pointers] = ((j * significand) >> exponent) + intercept + corrections_vector_residual[++correct_residual_pointers];
                }
            }

            auto end1 = std::chrono::high_resolution_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
            duration = duration1.count();
            return output;
        }

        static void set_correction(uint64_t* corrections_val, uint64_t n, uint64_t i, uint64_t value) {
            auto idx = get_correction_bit_offset(n, i);
            sdsl::bits::write_int(corrections_val + (idx >> 6), value, idx & 0x3F, BIT_WIDTH(Epsilon));
        }

        static uint64_t get_correction_bit_offset(uint64_t n, uint64_t i) {
            uint8_t bpc = BIT_WIDTH(Epsilon);
            if (i % extraction_density == 0)
                return bpc * (i / extraction_density);
            return bpc * (i + n / extraction_density - i / extraction_density);
        }

        static int64_t get_correction(const uint64_t* corrections_val, uint64_t n, uint64_t i, const sdsl::bit_vector& signs)  {
            auto idx = get_correction_bit_offset(n, i);
            uint64_t correction = sdsl::bits::read_int(corrections_val + (idx >> 6u), idx & 0x3F, BIT_WIDTH(Epsilon));
            return ((signs[i] == 0) ? correction : (-correction));
        }

        static int64_t get_correction(const uint64_t* corrections_val, uint64_t n, uint64_t i, const sdsl::rrr_vector<>& signs_rrr)  {
            auto idx = get_correction_bit_offset(n, i);
            uint64_t correction = sdsl::bits::read_int(corrections_val + (idx >> 6u), idx & 0x3F, BIT_WIDTH(Epsilon));
            return ((signs_rrr[i] == 0) ? correction : (-correction));
        }
     };

    #pragma pack(push, 1)

    // segment
    template <typename K, uint64_t Epsilon, uint64_t EpsilonRecursive, typename Floating>
    struct PGMIndex<K, Epsilon, EpsilonRecursive, Floating>::Segment {
        Covered_Value first;
        Intercept_Value intercept;
        uint8_t slope_exponent;
        uint64_t slope_significand;
        Covered_Value covered;

        Segment() = default;

        explicit Segment(position_type first):
            first(first),
            intercept(0),
            slope_exponent(0),
            slope_significand(0),
            covered(0) {} // covered is only used convenient for SIMD, actually there is no need to store it

        template <typename RandomIt>
        Segment(const canonical_segment& cs,
            uint64_t n,
            uint64_t i,
            uint64_t j,
            RandomIt data,
            uint64_t& errorPointCount,
            uint64_t* corrections_ptr,
            sdsl::bit_vector& signs,
            uint8_t bpc,
            position_type corrections_offset) {

            first = cs.get_first_x();
            if (first == n) {
                covered = 0;
                return;
            }

            auto [cs_significand, cs_exponent, cs_slope_float, cs_intercept] = cs.get_fixed_point_segment(first, j - i + 1); // 获取线段的斜率和截距，定点表示法
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max()) {
                throw std::overflow_error("Change the type of Segment::intercept to uint64");
            }
            if (cs_intercept < 0)
                throw std::overflow_error("Unexpected intercept < 0");

            slope_exponent = cs_exponent;
            slope_significand = cs_significand;
            intercept = cs_intercept;
            covered = (j - i);

            if (first == n - 1) {
                intercept = data[first];
                slope_exponent = 0;
                slope_significand = 0;
                auto error = static_cast<int64_t> (static_cast<int64_t> (data[first]) - approximate(first)); // __int128 ->  int64_t
                uint64_t correction_val = std::abs(error);
                set_correction(corrections_ptr, n, first, correction_val);
                signs[first] = 0;
                return;
            }

            for (auto p = first; p < j; p++) {
                auto error = static_cast<int64_t> (data[p]) - approximate(p); // __int128 ->  int64_t
                uint64_t correction = std::abs(error);
                set_correction(corrections_ptr, n, p, correction);
                bool sign = error > 0;
                signs[p] = sign ? 0 : 1;
                if (std::abs(error) > (Epsilon + 1)) {
                    errorPointCount++;
                }
            }
        }

        friend inline bool operator<(const Segment& s, const K& k) { return s.first < k; }
        friend inline bool operator<(const K& k, const Segment& s) { return k < s.first; }
        friend inline bool operator<(const Segment& s, const Segment& t) { return s.covered > t.covered; } // for sort

        operator K() { return first; };

        int64_t approximate(uint64_t i) const {
            return (int64_t(slope_significand * (i - first)) >> slope_exponent) + intercept;
        }

        inline int64_t operator()(const K& p) const { // 计算近似值
            return approximate(p);
        }
    };

    #pragma pack(pop)

}
