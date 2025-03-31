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

#include <climits>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <cmath>
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <tools.hpp>
#include <piecewise_linear_model.hpp>

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
        uint64_t Epsilon_Data;                ///< The epsilon value used to build the index.
        std::vector<uint64_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

        K first_pos;                     ///< The smallest element.
        std::vector<K> seg_first;
        std::vector<Intercept_Value> seg_intercept;
        std::vector<uint8_t> seg_slope_exponent;
        std::vector<int64_t> seg_slope_significand;
        std::vector<Covered_Value> seg_covered;
        sdsl::int_vector<64> seg_first_compress;
        sdsl::int_vector<64> seg_intercept_compress;
        sdsl::int_vector<64> seg_slope_exponent_compress;
        sdsl::int_vector<64> seg_slope_significand_compress;
        sdsl::int_vector<64> seg_covered_compress;
        std::vector<Segment> segments;      ///< The segments composing the index.


        static constexpr uint8_t bit_compress = Epsilon < 7 ? Epsilon < 3 ? 0 : 1 : 2;
        uint32_t epsilon_compress = Epsilon >> bit_compress;
        uint32_t bpc_compress = BIT_WIDTH(Epsilon) - bit_compress;
        uint32_t bpc_exception = BIT_WIDTH(Epsilon) + 1;
        sdsl::bit_vector signs_compress; // signs for saving
        sdsl::bit_vector signs_exception; // signs for saving
        uint64_t corrections_exception_num_write;
        uint64_t corrections_exception_num_read;
        sdsl::int_vector<64> corrections_compress; // corrections for saving, each value <= (epsilon / 4)
        sdsl::int_vector<64> corrections_exception; // corrections for saving, each (epsilon / 4) < value <= epsilon
        std::vector<Correction_Value> corrections_vector; // corrections for decode
        uint64_t segments_size;
        uint8_t bpc_first, bpc_covered, bpc_intercept, bpc_slope_significand, bpc_slope_exponent;

        /// Sentinel value to avoid bounds checking.
        static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity() : std::numeric_limits<K>::max();

        // static constexpr bool auto_bpc = false;
        // static constexpr uint64_t cache_line_bits = 64 * CHAR_BIT;
        // static constexpr uint64_t extraction_density = auto_bpc ? 1 : BIT_CEIL(4 * cache_line_bits / BIT_WIDTH(Epsilon));

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
            errorPointCount(0){
            Epsilon_Data = Epsilon;
            corrections_vector.resize(n);
            corrections_compress.resize(((n * bpc_compress + 63) >> 6)); // set 3 to avoid insEufficient memory
            signs_compress.resize(n);
            if (bit_compress > 0) {
                corrections_exception.resize(((n * bpc_exception + 63) >> 6));
                signs_exception.resize(n);
            }

            corrections_exception_num_write = 0;
            corrections_exception_num_read = 0;
            build(begin, end, Epsilon, EpsilonRecursive, segments, levels_offsets, errorPointCount, corrections_compress, corrections_exception,  signs_compress, signs_exception);

            if (bit_compress > 0) {
                corrections_exception.resize((corrections_exception_num_write * bpc_exception + 63) >> 6); // resize to the actual size
                signs_exception.resize(corrections_exception_num_write);
            }
            segments_compress();

            // varify_segments();
            // varify_corrections();
        }

        template <typename RandomIt>
        void build(RandomIt begin, RandomIt end,
                          uint64_t epsilon,
                          uint64_t epsilon_recursive, // we don't use this parameter to build the recursive PGM-index
                          std::vector<Segment>& segments,
                          std::vector<uint64_t>& levels_offsets,
                          uint64_t& errorPointCount,
                          sdsl::int_vector<64>& corrections_compress,
                          sdsl::int_vector<64>& corrections_exception,
                          sdsl::bit_vector& signs_compress,
                          sdsl::bit_vector& signs_exception) {

            auto n = (uint64_t) std::distance(begin, end);
            if (n == 0) return;

            levels_offsets.push_back(0);
            segments.reserve(n / (epsilon * epsilon));

            if (*std::prev(--end) == sentinel)
                throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

            std::vector<canonical_segment> canonical_segments;
            canonical_segments.reserve(epsilon > 0 ? n / (epsilon * epsilon) : n / 8);

            auto in_fun = [begin](auto i) { return std::pair<position_type, K>(i, begin[i]); };
            auto out_fun = [&canonical_segments](auto cs) { canonical_segments.push_back(cs); };
            auto n_segments = internal::make_segmentation_par(n, epsilon, in_fun, out_fun);
            auto last_n = n_segments;

            for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it) {
                auto i = it->get_first_x();
                auto j = std::next(it) != canonical_segments.end() ? std::next(it)->get_first_x() : n;
                build_segments(*it, n, i, j, begin, errorPointCount, corrections_compress.data(), corrections_exception.data(), signs_compress, signs_exception); // build the segment
                // segments.emplace_back(*it, n, i, j, begin, errorPointCount, corrections_compress.data(), corrections_exception.data(), signs_compress, signs_exception); // TODO: delete this line
            }
            levels_offsets.push_back(seg_first.size());
        }

        template <typename RandomIt>
        void build_segments(const canonical_segment& cs,
            uint64_t n, uint64_t i, uint64_t j,
            RandomIt data, uint64_t& errorPointCount,
            uint64_t* corrections_compress_ptr, uint64_t* corrections_exception_ptr,
            sdsl::bit_vector& signs_compress_ptr, sdsl::bit_vector& signs_exception_ptr) {

            uint32_t first = cs.get_first_x();
            if (first == n) return;

            seg_first.push_back(first);

            auto [cs_significand, cs_exponent, cs_intercept] = cs.get_fixed_point_segment(first, j - i + 1); // fixed point slope and intercept

            if (first == n - 1) {
                seg_intercept.push_back(data[first]);
                seg_slope_exponent.push_back(0);
                seg_slope_significand.push_back(0);
                seg_covered.push_back(1);
                set_correction(corrections_compress_ptr, corrections_exception_ptr, first, 0, signs_compress_ptr, signs_exception_ptr);
                return;
            }

            seg_slope_exponent.push_back(cs_exponent);
            seg_slope_significand.push_back(cs_significand);
            seg_intercept.push_back(cs_intercept);
            seg_covered.push_back(j - i);


            if (bit_compress > 0) {
                int64_t last_correction = 0;
                for (Covered_Value p = first; p < j; p++) {
                    int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                    corrections_vector[p] = error;
                    int64_t error_diff = error - last_correction;
                    set_correction(corrections_compress_ptr, corrections_exception_ptr, p, error_diff, signs_compress_ptr, signs_exception_ptr);
                    last_correction = error;
                }
            } else {
                for (Covered_Value p = first; p < j; p++) {
                    int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                    corrections_vector[p] = error;
                    uint8_t sign_value = error >= 0 ? 0 : 1;
                    error = std::abs(error);
                    if (error <= Epsilon)
                        set_correction_compress(corrections_compress_ptr, p, error, sign_value, signs_compress_ptr);
                    else if (error == Epsilon + 1 && sign_value == 0) // use -0 to represent error = epsilon
                        set_correction_compress(corrections_compress_ptr, p, 0, 1, signs_compress_ptr);
                    else
                        cerr << "Error: error = -epsilon: " << error << " " << sign_value << endl;
                }
            }
        }

        void segments_compress() {
            uint32_t max_fist = seg_first.back();
            // cerr << "max_fist: " << max_fist << endl;
            bpc_first = BIT_WIDTH(max_fist);
            auto max_iter_1 = std::max_element(seg_covered.begin(), seg_covered.end());
            bpc_covered = BIT_WIDTH(*max_iter_1);
            auto max_iter_2 = std::max_element(seg_intercept.begin(), seg_intercept.end());
            bpc_intercept = BIT_WIDTH(*max_iter_2);
            auto max_iter_3 = std::max_element(seg_slope_significand.begin(), seg_slope_significand.end());
            bpc_slope_significand = BIT_WIDTH(*max_iter_3);
            auto max_iter_4 = std::max_element(seg_slope_exponent.begin(), seg_slope_exponent.end());
            bpc_slope_exponent = BIT_WIDTH(*max_iter_4);

            segments_size = seg_first.size();
            seg_first_compress.resize((segments_size * bpc_first + 63) >> 6);
            seg_covered_compress.resize((segments_size * bpc_covered + 63) >> 6);
            seg_intercept_compress.resize((segments_size * bpc_intercept + 63) >> 6);
            seg_slope_significand_compress.resize((segments_size * bpc_slope_significand + 63) >> 6);
            seg_slope_exponent_compress.resize((segments_size * bpc_slope_exponent + 63) >> 6);

            for (int i = 0; i < seg_first.size(); i++) {
                auto first = seg_first[i];
                set_segment_compress(seg_first_compress.data(), i, first, bpc_first);
                auto covered = seg_covered[i];
                set_segment_compress(seg_covered_compress.data(), i, covered, bpc_covered);
                auto intercept = seg_intercept[i];
                set_segment_compress(seg_intercept_compress.data(), i, intercept, bpc_intercept);
                auto slope_significand = seg_slope_significand[i];
                set_segment_compress(seg_slope_significand_compress.data(), i, slope_significand, bpc_slope_significand);
                auto slope_exponent = seg_slope_exponent[i];
                set_segment_compress(seg_slope_exponent_compress.data(), i, slope_exponent, bpc_slope_exponent);
            }
        }

        void varify_segments() {
            if (segments_size != segments.size()) {
                cerr << "Error: seg_first.size != segments.size" << seg_first.size() << " " << segments.size() << endl;
            }
            for (int i = 0; i < segments_size; i++) {
                auto first = get_segment_compress(seg_first_compress.data(), i, bpc_first);
                auto covered = get_segment_compress(seg_covered_compress.data(), i, bpc_covered);
                auto intercept = get_segment_compress(seg_intercept_compress.data(), i, bpc_intercept);
                auto slope_significand = get_segment_compress(seg_slope_significand_compress.data(), i, bpc_slope_significand);
                auto slope_exponent = get_segment_compress(seg_slope_exponent_compress.data(), i, bpc_slope_exponent);
                auto segment = segments[i];
                if (first != segment.first) {
                    cerr << "Error: first != segment.first: " << first << " " << segment.first << endl;
                }
                if (covered != segment.covered) {
                    cerr << "Error: covered != segment.covered: " << covered << " " << segment.covered << endl;
                }
                if (intercept != segment.intercept) {
                    cerr << "Error: intercept != segment.intercept: " << intercept << " " << segment.intercept << endl;
                }
                if (slope_significand != segment.slope_significand) {
                    cerr << "Error: slope_significand != segment.slope_significand: " << slope_significand << " " << segment.slope_significand << endl;
                }
                if (slope_exponent != segment.slope_exponent) {
                    cerr << "Error: slope_exponent != segment.slope_exponent: " << slope_exponent << " " << segment.slope_exponent << endl;
                }
            }
        }

        void varify_corrections() {
            corrections_exception_num_read = 0;
            if (bit_compress > 0) {
                for (int i = 0; i < seg_first.size(); i++) {
                    auto first = seg_first[i];
                    auto covered = seg_covered[i];
                    int32_t last_correction = 0;
                    for (int j = first; j < covered + first; j++) {
                        auto error1 = corrections_vector[j];
                        auto error2 = last_correction + get_correction(corrections_compress.data(), corrections_exception.data(), j, signs_compress, signs_exception);
                        last_correction = error2;
                        if (error1 != error2) {
                            cerr << "Error 1: error1 != error2: " << error1 << " " << error2 << " " << epsilon_compress << endl;
                        }
                    }
                }
            } else {
                for (int i = 0; i < seg_first.size(); i++) {
                    auto first = seg_first[i];
                    auto covered = seg_covered[i];
                    for (int j = first; j < covered + first; j++) {
                        auto error1 = corrections_vector[j];
                        auto error2 = get_correction_uncompress(corrections_compress.data(), j, signs_compress);
                        if (error1 != error2) {
                            cerr << "Error 2: error1 != error2: " << error1 << " " << error2 << " " << epsilon_compress << endl;
                        }
                    }
                }
            }
        }

        int64_t seg_approximate(uint32_t i, uint32_t first, uint8_t exponent, int64_t significand, int32_t intercept) const {
            return (int64_t(significand * (i - first)) >> exponent) + intercept;
        }

        uint64_t segments_count() const { return segments.empty() ? 0 : levels_offsets[1] - 1; }

        uint64_t height() const { return levels_offsets.size() - 1; }

        // std::vector<Segment> get_segments() const { return segments; }

        // auto segment_for_pos(const K& pos) const {
        //     return std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), pos));
        // }

        uint64_t segment_slope_significand_max() const {
            auto max_iter = std::max_element(seg_slope_significand.begin(), seg_slope_significand.end());
            uint64_t max_slope_significand = *max_iter;
            // for (auto it = segments.begin(); it < std::prev(segments.end()); ++it)
            //     if (it->slope_significand > max_slope_significand)
            //         max_slope_significand = it->slope_significand;
            return max_slope_significand;
        }

        uint32_t segment_slope_exponent_max() const {
            auto max_iter = std::max_element(seg_slope_exponent.begin(), seg_slope_exponent.end());
            uint32_t max_slope_exponent = *max_iter;
            // for (auto it = segments.begin(); it < std::prev(segments.end()); ++it)
            //     if (it->slope_exponent > max_slope_exponent)
            //         max_slope_exponent = it->slope_exponent;
            return max_slope_exponent;
        }

        uint64_t size() const {
            return n;
        }

        uint64_t total_size_in_bytes() const {
            return segment_size_in_bytes() + corrections_size_in_bytes() + signs_size_in_bytes();
        }

        uint64_t segment_size_in_bytes() const {
            // return segments_count() * (sizeof(Covered_Value) + sizeof(Intercept_Value) + sizeof(uint64_t) + sizeof(uint8_t)); // first + Intercept + slope_significant + slope_exponent
            return (seg_first_compress.bit_size() + seg_intercept_compress.bit_size() + seg_slope_significand_compress.bit_size() + seg_slope_exponent_compress.bit_size()) / CHAR_BIT;
        }

        uint64_t corrections_size_in_bytes() const {
            return (corrections_compress.bit_size() + corrections_exception.bit_size()) / CHAR_BIT;
        }

        uint64_t signs_size_in_bytes() const {
            return (signs_compress.bit_size() + signs_exception.bit_size()) / CHAR_BIT;
        }

        void reoprt_residual_random_segment(std::ofstream &file, int seg_idx) {
            if (bit_compress > 0) {
                corrections_exception_num_read = 0;
                for (int i = 0; i < segments_size; i++) {
                    if (i == seg_idx) {
                        auto first = seg_first[i];
                        auto covered = seg_covered[i];
                        int32_t last_correction = 0;
                        file << covered << std::endl;
                        for (int j = first; j < covered + first; j++) {
                            auto error1 = get_correction(corrections_compress.data(), corrections_exception.data(), j, signs_compress, signs_exception);
                            auto error2 = last_correction + error1;
                            last_correction = error2;
                            auto flag = std::abs(error1) <= epsilon_compress ? 0 : 1;
                            file << error2 << "\t" << error1 << "\t"<< flag << std::endl;
                        }
                    }
                }
            } else {
                throw std::invalid_argument("Not support uncompress residuals");
            }
        }

        // SP_Tree<int> first_pos_index;
        // void build_first_pos_index(const std::string& output_basename) {
        //     std::vector<int> first_pos_vector;
        //     for (auto it = segments.begin(); it < std::prev(segments.end()); ++it) {
        //         int first_value = it->intercept + get_correction(corrections.data(), n, it -> first, signs);
        //         first_pos_vector.push_back(first_value);
        //     }
        //     first_pos_index.build(first_pos_vector);
        //     first_pos_index.save_tree(output_basename);
        // }

        // void load_first_pos_index(const std::string& input_basename) {
        //     first_pos_index.load_tree(input_basename);
        // }

        void huffman_init() {
            corrections_vector.resize(n, 0);
            for(auto &s : segments) {
                auto covered = s.covered;
                auto first = s.first;
                int32_t last_correction = 0;
                for (Covered_Value j = 0; j < covered; ++j) {
                    int32_t correction_diff = get_correction(corrections_compress.data(), corrections_exception.data(), j, signs_compress, signs_exception);
                    corrections_vector[j] = last_correction + correction_diff;
                    last_correction = corrections_vector[j];
                }
            }
        }

        void normal_init(){
            seg_first.resize(segments_size);
            seg_covered.resize(segments_size);
            seg_intercept.resize(segments_size);
            seg_slope_significand.resize(segments_size);
            seg_slope_exponent.resize(segments_size);

            for (int i = 0; i < segments_size; i++) {
                seg_first[i] = get_segment_compress(seg_first_compress.data(), i, bpc_first);
                seg_covered[i] = get_segment_compress(seg_covered_compress.data(), i, bpc_covered);
                seg_intercept[i] = get_segment_compress(seg_intercept_compress.data(), i, bpc_intercept);
                seg_slope_significand[i] = get_segment_compress(seg_slope_significand_compress.data(), i, bpc_slope_significand);
                seg_slope_exponent[i] = get_segment_compress(seg_slope_exponent_compress.data(), i, bpc_slope_exponent);
            }

            corrections_vector.resize(n);
            corrections_exception_num_read = 0;
            if (bit_compress > 0) {
                for(int i = 0; i < seg_first.size(); i++) {
                    auto covered = seg_covered[i];
                    auto first = seg_first[i];
                    // int32_t last_correction = 0;
                    for (Covered_Value j = first; j < first + covered; ++j) {
                        int32_t correction_varify = get_correction(corrections_compress.data(), corrections_exception.data(), j, signs_compress, signs_exception);
                        // if (correction_varify != corrections_vector[j]) {
                        //     cerr << "Error: correction_varify != corrections_vector[j]: " << correction_varify << " " << corrections_vector[j] << " " << epsilon_compress << endl;
                        // }
                        corrections_vector[j] = correction_varify;
                        // last_correction = correction_varify;
                    }
                }
            } else {
                for(int i = 0; i < seg_first.size(); i++) {
                    auto covered = seg_covered[i];
                    auto first = seg_first[i];
                    for (Covered_Value j = first; j < first + covered; ++j) {
                        int32_t correction_varify = get_correction_uncompress(corrections_compress.data(), j, signs_compress);
                        // if (correction_varify != corrections_vector[j]) {
                        //     cerr << "Error: correction_varify != corrections_vector[j]: " << correction_varify << " " << corrections_vector[j] << " " << epsilon_compress << endl;
                        // }
                        corrections_vector[j] = correction_varify;
                    }
                }
            }


            // for(auto &s : segments) {
            //     auto covered = s.covered;
            //     auto first = s.first;
            //     int32_t last_correction = 0;
            //     for (Covered_Value j = first; j < first + covered; ++j) {
            //         int32_t correction_diff = get_correction(corrections_compress.data(), corrections_exception.data(), j, signs_compress, signs_exception);
            //         corrections_vector[j] = last_correction + correction_diff;
            //         last_correction = corrections_vector[j];
            //     }
            // }
        }

        std::vector<K> normal_decode() {
            std::vector<K> output;
            output.resize(n);
            uint32_t pointer = 0;
            for (int i= 0; i < seg_first.size(); i++) {
                auto first = seg_first[i];
                auto covered = seg_covered[i];
                auto significand = seg_slope_significand[i];
                auto exponent = seg_slope_exponent[i];
                auto intercept = seg_intercept[i];

                int32_t last_correction = intercept;
                for (int j = 0; j < covered; ++j) {
                    last_correction = last_correction + corrections_vector[pointer];
                    output[pointer] = ((j * significand) >> exponent) + last_correction;
                    pointer++;
                }
            }
            // for (auto it = segments.begin(); it != std::prev(segments.end()); ++it) {
            //     auto& s = *it;
            //     auto covered = s.covered;
            //     auto significand = s.slope_significand;
            //     auto exponent = s.slope_exponent;
            //     auto intercept = s.intercept;
            //     auto first = s.first;
            //     for (Covered_Value j = 0; j < covered; ++j) {
            //         output[j + first] = ((j * significand) >> exponent) + intercept + corrections_vector[j +first];
            //     }
            // }
            return output;
        }

        uint64_t get_correction_bit_offset(uint64_t i, uint32_t bpc) {
            return i * bpc;
        }

        void set_segment_compress(uint64_t* compress_val, uint64_t i, uint64_t value, uint8_t bpc) {
            uint64_t idx = get_correction_bit_offset(i, bpc);
            sdsl::bits::write_int(compress_val + (idx >> 6u), value, idx & 0x3f, bpc);
        }

        int64_t get_segment_compress(const uint64_t* compress_val, uint64_t i, uint8_t bpc) {
            uint64_t idx = get_correction_bit_offset(i, bpc);
            uint64_t val = sdsl::bits::read_int(compress_val + (idx >> 6u), idx & 0x3f, bpc);
            return val;
        }

        void set_correction_compress(uint64_t* corrections_compress_val, uint64_t i, uint64_t value, uint8_t sign_value, sdsl::bit_vector& signs_compress_ptr) {
            uint64_t idx = get_correction_bit_offset(i, bpc_compress);
            sdsl::bits::write_int(corrections_compress_val + (idx >> 6u), value, idx & 0x3f, bpc_compress);
            signs_compress_ptr[i] = sign_value & 0x01;
        }

        void set_correction_exception(uint64_t* corrections_exception_val, uint64_t value, uint8_t sign_value, sdsl::bit_vector& signs_exception_ptr) {
            uint64_t idx = get_correction_bit_offset(corrections_exception_num_write, bpc_exception);
            sdsl::bits::write_int(corrections_exception_val + (idx >> 6u), value, idx & 0x3f, bpc_exception);
            signs_exception_ptr[corrections_exception_num_write++] = sign_value & 0x01;
        }

        void set_correction(uint64_t* corrections_compress_val, uint64_t* corrections_exception_val, uint64_t i, int64_t value, sdsl::bit_vector& signs_compress_ptr, sdsl::bit_vector& signs_exception_ptr) {
            uint64_t abs_value = std::abs(value);
            if (abs_value <= epsilon_compress) {
                set_correction_compress(corrections_compress_val, i, abs_value, value >= 0 ? 0 : 1, signs_compress_ptr);
            } else {
                set_correction_exception(corrections_exception_val, abs_value, value >= 0 ? 0 : 1, signs_exception_ptr);
                set_correction_compress(corrections_compress_val, i, 0, 1, signs_compress_ptr);
            }
        }

        int64_t get_correction(const uint64_t* corrections_compress_val, const uint64_t* corrections_exception_val, int64_t i, const sdsl::bit_vector& signs_compress_ptr, const sdsl::bit_vector& signs_exception_ptr)  {
            uint64_t idx = get_correction_bit_offset(i, bpc_compress);
            uint64_t correction = sdsl::bits::read_int(corrections_compress_val + (idx >> 6u), idx & 0x3f, bpc_compress);
            if (correction == 0 && signs_compress_ptr[i] == 1) {
                idx = get_correction_bit_offset(corrections_exception_num_read, bpc_exception);
                correction = sdsl::bits::read_int(corrections_exception_val + (idx >> 6u), idx & 0x3f, bpc_exception);
                return signs_exception_ptr[corrections_exception_num_read++] == 0 ? correction : -correction;
            }
            return signs_compress_ptr[i] == 0 ? correction : -correction;
        }

        int64_t get_correction_uncompress(const uint64_t* corrections_compress_val, int64_t i, const sdsl::bit_vector& signs_compress_ptr)  {
            uint64_t idx = get_correction_bit_offset(i, bpc_compress);
            uint64_t correction = sdsl::bits::read_int(corrections_compress_val + (idx >> 6u), idx & 0x3f, bpc_compress);
            if (correction == 0 && signs_compress_ptr[i] == 1)    return Epsilon + 1;
            return signs_compress_ptr[i] == 0 ? correction : -correction;
        }

     };

    #pragma pack(push, 1)

    // segment
    template <typename K, uint64_t Epsilon, uint64_t EpsilonRecursive, typename Floating>
    struct PGMIndex<K, Epsilon, EpsilonRecursive, Floating>::Segment {
        K first;
        Intercept_Value intercept; // 32 bits
        uint8_t slope_exponent;
        uint64_t slope_significand;
        Covered_Value covered; // covered is not necessary, just convenient for SIMD

        Segment() = default;

        explicit Segment(position_type first):
            first(first),
            intercept(0),
            slope_exponent(0),
            slope_significand(0),
            covered(0) {} // covered is only used convenient for SIMD decode, actually there is no need to store it

        template <typename RandomIt>
        Segment(const canonical_segment& cs,
            uint64_t n,
            uint64_t i,
            uint64_t j,
            RandomIt data,
            uint64_t& errorPointCount,
            uint64_t* corrections_compress_ptr,
            uint64_t* corrections_exception_ptr,
            sdsl::bit_vector& signs_compress_ptr,
            sdsl::bit_vector& signs_exception_ptr) {

            first = cs.get_first_x();
            if (first == n) {
                covered = 0;
                intercept = 0;
                slope_exponent = 0;
                slope_significand = 0;
                return;
            }

            auto [cs_significand, cs_exponent, cs_intercept] = cs.get_fixed_point_segment(first, j - i + 1); // 获取线段的斜率和截距，定点表示法
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max()) {
                throw std::overflow_error("Change the type of Segment::intercept to uint64");
            }
            else if(cs_intercept < 0)
                throw std::overflow_error("Unexpected intercept < 0");

            slope_exponent = cs_exponent;
            slope_significand = cs_significand;
            intercept = cs_intercept;
            covered = (j - i);

            if (first == n - 1) {
                intercept = data[first];
                slope_exponent = 0;
                slope_significand = 0;
                // set_correction(corrections_compress_ptr, corrections_exception_ptr, first, 0, signs_compress_ptr, signs_exception_ptr);
                return;
            }

            // int64_t last_correction = static_cast<int64_t> (data[first]) - approximate(first);
            // set_correction(corrections_compress_ptr, corrections_exception_ptr, first, last_correction, signs_compress_ptr, signs_exception_ptr);

            // for (K p = first + 1; p < j; p++) {
            //     int64_t error = static_cast<int64_t> (data[p]) - approximate(p);
            //     int64_t error_diff = error - last_correction;
            //     set_correction(corrections_compress_ptr, corrections_exception_ptr, p, error_diff, signs_compress_ptr, signs_exception_ptr);
            //     last_correction = error;
            //     if (error_diff >= (Epsilon << 1)) {
            //         errorPointCount++;
            //     }
            // }
        }

        friend inline bool operator<(const Segment& s, const K& k) { return s.first < k; }
        friend inline bool operator<(const K& k, const Segment& s) { return k < s.first; }
        friend inline bool operator<(const Segment& s, const Segment& t) { return s.covered > t.covered; } // for sort

        operator K() { return first; };

        Intercept_Value get_intercept() const { return intercept; }

        int64_t approximate(uint64_t i) const {
            return (int64_t(slope_significand * (i - first)) >> slope_exponent) + intercept;
        }

        inline int64_t operator()(const K& p) const { // 计算近似值
            return approximate(p);
        }
    };

    #pragma pack(pop)

}