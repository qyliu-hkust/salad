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
// #include <typeinfo>
#include <cmath>
#include "sdsl/bits.hpp"
#include "sdsl/int_vector.hpp"
#include "tools.hpp"
#include "piecewise_linear_model.hpp"

namespace pgm
{
    #define BIT_CEIL(x) ((x) < 2 ? 1u : 1u << (64u - __builtin_clzll((x) - 1)))
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

    template <typename K, size_t Epsilon = 64, size_t EpsilonRecursive = 0, typename Floating = float>
    class PGMIndex
    {
        struct Segment;
        struct constant_bpc;
        struct variable_bpc;

        static constexpr bool auto_bpc = false;
        static constexpr size_t cache_line_bits = 64 * CHAR_BIT;
        static constexpr size_t extraction_density = auto_bpc ? 1 : BIT_CEIL(4 * cache_line_bits / BIT_WIDTH(Epsilon));

        using position_type = typename std::conditional_t<sizeof(K) <= 4, uint32_t, uint64_t>;
        using larger_signed_key_type = typename std::conditional_t<sizeof(K) <= 4, int64_t, __int128>;
        using canonical_segment = typename internal::OptimalPiecewiseLinearModel<position_type, K>::CanonicalSegment;

    public:
        // protected:
        template <typename, size_t, size_t, uint8_t, typename>
        friend class BucketingPGMIndex;

        template <typename, size_t, typename>
        friend class EliasFanoPGMIndex;

        static_assert(Epsilon > 0);

        size_t n; ///< The number of elements this index was built on.
        K first_pos; ///< The smallest element.
        std::vector<Segment> segments; ///< The segments composing the index.
        std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

        /// 存储误差数据
        sdsl::bit_vector signs;
        sdsl::int_vector<64> corrections;

        /// Sentinel value to avoid bounds checking.
        static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity() : std::numeric_limits<K>::max();


        template <typename RandomIt>
        static void build(RandomIt begin, RandomIt end,
                          size_t epsilon, size_t epsilon_recursive,
                          std::vector<Segment>& segments,
                          std::vector<size_t>& levels_offsets,
                          size_t& errorPointCount,
                          sdsl::int_vector<64>& corrections,
                          sdsl::bit_vector& signs) {

            auto n = (size_t)std::distance(begin, end) + 1;
            if (n == 0)
                return;

            levels_offsets.push_back(0);
            segments.reserve(n / (epsilon * epsilon));

            if (*std::prev(--end) == sentinel)
                throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

            size_t corrections_offset = 0;

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
                auto nextSegFirstPos = std::next(it) != canonical_segments.end() ? std::next(it)->get_first_x() : n;
                segments.emplace_back(*it, n, i, j, begin, errorPointCount, corrections.data(), signs, bpc, corrections_offset);
                corrections_offset += bpc * (j - i);
            }
            // segments contains every first pos , interception and slope
            segments.emplace_back(n);
            levels_offsets.push_back(segments.size());
        }


        // public:

        static constexpr size_t epsilon_value = Epsilon;

        size_t errorPointCount;

        PGMIndex(): n(0), errorPointCount(0) {
        }

        explicit PGMIndex(const std::vector<K>& data) : PGMIndex(data.begin(), data.end()) {
        }

        template <typename RandomIt>
        PGMIndex(RandomIt begin, RandomIt end):
            n(std::distance(begin, end) + 1),
            first_pos(n ? *begin : K(0)),
            segments(),
            levels_offsets(),
            errorPointCount(0),
            corrections(n * BIT_WIDTH(Epsilon) / 64),
            signs(n) {
            build(begin, end, Epsilon, EpsilonRecursive, segments, levels_offsets, errorPointCount, corrections, signs);
        }


        size_t segments_count() const { return segments.empty() ? 0 : levels_offsets[1] - 1; }

        size_t height() const { return levels_offsets.size() - 1; }

        std::vector<Segment> get_segments() const { return segments; }

        auto segment_for_pos(const K& pos) const {
            return std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), pos));
        }

        larger_signed_key_type search(const K& pos) const {
            auto p = std::max(first_pos, pos);
            auto it = segment_for_pos(p);
            auto approx = (*it)(p);
            return approx;
        }

        size_t indexSegments_bytes() const {
            return segments.size() * sizeof(Segment) + levels_offsets.size() * sizeof(size_t);
        }

        void decode(K* out) const {
            for (auto it = segments.begin(); it != std::prev(segments.end()); ++it) {
                auto& s = *it;
                auto covered = std::next(it)->first - s.first;
                auto significand = s.slope_significand;
                auto exponent = s.slope_exponent;
                auto intercept = s.intercept;
                #pragma omp simd
                for (size_t j = 0; j < covered; ++j)
                    out[j] = ((j * significand) >> exponent) + intercept;

                for (size_t j = 0; j < covered; ++j)
                    out[j] += s.get_correction(corrections.data(), n, j + s.first, signs);
                out += covered;
            }
        }

        // summary the counts
        std::vector<K> decode() const {
            std::vector<K> out(n);
            decode(out.data());
            return out;
        }

        size_t segment_size_in_bytes() const {
            return segments_count() * sizeof(Segment);
        }

        size_t corrections_size_in_bytes() const {
            return corrections.bit_size() / CHAR_BIT;
        }

        size_t signs_size_in_bytes() const {
            return signs.bit_size() / CHAR_BIT;
        }
    };

#pragma pack(push, 1)


    template <typename K, size_t Epsilon, size_t EpsilonRecursive, typename Floating>
    struct PGMIndex<K, Epsilon, EpsilonRecursive, Floating>::Segment {
        uint32_t first;
        K intercept;
        uint8_t slope_exponent;
        // uint64_t slope_significand;
        larger_signed_key_type slope_significand;

        Segment(){
        }

        explicit Segment(position_type first):
            first(first),
            intercept(std::numeric_limits<decltype(intercept)>::max()),
            slope_exponent(0),
            slope_significand(0) {
        }

        // construct segment：get slope and intercept
        template <typename RandomIt>
        Segment(const canonical_segment& cs, size_t n, size_t i, size_t j, RandomIt data, size_t& errorPointCount,
                uint64_t* corrections_ptr, sdsl::bit_vector& signs, uint8_t bpc, position_type corrections_offset) {
            first = cs.get_first_x();

            // uint64_t, uint8_t, SY=__int128
            // cout << "first" << first << " i and j" << i << " " << j << endl;
            auto [cs_significand, cs_exponent, cs_intercept] = cs.get_fixed_point_segment(first, j - i + 1);

            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max()) {
                throw std::overflow_error("Change the type of Segment::intercept to uint64");
            }
            if (cs_intercept < 0)
                throw std::overflow_error("Unexpected intercept < 0");

            intercept = cs_intercept;
            slope_exponent = cs_exponent;
            slope_significand = cs_significand;

            for (auto p = first; p < j; p++) {
                auto error = static_cast<int64_t>(data[p] - approximate(p)); // __int128 ->  int64_t
                uint64_t correction = std::abs(error);
                set_correction(corrections_ptr, n, p, correction);
                bool sign = error > 0;
                signs[p] = sign ? 0 : 1;
                if (abs(error) > (Epsilon + 1))
                {
                    errorPointCount++;
                }
            }
        }

        void set_correction(uint64_t* corrections, size_t n, size_t i, uint64_t value) {
            auto idx = get_correction_bit_offset(n, i);
            sdsl::bits::write_int(corrections + (idx >> 6), value, idx & 0x3F, BIT_WIDTH(Epsilon));
        }

        size_t get_correction_bit_offset(size_t n, size_t i) const {
            uint8_t bpc = BIT_WIDTH(Epsilon);
            if (i % extraction_density == 0)
                return bpc * (i / extraction_density);
            return bpc * (i + n / extraction_density - i / extraction_density);
        }

        int64_t get_correction(const uint64_t* corrections, size_t n, size_t i, const sdsl::bit_vector& signs) const {
            auto idx = get_correction_bit_offset(n, i);
            uint64_t correction = sdsl::bits::read_int(corrections + (idx >> 6u), idx & 0x3F, BIT_WIDTH(Epsilon));
            return ((signs[i] == 0) ? correction : (-correction));
        }

        friend inline bool operator<(const Segment& s, const K& k) { return s.first < k; }
        friend inline bool operator<(const K& k, const Segment& s) { return k < s.first; }
        friend inline bool operator<(const Segment& s, const Segment& t) { return s.first < t.first; }

        operator K() { return first; };

        larger_signed_key_type approximate(size_t i) const {
            return ((larger_signed_key_type(slope_significand) * (i - first)) >> slope_exponent) + intercept;
        }

        inline larger_signed_key_type operator()(const K& p) const {
            return ((larger_signed_key_type(slope_significand) * (p - first)) >> slope_exponent) + intercept;
        }
    };

#pragma pack(pop)
}
