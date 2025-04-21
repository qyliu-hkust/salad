#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <vp2intersect.h>
#include <pgm_index.hpp>
#include <pgm_index_enumerate.hpp>
#include <vp2union.hpp>

#include "../external/mm_file/include/mm_file/mm_file.hpp"

namespace pgm_sequence {
    template <typename K, uint64_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t, Floating is unused
    class pgm_querier{
        using PGMIndexVariant = std::variant<
            pgm::PGMIndex<K, 1, 0, Floating>,
            pgm::PGMIndex<K, 3, 0, Floating>,
            pgm::PGMIndex<K, 7, 0, Floating>,
            pgm::PGMIndex<K, 15, 0, Floating>,
            pgm::PGMIndex<K, 31, 0, Floating>,
            pgm::PGMIndex<K, 63, 0, Floating>,
            pgm::PGMIndex<K, 127, 0, Floating>,
            pgm::PGMIndex<K, 255, 0, Floating>,
            pgm::PGMIndex<K, 511, 0, Floating>,
            pgm::PGMIndex<K, 1023, 0, Floating>,
            pgm::PGMIndex<K, 2047, 0, Floating>,
            pgm::PGMIndex<K, 4095, 0, Floating>>;
    protected:
        uint64_t data_size = 0;
        uint32_t query_num = 0;
        uint64_t query_time = 0;
        std::string input_basename = "";

        std::vector<std::vector<uint32_t>> read_query(const std::string& filename) {
            std::vector<std::vector<uint32_t>> idLists;
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return idLists; // Return an empty vector if the file could not be opened.
            }
            std::string line;
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                std::vector<uint32_t> ids;
                uint32_t id;
                while (iss >> id) {                 // Extract uint32_t from the line until no more can be found.
                    ids.push_back(id);
                }
                idLists.push_back(ids);
            }
            query_num = idLists[0].size();
            std::cout << "Total query sequences: " << idLists.size() << " Query num: " << query_num << std::endl;
            file.close();
            return idLists;
        }

        // std::vector<std::vector<K>> ori_data;
        // void load_orignal_data() {
        //     mm::file_source<K> input("/mnt/home/xyzhu/datasets/cw12b-1M/cw12b-1M.docs", mm::advice::sequential);
        //     K const* data = input.data();
        //     assert(data[0] == 1);
        //     std::cerr << "Universe Size: " << data[1] << std::endl;
        //     for (size_t i = 2; i < input.size();){
        //         uint64_t n = data[i];
        //         std::vector<K> sequence(data + i + 1, data + i + n);
        //         ori_data.push_back(sequence);
        //         data_size += n;
        //         i += n + 1;
        //     }
        //     input.close();
        // }

        std::vector<pgm_enumerator<K>> load_model(std::vector<uint32_t> idx_list) {
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return {};
            }
            std::vector<pgm_enumerator<K>> index_sequences;
            K index_num = 0;
            std::ifstream in(input_basename + "idx.size", std::ios::binary);
            in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            in.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));
            in.close();
            for (auto & index_count : idx_list) {
                if (index_count >= index_num) {
                    std::cerr << "Index count out of range: " << index_count << std::endl;
                    return {};
                }
                string input_filename = input_basename + std::to_string(index_count);
                std::ifstream in(input_filename + ".idx", std::ios::binary);
                uint64_t Epsilon_Data = 0;
                while (in.peek() != EOF) {
                    { // this {} is used to destroy index every loop obviously
                        PGMIndexVariant variant_index;
                        // epsilon
                        in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));
                        switch (Epsilon_Data) {
                            case 1: variant_index = pgm::PGMIndex<K, 1, 0, Floating>(); break;
                            case 3: variant_index = pgm::PGMIndex<K, 3, 0, Floating>(); break;
                            case 7: variant_index = pgm::PGMIndex<K, 7, 0, Floating>(); break;
                            case 15: variant_index = pgm::PGMIndex<K, 15, 0, Floating>(); break;
                            case 31: variant_index = pgm::PGMIndex<K, 31, 0, Floating>(); break;
                            case 63: variant_index = pgm::PGMIndex<K, 63, 0, Floating>(); break;
                            case 127: variant_index = pgm::PGMIndex<K, 127, 0, Floating>(); break;
                            case 255: variant_index = pgm::PGMIndex<K, 255, 0, Floating>(); break;
                            case 511: variant_index = pgm::PGMIndex<K, 511, 0, Floating>(); break;
                            case 1023: variant_index = pgm::PGMIndex<K, 1023, 0, Floating>(); break;
                            case 2047: variant_index = pgm::PGMIndex<K, 2047, 0, Floating>(); break;
                            case 4095: variant_index = pgm::PGMIndex<K, 4095, 0, Floating>(); break;
                            default: std::cerr << "Unsupported Epsilon Value: " << Epsilon_Data << std::endl; break;
                        }
                        std::visit([&in, &Epsilon_Data](auto &index) {
                            // Epsilon_Data
                            index.Epsilon_Data = Epsilon_Data;

                            //  n
                            in.read(reinterpret_cast<char*>(&index.n), sizeof(index.n));

                            //  first_pos
                            in.read(reinterpret_cast<char*>(&index.first_pos), sizeof(index.first_pos));

                            //  levels_offsets
                            size_t levels_offsets_size;
                            in.read(reinterpret_cast<char*>(&levels_offsets_size), sizeof(size_t));
                            index.levels_offsets.resize(levels_offsets_size);
                            in.read(reinterpret_cast<char*>(index.levels_offsets.data()), levels_offsets_size * sizeof(size_t));

                            //  segments
                            in.read(reinterpret_cast<char*>(&index.segments_size), sizeof(size_t));
                            in.read(reinterpret_cast<char*>(&index.bpc_first), sizeof(uint8_t));
                            in.read(reinterpret_cast<char*>(&index.bpc_covered), sizeof(uint8_t));
                            in.read(reinterpret_cast<char*>(&index.bpc_intercept), sizeof(uint8_t));
                            in.read(reinterpret_cast<char*>(&index.bpc_slope_exponent), sizeof(uint8_t));
                            in.read(reinterpret_cast<char*>(&index.bpc_slope_significand), sizeof(uint8_t));

                            size_t segments_size;
                            in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                            index.seg_first_compress.resize(segments_size);
                            in.read(reinterpret_cast<char*>(index.seg_first_compress.data()), segments_size * sizeof(uint64_t));
                            in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                            index.seg_covered_compress.resize(segments_size);
                            in.read(reinterpret_cast<char*>(index.seg_covered_compress.data()), segments_size * sizeof(uint64_t));
                            in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                            index.seg_intercept_compress.resize(segments_size);
                            in.read(reinterpret_cast<char*>(index.seg_intercept_compress.data()), segments_size * sizeof(uint64_t));
                            in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                            index.seg_slope_exponent_compress.resize(segments_size);
                            in.read(reinterpret_cast<char*>(index.seg_slope_exponent_compress.data()), segments_size * sizeof(uint64_t));
                            in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                            index.seg_slope_significand_compress.resize(segments_size);
                            in.read(reinterpret_cast<char*>(index.seg_slope_significand_compress.data()), segments_size * sizeof(uint64_t));

                            //  signs
                            size_t signs_size;
                            in.read(reinterpret_cast<char*>(&signs_size), sizeof(size_t));
                            index.signs_compress.resize(signs_size);
                            in.read(reinterpret_cast<char*>(index.signs_compress.data()), (signs_size + 7) / 8);
                            in.read(reinterpret_cast<char*>(&signs_size), sizeof(size_t));
                            index.signs_exception.resize(signs_size);
                            in.read(reinterpret_cast<char*>(index.signs_exception.data()), (signs_size + 7) / 8);

                            //  corrections
                            size_t corrections_size;
                            in.read(reinterpret_cast<char*>(&corrections_size), sizeof(size_t));
                            index.corrections_compress.resize(corrections_size);
                            in.read(reinterpret_cast<char*>(index.corrections_compress.data()), corrections_size * sizeof(uint64_t));
                            in.read(reinterpret_cast<char*>(&corrections_size), sizeof(size_t));
                            index.corrections_exception.resize(corrections_size);
                            in.read(reinterpret_cast<char*>(index.corrections_exception.data()), corrections_size * sizeof(uint64_t));
                        }, variant_index);
                        pgm_enumerator<K> enumerator_tmp;
                        std::visit([&enumerator_tmp](auto &index) {
                            index.normal_init();
                            for (int i = 0; i < index.segments_size; i++) {
                                auto first = index.seg_first[i];
                                auto intercept = index.seg_intercept[i];
                                auto slope_exponent = index.seg_slope_exponent[i];
                                auto slope_significand = index.seg_slope_significand[i];
                                auto covered = index.seg_covered[i];
                                enumerator_tmp.segments.emplace_back(first, intercept, slope_exponent, slope_significand, covered);
                            }
                            enumerator_tmp.load_copy(index.n, index.corrections_vector);
                        }, variant_index);
                        index_sequences.push_back(enumerator_tmp);
                    }
                }
                in.close();
            }
            return index_sequences;
        }

        uint32_t intersect_u32_normal(uint32_t const *a, uint32_t const *b,uint32_t a_length, uint32_t b_length, uint32_t *out) {
            uint32_t* const intersect_start = out;

            const uint32_t* a_end = a + a_length;
            const uint32_t* b_end = b + b_length;

            while (a < a_end && b < b_end) {
                if (*a < *b) {
                    ++a;
                } else if (*a > *b) {
                    ++b;
                } else {
                    // Found a match
                    *out++ = *a;
                    ++a;
                    ++b;
                }
            }
            return out - intersect_start;
        }

        uint32_t intersect_u32_simd_basic(uint32_t const* shorter, uint32_t const* longer, uint32_t shorter_length, uint32_t longer_length, uint32_t * out) {
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
                bool equal_count = equal_mask != 0;
                if (equal_count) {
                    out[intersection_count++] = shorter_member;
                }

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
            return intersection_count;
        }

        uint32_t intersect_u32_simd(uint32_t const* a, uint32_t const* b, uint32_t a_length, uint32_t b_length, uint32_t *out) {
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
                _mm512_mask_compressstoreu_epi32(out + c, a_matches, a_vec.zmm);


                c += _mm_popcnt_u32(a_matches); // The `_popcnt32` symbol isn't recognized by MSVC

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
            c += intersect_u32_normal(a, b, a_end - a, b_end - b, out + c);
            // result += c; // And merge it with the main body result
            return c;
        }

        // uint64_t union_u32_normal(const uint32_t *a, const uint32_t *b, size_t a_size, size_t b_size, uint32_t *out) {
        //     const uint32_t *a_end = a + a_size;
        //     const uint32_t *b_end = b + b_size;
        //     uint32_t *out_end = out;
        //     while (a != a_end && b != b_end) {
        //         bool le = (*a <= *b);
        //         bool ge = (*a >= *b);
        //         *out_end = le ? *a : *b;
        //         a += le;
        //         b += ge;
        //         out_end++;
        //     }
        //     // cerr << "Debug 1" << endl;
        //     std::memcpy(out_end, a, (a_end - a) * sizeof(uint32_t));
        //     out_end += a_end - a;
        //     std::memcpy(out_end, b, (b_end - b) * sizeof(uint32_t));
        //     out_end += b_end - b;
        //     return out_end - out;
        // }


        long double avg_skip = 0;
        long double avg_query_total_size = 0;
        long double avg_query_real_size = 0;
        void query_test_intersection(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            double avg_time = 0;
            K repeat_num = 1;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                for (auto &query : query_list) {
                    std::vector<pgm_enumerator<K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const pgm_enumerator<K> &a, const pgm_enumerator<K> &b) {return a.n < b.n;});

                    int branch_num = query.size() == 2 ? 2 : 3;

                    for (int i = 0; i < branch_num; i++) {
                        index_sequences[i].query_init(decode_type, "intersection");
                    }

                    uint32_t query_id_idx = 0;
                    uint32_t candidate_posting = 0;
                    uint32_t equal_result = 0;
                    uint32_t candidate_posting_tmp = 0;
                    int intersection_size = query.size() == 2 ? index_sequences[1].n + 1 : index_sequences[2].n + 1;
                    std::vector<K, HugePageAllocator<K>> intersection_result_1(intersection_size);
                    std::vector<K, HugePageAllocator<K>> intersection_result_2(intersection_size);
                    K *intersection_result_p1 = intersection_result_1.data(), *intersection_result_p2 = intersection_result_2.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < intersection_size; i++)
                            intersection_result_p2[i] = 0;
                        for (auto i = 0;i < intersection_size; i++)
                            intersection_result_p1[i] = 0;
                    }
                    // index_sequences[1].warm_up();
                    // index_sequences[0].warm_up();


                    auto start = std::chrono::high_resolution_clock::now();
                    if (query.size() == 2) {
                        // for (auto &enumerator : index_sequences) {
                        //     // enumerator.decode_query(decode_type);
                        // }
                        index_sequences[0].simd_decode_512i(intersection_result_p1);
                        index_sequences[1].simd_decode_512i(intersection_result_p2);
                        // equal_result = intersect_u32_simd(index_sequences[0].current_value_vector.data(), index_sequences[1].current_value_vector.data(), index_sequences[0].n, index_sequences[1].n, intercection_result.data());
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, index_sequences[0].n,  index_sequences[1].n, intersection_result_p1);
                    } else {
                        // for (int i = 0; i < 2; i++) {
                            // index_sequences[i].decode_query(decode_type);
                        // }
                        // equal_result = intersect_u32_simd(index_sequences[0].current_value_vector.data(), index_sequences[1].current_value_vector.data(), index_sequences[0].n, index_sequences[1].n, intercection_result.data());

                        index_sequences[0].simd_decode_512i(intersection_result_p1);
                        index_sequences[1].simd_decode_512i(intersection_result_p2);
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, index_sequences[0].n, index_sequences[1].n, intersection_result_p1);

                        index_sequences[2].simd_decode_512i(intersection_result_p2);
                        // equal_result = intersect_u32_simd(intercection_result.data(), index_sequences[2].current_value_vector.data(), equal_result, index_sequences[2].n, intercection_result.data());
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, equal_result, index_sequences[2].n, intersection_result_p1);

                        if (index_sequences.size() > branch_num) {
                            uint32_t intersection_length = equal_result;
                            equal_result = 0;
                            uint32_t intersection_idx = 0;
                            while (intersection_idx < intersection_length) {
                                candidate_posting = intersection_result_p1[intersection_idx];
                                for (query_id_idx = 3; query_id_idx < index_sequences.size(); query_id_idx++) {
                                    candidate_posting_tmp = index_sequences[query_id_idx].nextgeq(candidate_posting);
                                    if (candidate_posting_tmp != candidate_posting) {
                                        candidate_posting = candidate_posting_tmp;
                                        break;
                                    }
                                }
                                if (query_id_idx == index_sequences.size()) {
                                    intersection_result_p1[equal_result++] = candidate_posting;
                                }
                                intersection_idx++;
                            }
                        }
                    }
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

                    if (repeat == 0) {
                        total += equal_result;

                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;
                    }

                    if (repeat > 0)
                        avg_time += duration.count();
                    if (decode_type == "simd") {
                        for (auto &enumerator : index_sequences) {
                            enumerator.free_memory();
                        }
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
            }

            std::cerr << "Average query time: " <<  avg_time / query_list.size() / double(repeat_num) << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << endl;
        }

        void query_test_union(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            double avg_time = 0;
            K repeat_num = 1;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                for (auto &query : query_list) {
                    std::vector<pgm_enumerator<K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const pgm_enumerator<K> &a, const pgm_enumerator<K> &b) {return a.n < b.n;});

                    uint64_t union_size = 0, union_list_size = 0;
                    for (int i = 0; i < query.size(); i++) {
                        index_sequences[i].query_init(decode_type, "union");
                        union_size += index_sequences[i].n;
                        union_list_size = index_sequences[i].n > union_list_size ? index_sequences[i].n : union_list_size;
                    }

                    uint64_t equal_result = 0;
                    std::vector<K, HugePageAllocator<K>>  union_result_1(union_list_size);
                    std::vector<K, HugePageAllocator<K>>  union_result_2(union_size);
                    std::vector<K, HugePageAllocator<K>>  union_result_3(union_size);

                    K *union_result_p1 = union_result_1.data(), *union_result_p2 = union_result_2.data(), *union_result_p3 = union_result_3.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < union_size; i++)
                            union_result_p3[i] = 0;
                        for (auto i = 0;i < union_size; i++)
                            union_result_p2[i] = 0;
                        for (auto i = 0;i < union_list_size; i++)
                            union_result_p1[i] = 0;
                    }

                    auto start = std::chrono::high_resolution_clock::now();
                    if (query.size() == 2) {
                        // for (auto &enumerator : index_sequences) {
                        //     enumerator.decode_query(decode_type);
                        // }
                        // equal_result = union_u32_normal(index_sequences[0].current_value_vector.data(), index_sequences[1].current_value_vector.data(), index_sequences[0].n, index_sequences[1].n, union_result_p1);
                        // equal_result = union_u32_simd(index_sequences[0].current_value_vector.data(), index_sequences[1].current_value_vector.data(), index_sequences[0].n, index_sequences[1].n, union_result_p1);

                        // end_union_result = std::set_union(index_sequences[0].current_value_vector.data(), index_sequences[0].current_value_vector.data()+index_sequences[0].n, index_sequences[1].current_value_vector.data(), index_sequences[1].current_value_vector.data()+index_sequences[1].n, start_union_result);
                        // equal_result = end_union_result - start_union_result;
                        index_sequences[0].simd_decode_512i(union_result_p1);
                        index_sequences[1].simd_decode_512i(union_result_p2);
                        equal_result = union_u32_simd(union_result_p1, union_result_p2, index_sequences[0].n, index_sequences[1].n, union_result_p3);
                    } else {
                        // for (int i = 0; i < 2; i++) {
                            // index_sequences[i].decode_query(decode_type);
                        // }
                        // equal_result = union_u32_normal(index_sequences[0].current_value_vector.data(), index_sequences[1].current_value_vector.data(), index_sequences[0].n, index_sequences[1].n, union_result_p1);
                        // equal_result = union_u32_simd(index_sequences[0].current_value_vector.data(), index_sequences[1].current_value_vector.data(), index_sequences[0].n, index_sequences[1].n, union_result_p1);

                        // end_union_result = std::set_union(index_sequences[0].current_value_vector.data(), index_sequences[0].current_value_vector.data()+index_sequences[0].n, index_sequences[1].current_value_vector.data(), index_sequences[1].current_value_vector.data()+index_sequences[1].n, start_union_result);
                        // equal_result = end_union_result - start_union_result;
                        index_sequences[0].simd_decode_512i(union_result_p1);
                        index_sequences[1].simd_decode_512i(union_result_p2);
                        equal_result = union_u32_simd(union_result_p1, union_result_p2, index_sequences[0].n, index_sequences[1].n, union_result_p3);

                        for (int i = 2; i < query.size(); i++) {
                            index_sequences[i].simd_decode_512i(union_result_p1);
                            // index_sequences[i].decode_query(decode_type);
                            // equal_result = union_u32_normal(union_result_p1, index_sequences[i].current_value_vector.data(), equal_result, index_sequences[i].n, union_result_p2);
                            equal_result = union_u32_simd(union_result_p3, union_result_p1, equal_result, index_sequences[i].n, union_result_p2);
                            K *union_result_tmp = union_result_p2;
                            union_result_p2 = union_result_p3;
                            union_result_p3 = union_result_tmp;
                            // end_union_result = std::set_union(start_union_result, start_union_result + equal_result, index_sequences[i].current_value_vector.data(), index_sequences[i].current_value_vector.data()+index_sequences[i].n, start_union_result);
                            // equal_result = end_union_result - start_union_result;
                        }
                    }
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

                    if (repeat == 0) {
                        total += equal_result;

                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;
                    }

                    if (repeat > 0)
                        avg_time += duration.count();
                    if (decode_type == "simd") {
                        for (auto &enumerator : index_sequences) {
                            enumerator.free_memory();
                        }
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
            }

            std::cerr << "Average query time: " <<  avg_time / query_list.size() / double(repeat_num) << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << endl;
        }

    public:
        void test_query(const std::string input_filename, const std::string &decode_type, const std::string &query_filename, const std::string query_type) {
            std::vector<std::vector<uint32_t>> query_list = read_query(query_filename);
            input_basename = input_filename;
            if (query_type == "AND")
                query_test_intersection(query_list, decode_type);
            else if (query_type == "OR")
                query_test_union(query_list, decode_type);
        }
    };
}
