#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"

namespace pgm_sequence {
    template <typename K, size_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t
    class pgm_decoder{
    public:
        uint64_t data_size = 0;
        long double avg_decode_time1 = 0;
        long double avg_decode_time2 = 0;
        long double per_integer_time1 = 0;
        long double per_integer_time2 = 0;
        long double per_integer_time_tmp1 = 0;
        long double per_integer_time_tmp2 = 0;
        long double size_tmp = 0;
        uint64_t total_list = 0;
        uint64_t max_decode_time1 = 0;
        uint64_t max_decode_time2 = 0;
        uint64_t min_decode_time1 = UINT64_MAX - 1;
        uint64_t min_decode_time2 = UINT64_MAX - 1;

        void decode_test(pgm::PGMIndex<K, epsilon, 0, Floating> &index, std::string decode_type) {
            total_list++;
            if (decode_type == "simd" || decode_type == "simd_simple" || decode_type == "all") {
                if (decode_type == "simd" || "all") {
                    index.simd_init(false);
                    index.create_corrections();
                    index.create_corrections_residual();
                    std::vector<K> result1 = index.simd_decode_512i();
                } else {
                    index.create_corrections_simple();
                    index.create_corrections_residual_simple();
                    std::vector<K> result1 = index.simd_decode_512i_simple();
                }
                size_tmp = index.n;
                per_integer_time_tmp1 = index.duration;
                per_integer_time_tmp1 /= size_tmp;
                per_integer_time1 += per_integer_time_tmp1 * (size_tmp / data_size);

                avg_decode_time1 += index.duration;
                if (index.duration > max_decode_time1)
                    max_decode_time1 = index.duration;
                if (index.duration < min_decode_time1)
                    min_decode_time1 = index.duration;
                index.free_memory();
            }

            if (decode_type == "normal" || decode_type == "all") {
                index.normal_init();
                std::vector<K> result2 = index.normal_decode();

                size_tmp = index.n;
                per_integer_time_tmp2 = index.duration;
                per_integer_time_tmp2 /= size_tmp;
                per_integer_time2 += per_integer_time_tmp2 * (size_tmp / data_size);

                avg_decode_time2 += index.duration;
                if (index.duration > max_decode_time2)
                    max_decode_time2 = index.duration;
                if (index.duration < min_decode_time2)
                    min_decode_time2 = index.duration;
            }
        }

        void result_statistic(std::string decode_type) {
            std::cerr << "Total list: " << total_list << std::endl;
            if (decode_type == "simd" || decode_type == "all" || decode_type == "simd_simple") {
                avg_decode_time1 -= max_decode_time1;
                avg_decode_time1 -= min_decode_time1;
                std::cerr << "Decode time 1 simd, average: " << avg_decode_time1 / (total_list - 2)  << ", max: " << max_decode_time1 << ", min: " << min_decode_time1 << " microseconds" <<std::endl;
                std::cerr << "Decode per integer: " << per_integer_time1 * 1000.0 << " nanoseconds" << std::endl;
            }
            if (decode_type == "normal" || decode_type == "all") {
                avg_decode_time2 -= max_decode_time2;
                avg_decode_time2 -= min_decode_time2;
                std::cerr << "Decode time 2 normal, average: " << avg_decode_time2 / (total_list - 2)  << ", max: " << max_decode_time2 << ", min: " << min_decode_time2 << " microseconds" <<std::endl;
                std::cerr << "Decode per integer: " << per_integer_time2 * 1000.0 << " nanoseconds" << std::endl;
            }
            std::cerr << std::endl;
        }

        void test_model(std::string input_filename, std::string decode_type) {
            std::cerr << "Load index from: " << input_filename << std::endl;
            std::ifstream in(input_filename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Failed to open file for reading: " + input_filename);
            }
            std::cerr << "Decode type: " << decode_type << std::endl;
            in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            while (in.peek() != EOF) {
                { // this {} is used to destroy index every loop obviously
                    pgm::PGMIndex<K, epsilon, 0, Floating> index;
                    in.read(reinterpret_cast<char*>(&index.n), sizeof(index.n));
                    in.read(reinterpret_cast<char*>(&index.errorPointCount), sizeof(size_t));
                    in.read(reinterpret_cast<char*>(&index.first_pos), sizeof(index.first_pos));
                    // read levels_offsets
                    size_t levels_offsets_size;
                    in.read(reinterpret_cast<char*>(&levels_offsets_size), sizeof(size_t));
                    index.levels_offsets.resize(levels_offsets_size);
                    in.read(reinterpret_cast<char*>(index.levels_offsets.data()), levels_offsets_size * sizeof(size_t));
                    // read segments
                    size_t segments_size;
                    in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                    index.segments.resize(segments_size);
                    in.read(reinterpret_cast<char*>(index.segments.data()), segments_size * sizeof(typename pgm::PGMIndex<K, epsilon, 0, Floating>::Segment));
                    // read signs
                    size_t signs_size;
                    in.read(reinterpret_cast<char*>(&signs_size), sizeof(size_t));
                    index.signs.resize(signs_size);
                    in.read(reinterpret_cast<char*>(index.signs.data()), (signs_size + 7) / 8);
                    // index.signs_rrr.load(in);
                    // read corrections
                    size_t corrections_size;
                    in.read(reinterpret_cast<char*>(&corrections_size), sizeof(size_t));
                    index.corrections.resize(corrections_size);
                    in.read(reinterpret_cast<char*>(index.corrections.data()), corrections_size * sizeof(uint64_t));
                    // test PGMIndex
                    decode_test(index, decode_type);
                }
            }
            in.close();
            result_statistic(decode_type);
        }
    };
}
