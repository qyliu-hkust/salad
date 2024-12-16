#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"

namespace pgm_sequence {
    template <typename K, size_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t
    class pgm_decoder{
        using PGMIndexVariant = std::variant<
            pgm::PGMIndex<K, 1, 0, Floating>,
            pgm::PGMIndex<K, 3, 0, Floating>,
            pgm::PGMIndex<K, 7, 0, Floating>,
            pgm::PGMIndex<K, 15, 0, Floating>,
            pgm::PGMIndex<K, 16, 0, Floating>,
            pgm::PGMIndex<K, 31, 0, Floating>,
            pgm::PGMIndex<K, 32, 0, Floating>,
            pgm::PGMIndex<K, 63, 0, Floating>,
            pgm::PGMIndex<K, 64, 0, Floating>,
            pgm::PGMIndex<K, 126, 0, Floating>,
            pgm::PGMIndex<K, 127, 0, Floating>,
            pgm::PGMIndex<K, 128, 0, Floating>,
            pgm::PGMIndex<K, 255, 0, Floating>,
            pgm::PGMIndex<K, 256, 0, Floating>,
            pgm::PGMIndex<K, 511, 0, Floating>,
            pgm::PGMIndex<K, 512, 0, Floating>,
            pgm::PGMIndex<K, 1023, 0, Floating>,
            pgm::PGMIndex<K, 1024, 0, Floating>,
            pgm::PGMIndex<K, 2047, 0, Floating>,
            pgm::PGMIndex<K, 2048, 0, Floating>,
            pgm::PGMIndex<K, 4095, 0, Floating>,
            pgm::PGMIndex<K, 4096, 0, Floating>,
            pgm::PGMIndex<K, 8191, 0, Floating>,
            pgm::PGMIndex<K, 8192, 0, Floating>,
            pgm::PGMIndex<K, 16383, 0, Floating>,
            pgm::PGMIndex<K, 16384, 0, Floating>,
            pgm::PGMIndex<K, 32767, 0, Floating>,
            pgm::PGMIndex<K, 32768, 0, Floating>,
            pgm::PGMIndex<K, 65535, 0, Floating>,
            pgm::PGMIndex<K, 65536, 0, Floating>,
            pgm::PGMIndex<K, 131071, 0, Floating>,
            pgm::PGMIndex<K, 131072, 0, Floating>,
            pgm::PGMIndex<K, 262143, 0, Floating>,
            pgm::PGMIndex<K, 262144, 0, Floating>>;
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
        uint64_t total_calculated = 0;
        uint64_t total_calculated_add = 0;
        uint64_t total_conversion_time = 0;

        void decode_test(PGMIndexVariant &variant_index, std::string decode_type) {
            total_list++;
            std::visit([&](auto &index) {
                if (decode_type == "simd" || decode_type == "simd_simple") {
                    if (decode_type == "simd") {
                        index.simd_init(false);
                        index.create_corrections();
                        index.create_corrections_residual();
                        std::vector<K> result1 = index.simd_decode_512i();
                    } else {
                        index.create_corrections_simple();
                        index.create_corrections_residual_simple();
                        std::vector<K> result1 = index.simd_decode_512i_simple();
                    }
                    total_calculated += index.total_calculated;
                    total_calculated_add += index.total_calculated_add;
                    total_conversion_time += index.conversion_time;

                    size_tmp = index.n;
                    per_integer_time_tmp1 = index.duration;
                    per_integer_time_tmp1 /= size_tmp;
                    per_integer_time1 += per_integer_time_tmp1 * (size_tmp / data_size);

                    avg_decode_time1 += index.duration;
                    if (index.duration > max_decode_time1)
                        max_decode_time1 = index.duration;
                    if (index.duration < min_decode_time1)
                        min_decode_time1 = index.duration;
                    index.free_memory(decode_type);
                }

                if (decode_type == "normal") {
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
            }, variant_index);

        }

        void result_statistic(std::string decode_type) {
            std::cerr << "Total list: " << total_list << std::endl;
            if (decode_type == "simd" || decode_type == "simd_simple") {
                avg_decode_time1 -= max_decode_time1;
                avg_decode_time1 -= min_decode_time1;
                std::cerr << "Decode time 1 simd, average: " << avg_decode_time1 / (total_list - 2)  << ", max: " << max_decode_time1 << ", min: " << min_decode_time1 << " microseconds" <<std::endl;
                std::cerr << "Decode per integer: " << per_integer_time1 * 1000.0 << " nanoseconds" << std::endl;
                std::cerr << "Total calculated: " << total_calculated << ", Total calculated add: " << total_calculated_add << ", Total Conversion time: " << total_conversion_time << std::endl;
            }
            if (decode_type == "normal") {
                avg_decode_time2 -= max_decode_time2;
                avg_decode_time2 -= min_decode_time2;
                std::cerr << "Decode time 2 normal, average: " << avg_decode_time2 / (total_list - 2)  << ", max: " << max_decode_time2 << ", min: " << min_decode_time2 << " microseconds" <<std::endl;
                std::cerr << "Decode per integer: " << per_integer_time2 * 1000.0 << " nanoseconds" << std::endl;
            }
            std::cerr << std::endl;
        }

        void test_model(std::string input_basename, std::string decode_type) {
            std::cerr << "Load index from: " << input_basename << std::endl;
            std::ifstream in(input_basename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Failed to open file for reading: " + input_basename);
            }
            std::cerr << "Decode type: " << decode_type << std::endl;
            in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            uint64_t Epsilon_Data = 0;
            while (in.peek() != EOF) {
                {
                    // this {} is used to destroy index every loop obviously
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
                        size_t segments_size;
                        in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                        index.segments.resize(segments_size);
                        in.read(reinterpret_cast<char*>(index.segments.data()), segments_size * sizeof(typename pgm::PGMIndex<K, epsilon, 0, Floating>::Segment));
                        //  signs
                        size_t signs_size;
                        in.read(reinterpret_cast<char*>(&signs_size), sizeof(size_t));
                        index.signs.resize(signs_size);
                        in.read(reinterpret_cast<char*>(index.signs.data()), (signs_size + 7) / 8);
                        //  corrections
                        size_t corrections_size;
                        in.read(reinterpret_cast<char*>(&corrections_size), sizeof(size_t));
                        index.corrections.resize(corrections_size);
                        in.read(reinterpret_cast<char*>(index.corrections.data()), corrections_size * sizeof(uint64_t));
                    }, variant_index);
                    // test PGMIndex
                    decode_test(variant_index, decode_type);
                }
            }
            in.close();
            result_statistic(decode_type);
        }
    };
}
