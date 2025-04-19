#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <pgm_index.hpp>
#include <pgm_index_enumerate.hpp>

namespace pgm_sequence {
    template <typename K, size_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t
    class pgm_decoder{
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
        uint64_t total_unequal = 0;
        std::vector<std::string> perf_header;
        std::vector<std::vector<double>> perf_data;

        void decode_test(pgm_enumerator<K> &index, std::string decode_type, bool warm_up = true) {
            total_list++;
            uint64_t unequal = 0;

            if (decode_type == "simd") {
                index.simd_init();

                // std::vector<K> result1(index.n);
                std::vector<K, HugePageAllocator<K>> result1(index.n);
                if (warm_up) {
                    for (int warm_time = 0; warm_time < 5; warm_time++) {
                        for (int i = 0; i < index.n; i ++) {
                            result1[i] = 0;
                        }
                    }
                }
                PerfEvent perf_event;
                perf_event.startCounters();
                index.simd_decode_512i(result1.data());
                perf_event.stopCounters();
                std::stringstream header_out;
                std::stringstream data_out;
                PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
                perf_event.printReport(header_out, data_out, 1);
                std::vector<double> perf_data_tmp = parseDoubles(data_out.str(), ',');
                for (int i = 0; i < perf_data_tmp.size(); i++) {
                    perf_data[i].push_back(perf_data_tmp[i]);
                }

                // std::vector<K> result2(index.n);
                // index.normal_decode(result2.data());
                // for (int i = 0; i < result1.size(); i++) {
                //     // std::cerr << "Error: result1 != result2: " << result1[i] << " " << result2[i] << " " << i << " " << result1.size()  << " " << int(result1[i]) - int(result2[i]) << std::endl;
                //     if (result1[i] != result2[i]) {
                //         unequal++;
                //         // std::cerr << "Error: result1 != result2: " << result1[i] << " " << result2[i] << " " << i << " " << result1.size()  << " " << int(result1[i]) - int(result2[i]) << std::endl;
                //     }
                // }

                std::vector<K, HugePageAllocator<K>> ().swap(result1);
                total_unequal += unequal;
                total_calculated += index.total_calculated;
                total_calculated_add += index.total_calculated_add;
                total_conversion_time += index.conversion_time;

                size_tmp = index.n;
                per_integer_time_tmp1 = index.total_duration;
                per_integer_time_tmp1 /= size_tmp;
                per_integer_time1 += per_integer_time_tmp1 * (size_tmp / data_size);

                avg_decode_time1 += index.total_duration;
                if (index.total_duration > max_decode_time1)
                    max_decode_time1 = index.total_duration;
                if (index.total_duration < min_decode_time1)
                    min_decode_time1 = index.total_duration;
                index.free_memory(decode_type);
            } else if (decode_type == "normal") {
                std::vector<K, HugePageAllocator<K>> result2(index.n);
                // std::vector<K> result2(index.n);
                if (warm_up) {
                    for (int warm_time = 0; warm_time < 5; warm_time++) {
                        for (int i = 0; i < index.n; i ++) {
                            result2[i] = 0;
                        }
                    }
                }
                PerfEvent perf_event;
                perf_event.startCounters();
                index.normal_decode(result2.data());
                perf_event.stopCounters();
                std::stringstream header_out;
                std::stringstream data_out;
                PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
                perf_event.printReport(header_out, data_out, 1);
                std::vector<double> perf_data_tmp = parseDoubles(data_out.str(), ',');
                for (int i = 0; i < perf_data_tmp.size(); i++) {
                    perf_data[i].push_back(perf_data_tmp[i]);
                }

                std::vector<K, HugePageAllocator<K>> ().swap(result2);
                // std::vector<K> ().swap(result2);
                size_tmp = index.n;
                per_integer_time_tmp2 = index.total_duration;
                per_integer_time_tmp2 /= size_tmp;
                per_integer_time2 += per_integer_time_tmp2 * (size_tmp / data_size);

                avg_decode_time2 += index.total_duration;
                if (index.total_duration > max_decode_time2)
                    max_decode_time2 = index.total_duration;
                if (index.total_duration < min_decode_time2)
                    min_decode_time2 = index.total_duration;
            }
        }

        void result_statistic(std::string decode_type) {
            std::cerr << "Total list: " << total_list << std::endl;
            if (decode_type == "simd" || decode_type == "simd_simple") {
                avg_decode_time1 -= max_decode_time1;
                avg_decode_time1 -= min_decode_time1;
                std::cerr << "Decode time 1 simd, average: " << avg_decode_time1 / (total_list - 2)  << ", max: " << max_decode_time1 << ", min: " << min_decode_time1 << " microseconds" <<std::endl;
                std::cerr << "Decode per integer: " << per_integer_time1 * 1000.0 << " nanoseconds" << std::endl;
                std::cerr << "Total calculated: " << total_calculated << " (" << double(total_calculated) / data_size << ") , Total calculated add: " << total_calculated_add << ", Total Conversion time: " << total_conversion_time << ", Total Unequal: " << total_unequal << std::endl;
            }
            if (decode_type == "normal") {
                avg_decode_time2 -= max_decode_time2;
                avg_decode_time2 -= min_decode_time2;
                std::cerr << "Decode time 2 normal, average: " << avg_decode_time2 / (total_list - 2)  << ", max: " << max_decode_time2 << ", min: " << min_decode_time2 << " microseconds" <<std::endl;
                std::cerr << "Decode per integer: " << per_integer_time2 * 1000.0 << " nanoseconds" << std::endl;
            }

            std::cerr << "Performance: " << std::endl;
            for (const auto &header : perf_header) {
                std::cerr << std::setw(15) << header; // 每个字段占 15 个字符宽度
            }
            std::cerr << std::endl;

            // 输出均值
            for (const auto &data : perf_data) {
                double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
                std::cerr << std::setw(15) << std::fixed << std::setprecision(2) << mean; // 定宽输出，保留两位小数
            }
            std::cerr << std::endl << std::endl;

        }

        void test_model(std::string input_basename, std::string decode_type) {
            // init perf_header
            PerfEvent perf_event;
            perf_event.startCounters();
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }
            perf_event.stopCounters();
            std::stringstream header_out;
            std::stringstream data_out;
            PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
            perf_event.printReport(header_out, data_out, 1);
            perf_header = split_str(header_out.str(), ',');
            perf_data.resize(perf_header.size());

            std::cerr << "Load index from: " << input_basename << std::endl;
            std::cerr << "Decode type: " << decode_type << std::endl;
            K index_num = 0;
            std::ifstream in(input_basename + "idx.size", std::ios::binary);
            in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            in.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));
            in.close();
            for (K index_count = 0; index_count < index_num; index_count++) {
                std::ifstream in(input_basename + std::to_string(index_count) + ".idx", std::ios::binary);
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
                        // std::vector<K> result1 = ori_data[index_count];
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
                        // test PGMIndex
                        decode_test(enumerator_tmp, decode_type);
                    }
                }
                in.close();
            }
            result_statistic(decode_type);
        }
    };
}
