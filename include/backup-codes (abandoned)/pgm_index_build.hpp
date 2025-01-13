#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"
#include "../external/mm_file/include/mm_file/mm_file.hpp"

namespace pgm_sequence {
    template <typename K, size_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t
    class pgm_indexes{
    public:
        std::vector<pgm::PGMIndex<K, epsilon, 0, Floating>> index_sequences;
        // std::vector<std::vector<K>> data_sequences;
        // std::vector<unordered_map<K, K>> maps_sequences;
        std::vector<std::vector<uint32_t>> maps_sequences;
        // std::vector<std::vector<K>> decode_results;
        uint64_t data_size = 0;
        uint64_t data_unequal = 0;
        uint32_t exp_rounds = 10;
        bool normal_init_flag = false;

        void build_model(std::string input_basename) {
            std::cerr << std::endl << "Epsilon: " << epsilon << std::endl;
            std::cerr << "Read File: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            std::cerr << "Universe Size: " << data[1] << std::endl;
            for (size_t i = 2; i < input.size();){
                uint64_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n);
                // data_sequences.emplace_back(sequence);
                index_sequences.push_back(pgm::PGMIndex<K, epsilon, 0, Floating>(sequence));
                data_size += n;
                i += n + 1;
            }
            input.close();
            normal_init(); // 初始化残差，转为int64
            data_test(input_basename); // 测试数据是否正确
        }

        void normal_init() {
            // if (normal_init_flag)
            //     return;
            std::cerr << "Normal Decode Init" << std::endl;
            // decode_results.resize(index_sequences.size());
            // for (K i = 0; i < index_sequences.size(); i++)
                // decode_results[i].resize(index_sequences[i].n);

            for (auto& index : index_sequences) {
                index.decode_init();
            }
            // normal_init_flag = true;
        }

        void simd_init(bool use_max=true) {
            std::cerr << "Simd init whether use max_length: " << (use_max == 0 ? "false" :"true") << std::endl;
            // if (!normal_init_flag)
            // decode_results.resize(index_sequences.size());
            // for (K i = 0; i < index_sequences.size(); i++)
                // decode_results[i].resize(index_sequences[i].n);

            // for (auto& index : index_sequences)
                // index.decode_init();
            for (auto& index : index_sequences) {
//                index.decode_init();
                index.simd_init_align(use_max);
                index.create_correction();
                // index.corrections_residual_init();
                // maps_sequences.emplace_back(index.result_map_init());
            }
//            for (auto& index : index_sequences) {
//                maps_sequences.emplace_back(index.result_map_init());
//            }

        }

        void simd_free() {
            std::cerr << "Free Memory " << std::endl;
            for (auto& index : index_sequences) {
                index.free_memory();
            }
        }

        void simd_decode() {
            // simd_init(false);
            std::cerr << "Simd Decode" << std::endl;
            double avg_decode_time2 = 0;
            uint64_t max_decode_time2 = 0;
            uint64_t min_decode_time2 = UINT64_MAX;
            // double calculate_rate = 0;
            // data_unequal = 0;
            for (int tim = 0; tim < exp_rounds + 1; tim++) {
                uint64_t total_decode_time2 = 0;
                for (K i = 0; i < index_sequences.size(); i++) {
                    auto& index = index_sequences[i];
                    // auto start2 = std::chrono::high_resolution_clock::now();
                    std::vector<K> result2 = index.simd_decode_512i();
                    // auto end2 = std::chrono::high_resolution_clock::now();
                    // auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
                    total_decode_time2 += index.duration;

                    // std::vector<K> result1 = index.decode();
                    // auto& result_map = maps_sequences[i];
                    // if (tim == 0) {
                    // // calculate_rate += index.total_calculated;
                    //     if (result1.size() != result2.size())
                    //     std::cerr << "ERROR: decode size error " << result2.size() << " " << result1.size() << std::endl;
                    //     for (auto j = 0; j < result2.size(); j++) {
                    //         if (result1[j] != result2[result_map[j]]) {
                    //         // std::cerr << "ERROR: decode error " << result1[j] << " " << result2[j] << std::endl;
                    //             data_unequal++;
                    //         }
                    //     }
                    // }
                }
                if (tim != 0) {
                    avg_decode_time2 += total_decode_time2;
                    if (total_decode_time2 > max_decode_time2)
                        max_decode_time2 = total_decode_time2;
                    if (total_decode_time2 < min_decode_time2)
                        min_decode_time2 = total_decode_time2;
                    cerr << "Simd test round " << tim << ", decode time simd: " << total_decode_time2 << " microseconds" << std::endl;
                }
            }
            avg_decode_time2 -= max_decode_time2;
            avg_decode_time2 -= min_decode_time2;
            std::cerr << "Decode time 2 simd, average: " << avg_decode_time2 / (exp_rounds - 2)  << ", max: " << max_decode_time2 << ", min: " << min_decode_time2 << " microseconds" <<std::endl;
            std::cerr << "Decode per integer: " << avg_decode_time2 / (exp_rounds - 2) * 1000.0 / data_size << " nanoseconds" << std::endl;
            // std::cerr << "Calculate rate: " << calculate_rate << " " << data_size << " " << calculate_rate / data_size << std::endl;
            // std::cerr << "Unequal postings: " << data_unequal << std::endl;
            simd_free();
        }

        void normal_decode() {
            // normal_init();
            std::cerr << "Normal Decode" << std::endl;
            double avg_decode_time1 = 0;
            uint64_t max_decode_time1 = 0;
            uint64_t min_decode_time1 = UINT64_MAX;
            // data_unequal = 0;
            for (int tim = 0; tim < exp_rounds + 1; tim++) {
                uint64_t total_decode_time1 = 0;
                for (K i = 0; i < index_sequences.size(); i++) {
                    auto& index = index_sequences[i];
                    // std::vector<K> result1 (index.n);
                    // auto start1 = std::chrono::high_resolution_clock::now();
                    // index.decode(decode_results[i]);
                    std::vector<K> result1 = index.decode();
                    // auto end1 = std::chrono::high_resolution_clock::now();
                    // auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
                    // total_decode_time1 += duration1.count();
                    total_decode_time1 += index.duration;
                }
                if (tim != 0) {
                    avg_decode_time1 += total_decode_time1;
                    if (total_decode_time1 > max_decode_time1)
                        max_decode_time1 = total_decode_time1;
                    if (total_decode_time1 < min_decode_time1)
                        min_decode_time1 = total_decode_time1;
                    std::cerr << "Normal test round " << tim << ", decode time single_floop: " << total_decode_time1 << " microseconds" << std::endl;
                }
            }
            avg_decode_time1 -= max_decode_time1;
            avg_decode_time1 -= min_decode_time1;
            std::cerr << "Decode time 1 normal, average: " << avg_decode_time1 / (exp_rounds - 2)  << ", max: " << max_decode_time1 << ", min: " << min_decode_time1 << " microseconds" <<std::endl;
            std::cerr << "Decode per integer: " << avg_decode_time1 / (exp_rounds - 2) * 1000.0 / data_size << " nanoseconds" << std::endl;
        }

        void statistic() {
            double segments_count = 0;
            uint64_t segments_size = 0;
            uint64_t corrections_size = 0;
            uint64_t signs_size = 0;
            uint64_t errorpoint_size = 0;
            double avg_covered = 0;
            for (auto& index : index_sequences) {
                for (auto& segment : index.segments) {
                    avg_covered += segment.covered;
                }
                segments_count += index.segments.size();
                segments_size += index.segment_size_in_bytes();
                corrections_size += index.corrections_size_in_bytes();
                signs_size += index.signs_size_in_bytes();
                errorpoint_size += index.errorPointCount;
            }
            double ratio = (segments_size + corrections_size + signs_size) / double(data_size) / double(sizeof(K));
            std::cerr << "Epsilon: " << epsilon << std::endl;
            std::cerr << "Integer Count: " << data_size << std::endl;
            std::cerr << "Segments count: " << segments_count << std::endl;
            std::cerr << "Error Point Count: " << errorpoint_size << std::endl;
            std::cerr << "Decode Unequal: " << data_unequal << std::endl;
            std::cerr << "Average Covered: " << avg_covered / segments_count << std::endl;
            std::cerr << "Average Length: " << segments_count / index_sequences.size() << std::endl;
            std::cerr << "Segment Size: " << segments_size << " byte" << std::endl;
            std::cerr << "Corrections Size: " << corrections_size << " byte" << std::endl;
            std::cerr << "Signs Size: " << signs_size << " byte" << std::endl;
            std::cerr << "Compression Ratio: " << ratio << std::endl;
            long double total_size_in_bytes = segments_size + corrections_size + signs_size;
            long double total_size_in_gib = total_size_in_bytes / 1024.0 / 1024.0 / 1024.0;
            std::cerr << "Total Size: " << total_size_in_bytes << " in bytes, " << total_size_in_gib << " in GiB, " << total_size_in_bytes / data_size << " bytes per int " << std::endl;
        }

        void save_covered(const std::string output_filename) {
            std::ofstream file(output_filename + ".covered.txt");
            if (file.is_open()) {
                for (K i = 0; i < index_sequences.size(); i++) {
                    auto index = index_sequences[i];
                    for (auto& segment : index.segments) {
                        file << segment.covered << std::endl;
                    }
                }
                file.close();
            }
            else {
                std::cerr << "Couldn't open the log_file" << std::endl;
            }
        }

        void save_residual(const std::string output_filename) {
            std::ofstream file(output_filename + ".residual.txt");
            // 随机取生成8个数
            std::vector<K> random_index;
            for (K i = 0; i < 8; i++) {
                random_index.push_back(rand() % index_sequences.size());
            }

            if (file.is_open()) {
                for (auto i : random_index) {
                    auto& index = index_sequences[i];
                    for (auto& correction : index.corrections) {
                        file << correction << std::endl;
                    }
                }
                file.close();
            }
            else {
                std::cerr << "Couldn't open the log_file" << std::endl;
            }
        }

        void data_test(std::string input_basename) {
            std::cerr << std::endl << "Epsilon: " << epsilon << std::endl;
            std::cerr << "Read File: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::cerr << "Universe Size: " << data[1] << std::endl;
            data_unequal = 0;
            K posi = 0;
            for (size_t i = 2; i < input.size();) {
                K n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n);
                auto& index = index_sequences[posi++];
                std::vector<K> result1 = index.decode();
                if (result1.size() != sequence.size())
                    std::cerr << "ERROR: decode size error " << result1.size() << " " << sequence.size() << std::endl;
                for (auto j = 0; j < result1.size(); j++) {
                    if (sequence[j] != result1[j]) {
                        // std::cerr << "ERROR: decode error " << sequence[j] << " " << result1[j] << " " << int32_t(result1[j]) - int32_t(sequence[j]) << " " << j << " " << result1.size() << std::endl;
                        data_unequal++;
                    }
                }
                i += n + 1;
            }
            std::cerr << "Unequal postings: " << data_unequal << std::endl;
        }

        void save_model(const std::string output_filename) {
            std::cerr << "Save index to: " << output_filename << std::endl;
            std::ofstream out(output_filename, std::ios::binary);
            out.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
            for (auto& index : index_sequences) {
                out.write(reinterpret_cast<const char*>(&index.n), sizeof(index.n)); // size_t n;
                out.write(reinterpret_cast<const char*>(&index.errorPointCount), sizeof(size_t));
                out.write(reinterpret_cast<const char*>(&index.first_pos), sizeof(index.first_pos)); // K first_pos;
                size_t segments_size = index.segments.size();
                out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                out.write(reinterpret_cast<const char*>(index.segments.data()), index.segments.size() * sizeof(typename pgm::PGMIndex<K, epsilon, 0, Floating>::Segment));
                size_t levels_offsets_size = index.levels_offsets.size();
                out.write(reinterpret_cast<const char*>(&levels_offsets_size), sizeof(size_t));
                out.write(reinterpret_cast<const char*>(index.levels_offsets.data()), index.levels_offsets.size() * sizeof(size_t));
                size_t signs_size = index.signs.size();
                out.write(reinterpret_cast<const char*>(&signs_size), sizeof(size_t)); // 写入大小
                out.write(reinterpret_cast<const char*>(index.signs.data()), (index.signs.size() + 7) / 8); // 写入数据
                size_t corrections_size = index.corrections.size();
                out.write(reinterpret_cast<const char*>(&corrections_size), sizeof(size_t)); // 写入大小
                out.write(reinterpret_cast<const char*>(index.corrections.data()), index.corrections.size() * sizeof(uint64_t)); // 写入数据
            }
            out.close();
        }

        void load_model(std::string input_filename) {
            std::cerr << "Load index from: " << input_filename << std::endl;
            std::ifstream in(input_filename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Failed to open file for reading: " + input_filename);
            }
            index_sequences.clear(); // 清空现有数据
            in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            while (in.peek() != EOF) {
                pgm::PGMIndex<K, epsilon, 0, Floating> index;
                in.read(reinterpret_cast<char*>(&index.n), sizeof(index.n));
                in.read(reinterpret_cast<char*>(&index.errorPointCount), sizeof(size_t));
                in.read(reinterpret_cast<char*>(&index.first_pos), sizeof(index.first_pos));
                size_t segments_size;
                in.read(reinterpret_cast<char*>(&segments_size), sizeof(size_t));
                index.segments.resize(segments_size);
                in.read(reinterpret_cast<char*>(index.segments.data()), segments_size * sizeof(typename pgm::PGMIndex<K, epsilon, 0, Floating>::Segment));
                size_t levels_offsets_size;
                in.read(reinterpret_cast<char*>(&levels_offsets_size), sizeof(size_t));
                index.levels_offsets.resize(levels_offsets_size);
                in.read(reinterpret_cast<char*>(index.levels_offsets.data()), levels_offsets_size * sizeof(size_t));
                size_t signs_size;
                in.read(reinterpret_cast<char*>(&signs_size), sizeof(size_t));
                index.signs.resize(signs_size);
                in.read(reinterpret_cast<char*>(index.signs.data()), (signs_size + 7) / 8);
                size_t corrections_size;
                in.read(reinterpret_cast<char*>(&corrections_size), sizeof(size_t));
                index.corrections.resize(corrections_size);
                in.read(reinterpret_cast<char*>(index.corrections.data()), corrections_size * sizeof(uint64_t));
                index_sequences.push_back(std::move(index));
            }
            in.close();
        }
    };
}
