#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"
#include "huffman_encode.hpp"
#include "../external/mm_file/include/mm_file/mm_file.hpp"

namespace pgm_sequence {
    template <typename K, size_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t
    class pgm_builder{
    public:
        std::vector<pgm::PGMIndex<K, epsilon, 0, Floating>> index_sequences;
        uint64_t data_size = 0;
        uint64_t data_unequal = 0;

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
            normal_init();
            data_test(input_basename); // test data
        }

        void normal_init() {
            // std::cerr << "Data Test Init" << std::endl;
            for (auto& index : index_sequences) {
                index.normal_init();
            }
        }

        void statistic(std::string output_filename="") {
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
                // signs_size += index.signs_rrr_size_in_bytes();
                errorpoint_size += index.errorPointCount;
            }

            double ratio = (segments_size + corrections_size + signs_size) / double(data_size) / double(sizeof(K));
            std::cerr << "Epsilon:\t" << epsilon << std::endl;
            std::cerr << "Integer Count:\t" << data_size << std::endl;
            std::cerr << "Segments count:\t" << segments_count << std::endl;
            std::cerr << "Error Point Count:\t" << errorpoint_size << std::endl;
            std::cerr << "Decode Unequal:\t" << data_unequal << std::endl;
            std::cerr << "Average Covered:\t" << avg_covered / segments_count << std::endl;
            std::cerr << "Average Length:\t" << segments_count / index_sequences.size() << std::endl;
            std::cerr << "Segment Size:\t" << segments_size << "\tbyte" << std::endl;
            std::cerr << "Corrections Size:\t" << corrections_size << "\tbyte" << std::endl;
            std::cerr << "Signs Size:\t" << signs_size << "\tbyte" << std::endl;
            std::cerr << "Compression Ratio:\t" << ratio << std::endl;
            long double total_size_in_bytes = segments_size + corrections_size + signs_size;
            long double total_size_in_gib = total_size_in_bytes / 1024.0 / 1024.0 / 1024.0;
            std::cerr << "Total Size:\t" << total_size_in_bytes << "\tin bytes,\t" << total_size_in_gib << "\tin GiB,\t" << total_size_in_bytes / data_size << "\tbytes per int\t" << std::endl;

            if (!output_filename.empty()){
                std::ofstream file(output_filename + ".statistic_log.txt");
                file << "Epsilon:\t" << epsilon << std::endl;
                file << "Integer Count:\t" << data_size << std::endl;
                file << "Segments count:\t" << segments_count << std::endl;
                file << "Error Point Count:\t" << errorpoint_size << std::endl;
                file << "Decode Unequal:\t" << data_unequal << std::endl;
                file << "Average Covered:\t" << avg_covered / segments_count << std::endl;
                file << "Average Length:\t" << segments_count / index_sequences.size() << std::endl;
                file << "Segment Size:\t" << segments_size << "\tbyte" << std::endl;
                file << "Corrections Size:\t" << corrections_size << "\tbyte" << std::endl;
                file << "Signs Size:\t" << signs_size << "\tbyte" << std::endl;
                file << "Compression Ratio:\t" << ratio << std::endl;
                file << "Total Size:\t" << total_size_in_bytes << "\tin bytes,\t" << total_size_in_gib << "\tin GiB,\t" << total_size_in_bytes / data_size << "\tbytes per int\t" << std::endl;
            }
        }

        void huffman_encode() {
            std::cerr << "Huffman Encode Epsilon: " << epsilon << std::endl;
            uint64_t total_huffman_size = 0;
            uint64_t index_count = 0;
            HuffmanEncoder huffman_encoder;
            for (auto& index : index_sequences) {
                index.normal_init();
                huffman_encoder = HuffmanEncoder(index.corrections_vector);
                total_huffman_size += huffman_encoder.calculateSpaceUsage(index.corrections_vector);
                // index_count++;
                // std::cerr << "Index " << index_count << " Total size " << total_huffman_size << std::endl;
            }
            std::cerr << "Total Huffman Size: " << total_huffman_size << " bytes" << std::endl << std::endl;
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
            // random select 8 index
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
            std::cerr << std::endl << "Data Test Epsilon: " << epsilon << std::endl;
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
                std::vector<K> result1 = index.normal_decode();
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
            std::cerr << "Unequal postings: " << data_unequal << std::endl << std::endl;
        }

        void save_model(const std::string output_filename) {
            std::cerr << "Save index to: " << output_filename << std::endl;
            std::ofstream out(output_filename, std::ios::binary);
            out.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size)); // 数据大小
            for (auto& index : index_sequences) {
                out.write(reinterpret_cast<const char*>(&index.n), sizeof(index.n)); // size_t n;
                out.write(reinterpret_cast<const char*>(&index.errorPointCount), sizeof(size_t)); // 错误点数
                out.write(reinterpret_cast<const char*>(&index.first_pos), sizeof(index.first_pos)); // K first_pos;
                //  levels_offsets
                size_t levels_offsets_size = index.levels_offsets.size();
                out.write(reinterpret_cast<const char*>(&levels_offsets_size), sizeof(size_t));
                out.write(reinterpret_cast<const char*>(index.levels_offsets.data()), index.levels_offsets.size() * sizeof(size_t));
                //  segments
                size_t segments_size = index.segments.size();
                out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                out.write(reinterpret_cast<const char*>(index.segments.data()), index.segments.size() * sizeof(typename pgm::PGMIndex<K, epsilon, 0, Floating>::Segment));
                //  bit_vector, signs
                size_t signs_size = index.signs.size();
                out.write(reinterpret_cast<const char*>(&signs_size), sizeof(size_t)); // 大小
                out.write(reinterpret_cast<const char*>(index.signs.data()), (index.signs.size() + 7) / 8); // 数据
                //  rrr_vector, signs_rrr
                // serialize(index.signs_rrr, out);
                //  int_vector, corrections
                size_t corrections_size = index.corrections.size();
                out.write(reinterpret_cast<const char*>(&corrections_size), sizeof(size_t)); // 大小
                out.write(reinterpret_cast<const char*>(index.corrections.data()), index.corrections.size() * sizeof(uint64_t)); // 数据
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
                //  n
                in.read(reinterpret_cast<char*>(&index.n), sizeof(index.n));
                //  errorPointCount
                in.read(reinterpret_cast<char*>(&index.errorPointCount), sizeof(size_t)); // 错误点数
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
                //  signs_rrr
                // index.signs_rrr.load(in);
                //  corrections
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
