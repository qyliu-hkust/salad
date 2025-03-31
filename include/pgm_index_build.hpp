#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <pgm_index.hpp>
#include <huffman_encode.hpp>
#include <optPFD_encode.hpp>
#include <../external/mm_file/include/mm_file/mm_file.hpp>

namespace pgm_sequence {

    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))
    template <typename K, size_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t
    class pgm_builder{
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
        std::vector<PGMIndexVariant> index_sequences;
        uint64_t data_size = 0;
        uint64_t data_unequal = 0;
        uint64_t optPFD_size = 0;

        void build_model(std::string input_basename) {
            std::map<uint32_t, uint32_t> epsilon_stats;
            std::cerr << std::endl << "Epsilon: " << epsilon << std::endl;
            std::cerr << "Read File [Build]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            std::cerr << "Universe Size: " << data[1] << std::endl;
            for (size_t i = 2; i < input.size();){
                uint64_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n);

                // OptPFDCompressor compressor;
                // compressor.compress_uint(sequence);
                // optPFD_size += compressor.calculateSpaceUsage();

                if (epsilon == 1) {                     // adapt each list build
                    std::vector<PGMIndexVariant> build_index_sequences;
                    build_index_sequences.push_back(pgm::PGMIndex<K, 1, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 3, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 7, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 15, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 31, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 63, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 127, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 255, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 511, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 1023, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 2047, 0, Floating>(sequence));
                    build_index_sequences.push_back(pgm::PGMIndex<K, 4095, 0, Floating>(sequence));
                    // cerr << "Debug 2" << endl;
                    uint64_t min_size = UINT64_MAX;
                    uint32_t min_index = 0;
                    uint32_t index_count = 0;
                    for (auto &variant_index : build_index_sequences) {
                        std::visit([&min_size, &min_index, &index_count](auto &index) {
                            if (index.total_size_in_bytes() < min_size) {
                                min_size = index.total_size_in_bytes();
                                min_index = index_count;
                            }
                        }, variant_index);
                        index_count++;
                    }
                    index_sequences.push_back(build_index_sequences[min_index]);
                    epsilon_stats[min_index]++;
                }
                else {
                    index_sequences.push_back(pgm::PGMIndex<K, epsilon, 0, Floating>(sequence));
                }
                data_size += n;
                i += n + 1;
            }
            input.close();
            // std::cerr << "OptPFD Size: " << optPFD_size * 4 << " Data Size: " << data_size * 4 << endl;
            for (auto it = epsilon_stats.begin(); it != epsilon_stats.end(); it++) {
                std::cerr << "Epsilon:\t" << (pow(2, it->first + 1)) - 1 << "\tCount:\t" << it->second << std::endl;
            }
            // normal_init();
            std::cerr << "Finish building, then testing" << std::endl;
            data_test(input_basename); // test data
        }

        void normal_init_all() {
            for (auto& variant_index : index_sequences) {
                std::visit([](auto& index) {
                    index.normal_init();
                }, variant_index);
            }
        }

        void statistic_index(std::string output_basename="") {
            normal_init_all();
            double segments_count = 0;
            uint64_t segments_size = 0;
            uint64_t corrections_size = 0;
            uint64_t signs_size = 0;
            uint64_t errorpoint_size = 0;
            uint64_t calculate_by_total_size = 0;
            uint64_t slope_significand_max = 0;
            uint32_t slope_exponent_max = 0;
            double avg_covered = 0;
            for (const auto& variant_index : index_sequences) {
                std::visit([&avg_covered, &segments_count, &segments_size, &corrections_size, &signs_size, &errorpoint_size, &calculate_by_total_size, &slope_exponent_max, &slope_significand_max](auto &index) {
                    for (auto covered : index.seg_covered)
                        avg_covered += covered;
                    slope_exponent_max = std::max(slope_exponent_max, index.segment_slope_exponent_max());
                    slope_significand_max = std::max(slope_significand_max, index.segment_slope_significand_max());
                    calculate_by_total_size += index.total_size_in_bytes();
                    segments_count += index.segments_size;
                    segments_size += index.segment_size_in_bytes();
                    corrections_size += index.corrections_size_in_bytes();
                    signs_size += index.signs_size_in_bytes();
                    errorpoint_size += index.errorPointCount;
                }, variant_index);
            }
            double ratio = (segments_size + corrections_size + signs_size) / double(data_size) / double(sizeof(K));
            std::cerr << "Epsilon:\t" << epsilon << std::endl;
            std::cerr << "Integer Count:\t" << data_size << std::endl;
            std::cerr << "Segments count:\t" << segments_count << std::endl;
            std::cerr << "Error Point Count:\t" << errorpoint_size << std::endl;
            std::cerr << "Average Covered:\t" << avg_covered / segments_count << std::endl;
            std::cerr << "Average Length:\t" << segments_count / index_sequences.size() << std::endl;
            std::cerr << "Segment Size:\t" << segments_size << "\tbyte" << std::endl;
            std::cerr << "Corrections Size:\t" << corrections_size << "\tbyte" << std::endl;
            std::cerr << "Signs Size:\t" << signs_size << "\tbyte" << std::endl;
            std::cerr << "Compression Ratio:\t" << ratio << std::endl;
            long double total_size_in_bytes = segments_size + corrections_size + signs_size;
            long double total_size_in_gib = total_size_in_bytes / 1024.0 / 1024.0 / 1024.0;
            std::cerr << "Total Size:\t" << total_size_in_bytes << "\tin bytes,\t" << total_size_in_gib << "\tin GiB,\t" << total_size_in_bytes / data_size << "\tbytes per int\t" << std::endl;
            std::cerr << "Slope Significand Max:\t" << slope_significand_max << " (" << BIT_WIDTH(slope_significand_max) << ")\tSlope Exponent Max:\t" << slope_exponent_max << " (" << BIT_WIDTH(slope_exponent_max) << ")" << std::endl;
            if (!output_basename.empty()){
                std::ofstream file(output_basename + ".statistic_log.txt");
                file << "Epsilon:\t" << epsilon << std::endl;
                file << "Integer Count:\t" << data_size << std::endl;
                file << "Segments count:\t" << segments_count << std::endl;
                file << "Error Point Count:\t" << errorpoint_size << std::endl;
                file << "Average Covered:\t" << avg_covered / segments_count << std::endl;
                file << "Average Length:\t" << segments_count / index_sequences.size() << std::endl;
                file << "Segment Size:\t" << segments_size << "\tbyte" << std::endl;
                file << "Corrections Size:\t" << corrections_size << "\tbyte" << std::endl;
                file << "Signs Size:\t" << signs_size << "\tbyte" << std::endl;
                file << "Compression Ratio:\t" << ratio << std::endl;
                file << "Total Size:\t" << total_size_in_bytes << "\tin bytes,\t" << total_size_in_gib << "\tin GiB,\t" << total_size_in_bytes / data_size << "\tbytes per int\t" << std::endl;
                file << "Slope Significand Max:\t" << slope_significand_max << " (" << BIT_WIDTH(slope_significand_max) << ")\tSlope Exponent Max:\t" << slope_exponent_max << " (" << BIT_WIDTH(slope_exponent_max) << ")" << std::endl;
            }
        }

        void statistic_gap_list(std::string input_basename, std::string output_basename) {
            std::ofstream file_gap(output_basename + ".gaps.txt");
            std::ofstream file_gap_sta(output_basename + ".gaps_statistic.txt");
            std::cerr << "Read File [Gap]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::vector<double> gap_mean;
            std::vector<double> gap_variance;
            K idx = -1;
            bool save_gap = false;
            // random select 8 lists
            std::vector<K> random_index;
            for (K i = 0; i < 8; i++)
                random_index.push_back(rand() % index_sequences.size());
            for (size_t i = 2; i < input.size();) {
                K n = data[i];
                idx++;
                if (std::find(random_index.begin(), random_index.end(), idx) != random_index.end())
                    save_gap = true;
                else
                    save_gap = false;

                std::vector<K> sequence(data + i + 1, data + i + n);
                if (save_gap)
                    file_gap << sequence.size() - 1 << std::endl;
                // statistic gap mean and variance
                K last = sequence[0];
                std::vector<int> gaps;
                uint64_t sum = 0;
                bool flag = true;
                for (auto value : sequence) {
                    if (flag) {
                        flag = false;
                        continue;
                    }
                    int gap = value - last;
                    if (save_gap)
                        file_gap << gap << std::endl;
                    sum += gap;
                    gaps.push_back(gap);
                    last = value;
                }
                // std::cerr << n << " " <<  sequence.size() << " " << gaps.size() <<  std::endl;
                double mean = sum / double(gaps.size());
                gap_mean.push_back(mean);
                double variance = 0;
                for (auto& gap : gaps) {
                    variance += (gap - mean) * (gap - mean);
                }
                variance /= gaps.size();
                gap_variance.push_back(variance);
                file_gap_sta << "Mean:\t" << mean << "\tVariance:\t" << variance << std::endl;
                i += n + 1;
            }
            input.close();
            double gap_mean_sum = 0;
            double gap_variance_sum = 0;
            for (auto& mean : gap_mean)
                gap_mean_sum += mean;
            for (auto& variance : gap_variance)
                gap_variance_sum += variance;
            file_gap_sta << "Total Mean:\t" << gap_mean_sum / gap_mean.size() << "\tTotal Variance:\t" << gap_variance_sum / gap_variance.size() << std::endl;
            std::cerr << "Total Mean:\t" << gap_mean_sum / gap_mean.size() << "\tTotal Variance:\t" << gap_variance_sum / gap_variance.size() << std::endl;
            std::cerr << "Save Gaps List to: " << output_basename << ".gaps.txt" << std::endl;
            std::cerr << "Save Gaps Statistic to: " << output_basename << ".gaps_statistic.txt" << std::endl;
            file_gap.close();
            file_gap_sta.close();
        }

        void statistic_gap_segment(std::string output_basename) {
            std::ofstream file_gap(output_basename + ".seg_gaps.txt");
            std::ofstream file_gap_sta(output_basename + ".seg_gaps_statistic.txt");
            std::vector<K> random_index;
            for (K i = 0; i < 8; i++)
                random_index.push_back(rand() % index_sequences.size());
            for (auto i : random_index) {
                auto& variant_index = index_sequences[i];
                std::visit([&file_gap, &file_gap_sta, &i](auto &index) {
                    index.normal_init();
                    // file_gap << i << std::endl;
                    file_gap << index.segments.size() - 1 << std::endl;
                    for (K idx = 0; idx < index.segments.size() - 1; idx++) {
                        std::vector<K> segment_list = index.segment_decode(idx);
                        K n = segment_list.size();
                        file_gap << (n - 1) << std::endl;
                        K last = segment_list[0];
                        std::vector<int> gaps;
                        uint64_t sum = 0;
                        bool flag = true;
                        for (auto value : segment_list) {
                            if (flag) {
                                flag = false;
                                continue;
                            }
                            int gap = value - last;
                            file_gap << gap << std::endl;
                            sum += gap;
                            gaps.push_back(gap);
                            last = value;
                        }
                        // std::cerr <<  segment_list.size() << " " << gaps.size() <<  std::endl;
                        double mean = sum / double(gaps.size());
                        double variance = 0;
                        for (auto& gap : gaps) {
                            variance += (gap - mean) * (gap - mean);
                        }
                        variance /= gaps.size();
                        file_gap_sta << "Index:\t" << i << "\tEpsilon:\t" << index.Epsilon_Data << "\tSegment:\t" << idx << "\tMean:\t" << mean << "\tVariance:\t" << variance << std::endl;
                    }
                }, variant_index);
            }
            std::cerr << "Save Segment Gaps List to: " << output_basename << ".seg_gaps.txt" << std::endl;
            std::cerr << "Save Segment Gaps Statistic to: " << output_basename << ".seg_gaps_statistic.txt" << std::endl;
            file_gap.close();
            file_gap_sta.close();
        }

        void huffman_encode() {
            std::cerr << "Huffman Encode Epsilon: " << epsilon << std::endl;
            uint64_t total_huffman_size = 0;

            for (auto& variant_index : index_sequences) {
                std::visit([&total_huffman_size](auto &index) {
                    index.huffman_init();
                    HuffmanEncoder huffman_encoder = HuffmanEncoder(index.corrections_vector, BIT_WIDTH(index.Epsilon_Data));;
                    total_huffman_size += huffman_encoder.calculateSpaceUsage(index.corrections_vector);
                }, variant_index);
            }
            total_huffman_size = total_huffman_size / 8;
            std::cerr << "Total Huffman Size: " << total_huffman_size << " bytes" << std::endl << std::endl;
        }

        void optPFD_encode(){
            // std::cerr << "OptPFD Encoding: " << std::endl;
            uint64_t total_size = 0;
            for (auto& variant_index : index_sequences) {
                std::visit([&total_size](auto &index) {
                    index.normal_init();
                    OptPFDCompressor encoder;
                    encoder.compress_int(index.corrections_vector);
                    total_size += encoder.calculateSpaceUsage();
                }, variant_index);
            }
            total_size = total_size / 8;
            std::cerr << "Total OptPFD Size: " << total_size << " bytes" << std::endl << std::endl;
        }

        void save_covered(const std::string output_basename) {
            std::ofstream file(output_basename + ".covered.txt");
            if (file.is_open()) {
                for (auto& variant_index : index_sequences) {
                    std::visit([&file](auto &index) {
                        file << index.segments.size() << std::endl;
                        for (auto& segment : index.segments) {
                            file << segment.covered << std::endl;
                        }
                    }, variant_index);
                }
                file.close();
            }
            else {
                std::cerr << "Couldn't open the log_file" << std::endl;
            }
        }

        void save_residual(const std::string output_basename) {
            std::ofstream file(output_basename + ".residual.txt");
            // random select 8 index
            std::vector<K> random_index;
            for (K i = 0; i < 8; i++) {
                random_index.push_back(rand() % index_sequences.size());
            }
            if (file.is_open()) {
                for (auto i : random_index) {
                    auto& variant_index = index_sequences[i];
                    std::visit([&file](auto &index) {
                        index.normal_init();
                        file << index.n << std::endl;
                        for (auto& correction : index.corrections_vector) {
                            file << correction << std::endl;
                        }
                    }, variant_index);
                }
                file.close();
            }
            else {
                std::cerr << "Couldn't open the log_file" << std::endl;
            }
        }

        void save_residual_random_segment(const std::string output_basename) {
            std::ofstream file(output_basename + ".random_segment_residual.txt");
            std::vector<K> random_index;
            for (K i = 0; i < 16; i++) {
                random_index.push_back(rand() % index_sequences.size());
            }
            if (file.is_open()) {
                for (auto i : random_index) {
                    auto& variant_index = index_sequences[i];
                    std::visit([&file](auto &index) {
                        int save_seg_idx = rand() % (index.segments.size() + 5);
                        index.reoprt_residual_random_segment(file, save_seg_idx);
                    }, variant_index);
                }
                file.close();
                std::cerr << "Save Random Segment Residual to: " << output_basename << ".random_segment_residual.txt" << std::endl;
            }
            else {
                std::cerr << "Couldn't open the log_file: " << output_basename << std::endl;
            }
        }

        void data_test(std::string input_basename) {
            normal_init_all();
            std::cerr << std::endl << "Data Test Epsilon: " << epsilon << std::endl;
            std::cerr << "Read File [Test]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::cerr << "Universe Size: " << data[1] << std::endl;
            K data_unequal_test = 0;
            K posi = 0;
            for (size_t i = 2; i < input.size();) {
                K n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n);
                auto& variant_index = index_sequences[posi++];
                std::visit([&sequence, &data_unequal_test](auto &index) {
                    std::vector<K> result1 = index.normal_decode();
                    if (result1.size() != sequence.size())
                        std::cerr << "ERROR: decode size error " << result1.size() << " " << sequence.size() << std::endl;
                    for (auto j = 0; j < result1.size(); j++) {
                        if (sequence[j] != result1[j]) {
                            data_unequal_test++;
                        }
                    }
                }, variant_index);

                i += n + 1;
            }
            std::cerr << "Unequal postings: " << data_unequal_test << std::endl << std::endl;
        }

        // void init_first_pos_index(const std::string output_basename) {
        //     if (output_basename.back() != '/') {
        //         std::cerr << "Error: output_basename must end with '/'" << std::endl;
        //         return;
        //     }
        //     K index_count = 0;
        //     for (auto& variant_index : index_sequences) {
        //         string out_string = output_basename + std::to_string(index_count) + ".sp_tree";
        //         std::visit([&out_string](auto &index) {
        //             index.build_first_pos_index(out_string);
        //         }, variant_index);
        //         index_count++;
        //     }
        //     std::cerr << "Finish Build First Pos Index using S+ Tree" << std::endl;
        // }

        void save_model(const std::string output_basename) {
            //判断是否以/结尾
            if (output_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }
            std::cerr << "Save index to: " << output_basename << std::endl;
            std::ofstream out(output_basename + "idx.size", std::ios::binary);
            K index_num = index_sequences.size();
            out.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size)); // 数据大小
            out.write(reinterpret_cast<const char*>(&index_num), sizeof(index_num)); // 索引数量
            out.close();
            K index_count = 0;
            for (auto & variant_index: index_sequences) {
                std::ofstream out(output_basename + std::to_string(index_count) + ".idx", std::ios::binary);
                std::visit([&out](auto &index) {
                    // epsilon
                    out.write(reinterpret_cast<const char*>(&index.Epsilon_Data), sizeof(index.Epsilon_Data)); // epsilon
                    // n
                    out.write(reinterpret_cast<const char*>(&index.n), sizeof(index.n)); // size_t n;
                    // first_pos
                    out.write(reinterpret_cast<const char*>(&index.first_pos), sizeof(index.first_pos)); // K first_pos;
                    //  levels_offsets
                    size_t levels_offsets_size = index.levels_offsets.size();
                    out.write(reinterpret_cast<const char*>(&levels_offsets_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(index.levels_offsets.data()), index.levels_offsets.size() * sizeof(size_t));

                    //  segments
                    // size_t segments_size = index.segments.size();
                    // out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    // out.write(reinterpret_cast<const char*>(index.segments.data()), index.segments.size() * sizeof(typename pgm::PGMIndex<K, epsilon, 0, Floating>::Segment));
                    out.write(reinterpret_cast<const char*>(&index.segments_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(&index.bpc_first), sizeof(uint8_t));
                    out.write(reinterpret_cast<const char*>(&index.bpc_covered), sizeof(uint8_t));
                    out.write(reinterpret_cast<const char*>(&index.bpc_intercept), sizeof(uint8_t));
                    out.write(reinterpret_cast<const char*>(&index.bpc_slope_exponent), sizeof(uint8_t));
                    out.write(reinterpret_cast<const char*>(&index.bpc_slope_significand), sizeof(uint8_t));

                    size_t segments_size = index.seg_first_compress.size();
                    out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(index.seg_first_compress.data()), index.seg_first_compress.size() * sizeof(uint64_t));
                    segments_size = index.seg_covered_compress.size();
                    out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(index.seg_covered_compress.data()), index.seg_covered_compress.size() * sizeof(uint64_t));
                    segments_size = index.seg_intercept_compress.size();
                    out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(index.seg_intercept_compress.data()), index.seg_intercept_compress.size() * sizeof(uint64_t));
                    segments_size = index.seg_slope_exponent_compress.size();
                    out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(index.seg_slope_exponent_compress.data()), index.seg_slope_exponent_compress.size() * sizeof(uint64_t));
                    segments_size = index.seg_slope_significand_compress.size();
                    out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    out.write(reinterpret_cast<const char*>(index.seg_slope_significand_compress.data()), index.seg_slope_significand_compress.size() * sizeof(uint64_t));


                    // out.write(reinterpret_cast<const char*>(&segments_size), sizeof(size_t));
                    // out.write(reinterpret_cast<const char*>(index.seg_first.data()), index.seg_first.size() * sizeof(uint32_t));
                    // out.write(reinterpret_cast<const char*>(index.seg_covered.data()), index.seg_covered.size() * sizeof(uint32_t));
                    // out.write(reinterpret_cast<const char*>(index.seg_intercept.data()), index.seg_intercept.size() * sizeof(int32_t));
                    // out.write(reinterpret_cast<const char*>(index.seg_slope_exponent.data()), index.seg_slope_exponent.size() * sizeof(uint8_t));
                    // out.write(reinterpret_cast<const char*>(index.seg_slope_significand.data()), index.seg_slope_significand.size() * sizeof(int64_t));

                    //  bit_vector, signs_compress, signs_exception
                    size_t signs_size = index.signs_compress.size();
                    out.write(reinterpret_cast<const char*>(&signs_size), sizeof(size_t)); // 大小
                    out.write(reinterpret_cast<const char*>(index.signs_compress.data()), (index.signs_compress.size() + 7) / 8); // 数据
                    signs_size = index.signs_exception.size();
                    out.write(reinterpret_cast<const char*>(&signs_size), sizeof(size_t)); // 大小
                    out.write(reinterpret_cast<const char*>(index.signs_exception.data()), (index.signs_exception.size() + 7) / 8); // 数据

                    //  int_vector, corrections_compress, corrections_exception
                    size_t corrections_size = index.corrections_compress.size();
                    out.write(reinterpret_cast<const char*>(&corrections_size), sizeof(size_t)); // 大小
                    out.write(reinterpret_cast<const char*>(index.corrections_compress.data()), index.corrections_compress.size() * sizeof(uint64_t)); // 数据
                    corrections_size = index.corrections_exception.size();
                   out.write(reinterpret_cast<const char*>(&corrections_size), sizeof(size_t)); // 大小
                   out.write(reinterpret_cast<const char*>(index.corrections_exception.data()), index.corrections_exception.size() * sizeof(uint64_t)); // 数据

                    // index.build_first_pos_index(out_string);
                }, variant_index);
                out.close();
                index_count++;
            }
            out.close();
        }

        void load_model(std::string input_basename) {
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }
            std::cerr << "Load index from: " << input_basename << std::endl;
            index_sequences.clear(); // 清空现有数据
            K index_num = 0;
            std::ifstream in(input_basename + "idx.size", std::ios::binary);
            in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            in.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));
            in.close();
            for (K index_count = 0; index_count < index_num; index_count++) {
                std::ifstream in(input_basename + std::to_string(index_count) + ".idx", std::ios::binary);
                uint64_t Epsilon_Data = 0;
                while (in.peek() != EOF) {
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

                        // cerr << "Debug 1" << endl;
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
                    index_sequences.push_back(std::move(variant_index));
                }
                in.close();
            }
            // cerr << "Debug 6" << endl;
        }
    };
}
