#pragma once

// #include <pgm_index_enumerate.hpp>

#include "pgm_index_enumerate.hpp"

namespace pgm_sequence {
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll((x)))

    template <typename K=uint32_t, uint64_t epsilon = 1, typename Floating=double> // K is uint32_t or uint64_t, Floating is unused
    class pgm_index_residual_compressor {
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
        uint64_t total_residual_bit_size_distance = 0;
        uint64_t total_residual_bit_size_flag = 0;
        uint64_t original_residual_bit_size = 0;
        int compression_bit = 2;
        void compress_test(pgm_enumerator<K> &index, const std::string log_filename) {
            index.second_difference_compress();
            total_residual_bit_size_distance += index.total_residual_bit_size_distance;
            total_residual_bit_size_flag += index.total_residual_bit_size_flag;
            original_residual_bit_size += index.corrections_vector.size() * BIT_WIDTH(index.Epsilon_Data);
            compression_bit = index.compress_bit;
            std::ofstream out(log_filename + ".compress-residual-log.txt", std::ios_base::app);
            out << "Epsilon:\t" << index.Epsilon_Data << "\t(" << BIT_WIDTH(index.Epsilon_Data)  << ")";
            out << "\tResidual_Min:\t" << index.residual_min << "\t(" << BIT_WIDTH(abs(index.residual_min)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_min));
            out << "\tResidual_5%:\t" << index.residual_5 << "\t(" << BIT_WIDTH(abs(index.residual_5)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_5));
            out << "\tResidual_10%:\t" << index.residual_10 << "\t(" << BIT_WIDTH(abs(index.residual_10)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_10));
            out << "\tResidual_20%:\t" << index.residual_20 << "\t(" << BIT_WIDTH(abs(index.residual_20)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_20));
            out << "\tResidual_25%:\t" << index.residual_25 << "\t(" << BIT_WIDTH(abs(index.residual_25)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_25));
            out << "\tResidual_50%:\t" << index.residual_50 << "\t(" << BIT_WIDTH(abs(index.residual_50)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_50));
            out << "\tResidual_75%:\t" << index.residual_75 << "\t(" << BIT_WIDTH(abs(index.residual_75)) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(abs(index.residual_75));
            out << "\tResidual_80%:\t" << index.residual_80 << "\t(" << BIT_WIDTH(index.residual_80) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(index.residual_80);
            out << "\tResidual_90%:\t" << index.residual_90 << "\t(" << BIT_WIDTH(index.residual_90) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(index.residual_90);
            out << "\tResidual_95%:\t" << index.residual_95 << "\t(" << BIT_WIDTH(index.residual_95) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(index.residual_95);
            out << "\tResidual_Max:\t" << index.residual_max << "\t(" << BIT_WIDTH(index.residual_max) << ")\t" << BIT_WIDTH(index.Epsilon_Data) - BIT_WIDTH(index.residual_max);
            out << "\tOver_Rate:\t" << index.over_num << "\t" << index.n << "\t" << index.over_data;
            out << "\tMax_Distance:\t" << index.max_distance << "\t(" << BIT_WIDTH(index.max_distance) << ")";
            out << endl;
        }

        void test_model(std::string input_basename, const std::string log_filename) {
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }
            std::cerr << "Load index from: " << input_basename << std::endl;
            uint64_t data_size = 0;
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

                        pgm_enumerator<K> residual_compressor_tmp;
                        std::visit([&residual_compressor_tmp](auto &index) {
                            residual_compressor_tmp.Epsilon_Data = index.Epsilon_Data;
                            index.normal_init();
                            for (auto seg : index.segments) {
                                auto first = seg.first;
                                auto intercept = seg.intercept;
                                auto slope_significand = seg.slope_significand;
                                auto slope_exponent = seg.slope_exponent;
                                auto covered = seg.covered;
                                residual_compressor_tmp.segments.emplace_back(first, intercept, slope_exponent, slope_significand, covered);
                            }
                            residual_compressor_tmp.load_copy(index.n, index.corrections_vector);
                        }, variant_index);
                        compress_test(residual_compressor_tmp, log_filename);
                    }
                }
                in.close();
            }
            statistic();
        }

        void statistic() {
            cout << "Compression Bit:\t" << compression_bit;
            cout << "\tOriginal Bit Size:\t" << original_residual_bit_size << endl;
            cout << "\tTotal Bit Size Distance:\t" << total_residual_bit_size_distance;
            cout << "\tCompression Ratio:\t" << (double) total_residual_bit_size_distance / (double) original_residual_bit_size;
            cout << "\tTotal Bit Size Flag:\t" << total_residual_bit_size_flag;
            cout << "\tCompression Ratio:\t" << (double) total_residual_bit_size_flag / (double) original_residual_bit_size << endl;
        }
    };
}