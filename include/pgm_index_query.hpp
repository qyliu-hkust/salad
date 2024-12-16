#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include "pgm_index.hpp"

namespace pgm_sequence {
    template <typename K, uint64_t epsilon = 64, typename Floating=double> // K is uint32_t or uint64_t, Floating is unused
    class pgm_querier{
    protected:

        uint64_t data_size = 0;
        uint32_t query_num = 0;
        uint64_t query_time = 0;

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

        std::vector<pgm::PGMIndex<K, epsilon, 0, Floating>> load_model(std::string input_filename) {
            std::cerr << "Load index from: " << input_filename << std::endl;
            std::ifstream in(input_filename, std::ios::binary);
            if (!in) {
                throw std::runtime_error("Failed to open file for reading: " + input_filename);
            }
            std::vector<pgm::PGMIndex<K, epsilon, 0, Floating>> index_sequences;
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
                    index_sequences.push_back(std::move(index));
                }
            }
            in.close();
            return index_sequences;
        }

        void query_and_test(std::vector<pgm::PGMIndex<K, epsilon, 0, Floating>> &index_sequences, const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            std::vector<uint32_t> result;
            for (auto &query : query_list) {
                uint32_t query_id_idx = 0;
                std::vector<uint32_t> candidate_id_idx(query.size(), 0);
                for (auto &query_idx : query) {
                    index_sequences[query_idx].query_init();
                }
                uint32_t candidate_posting = 0;
                uint32_t equal_result = 0;

                while (true) {
                    uint32_t candidate_posting_tmp = index_sequences[query[query_id_idx]].nextgeq(candidate_posting, decode_type);
                    if (candidate_posting_tmp == INT_MAX) {
                        break;
                    } else if (candidate_posting == candidate_posting_tmp) {
                        equal_result++;
                    } else {
                        candidate_posting = candidate_posting_tmp;
                    }
                    query_id_idx = (query_id_idx + 1) % query.size();
                }
            }

        }



    public:
        void test_query(const std::string &input_filename, const std::string &decode_type, const std::string &query_filename, const std::string query_type) {

            std::vector<std::vector<uint32_t>> query_list = read_query(query_filename);
            // for (auto &query : query_list) {
            //     for (auto &id : query) {
            //         std::cout << id << " ";
            //     }
            //     std::cout << std::endl;
            // }
            std::vector<pgm::PGMIndex<K, epsilon, 0, Floating>> index_sequences = load_model(input_filename);
            if (query_type == "AND")
                query_and_test(index_sequences, query_list, decode_type);



        }


    };
}
