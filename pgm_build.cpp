#include <vector>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <unistd.h> // Linux
#include "pgm_index_build.hpp"

std::string log_filename = "";
std::string decode_type = "";
bool read_only = false;

template <uint64_t epsilon=64>
void create_collection_pgm(const string input_basename, const string output_basename) {

    typedef pgm_sequence::pgm_builder<uint32_t, epsilon> PGM_INDEX_BUILDER;
    PGM_INDEX_BUILDER index;
    if (!read_only){
        index.build_model(input_basename + ".docs");
        // index.normal_init();
        // index.data_test(input_basename + ".docs");
        index.save_model(output_basename);
        if(!log_filename.empty()) { // save covered and residual
            index.statistic_index(log_filename);
            // index.statistic_gap_list(input_basename + ".docs", log_filename);
            // index.statistic_gap_segment(log_filename);
            // index.save_residual(log_filename);
        }
    } else {
        index.load_model(output_basename);
        // index.optPFD_encode();
        // index.init_first_pos_index(output_basename); // build s+ tree for first_pos

        if(!log_filename.empty()) { // save covered and residual
            index.statistic_index(log_filename);
            // index.statistic_gap_list(input_basename + ".docs", log_filename);
            // index.statistic_gap_segment(log_filename);
            // index.save_residual(log_filename);

            index.save_residual_random_segment(log_filename);
        }
        // index.statistic_index(log_filename);
        // index.normal_init();
        index.data_test(input_basename + ".docs");
         // index.save_residual(log_filename);
        // index.huffman_encode();
    }
}

int main(int argc, const char** argv)
{
    ios::sync_with_stdio(0);
    int mandatory = 6;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <output_basename> <epsilon> <read_only> <decode_type> <test_log>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_basename = argv[3];
    const std::string epsilonstr = argv[4];
    const std::string read_only_str = argv[5];
    decode_type = argv[6];
    log_filename = argv[7];
    const uint64_t epsilon = static_cast<uint64_t>(std::stoi(epsilonstr));
    read_only = (read_only_str == "t");

    if (index_type == "pgm")  //compress docs
        switch (epsilon)
        {
            case 1: create_collection_pgm<1>(input_basename, output_basename); break;
            // case 15: create_collection_pgm<15>(input_basename, output_basename); break;
            // case 31: create_collection_pgm<31>(input_basename, output_basename); break;
            // case 63: create_collection_pgm<63>(input_basename, output_basename); break;
            // case 127: create_collection_pgm<127>(input_basename, output_basename); break;
            // case 255: create_collection_pgm<255>(input_basename, output_basename); break;
            // case 511: create_collection_pgm<511>(input_basename, output_basename); break;
            // case 1023: create_collection_pgm<1023>(input_basename, output_basename); break;
            // case 2047: create_collection_pgm<2047>(input_basename, output_basename); break;
            // case 4095: create_collection_pgm<4095>(input_basename, output_basename); break;
            default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
        }
    else
        std::cerr << "ERROR: only support pgm index " << std::endl;

    return 0;
}
