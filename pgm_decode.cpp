#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h> // Linux
#include <sys/types.h>
#include "pgm_index_decode.hpp"
// using namespace std;

std::string log_filename = "";
std::string decode_type = "";
bool read_only = false;

template <uint64_t epsilon=64>
void test_collection_pgm(const std::string input_basename, const std::string output_basename) {

    typedef pgm_sequence::pgm_decoder<uint32_t, epsilon> PGM_INDEX_DECODER;
    PGM_INDEX_DECODER index;

    index.test_model(output_basename, decode_type);
}

int main(int argc, const char** argv)
{
//    cpu_set_t mask;
//    CPU_ZERO(&mask); // clear CPU mask
//    CPU_SET(16, &mask); // set CPU 16
//
//    if (sched_setaffinity(0, sizeof(mask), &mask) == -1) { // 0 is the calling process
//        perror("sched_setaffinity");
//        return 1;
//    }

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
    uint32_t epsilon = static_cast<uint32_t>(std::stoi(epsilonstr));
    read_only = (read_only_str == "t");

    if (index_type == "pgm")
        switch (epsilon)
        {
            case 1: test_collection_pgm<1>(input_basename, output_basename); break;
            case 15: test_collection_pgm<15>(input_basename, output_basename); break;
            // case 16: test_collection_pgm<16>(input_basename, output_basename); break;
            case 31: test_collection_pgm<31>(input_basename, output_basename); break;
            // case 32: test_collection_pgm<32>(input_basename, output_basename); break;
            case 63: test_collection_pgm<63>(input_basename, output_basename); break;
            // case 64: test_collection_pgm<64>(input_basename, output_basename); break;
            // case 126: test_collection_pgm<126>(input_basename, output_basename); break;
            case 127: test_collection_pgm<127>(input_basename, output_basename); break;
            // case 128: test_collection_pgm<128>(input_basename, output_basename); break;
            case 255: test_collection_pgm<255>(input_basename, output_basename); break;
            // case 256: test_collection_pgm<256>(input_basename, output_basename); break;
            case 511: test_collection_pgm<511>(input_basename, output_basename); break;
            // case 512: test_collection_pgm<512>(input_basename, output_basename); break;
            case 1023: test_collection_pgm<1023>(input_basename, output_basename); break;
            // case 1024: test_collection_pgm<1024>(input_basename, output_basename); break;
            case 2047: test_collection_pgm<2047>(input_basename, output_basename); break;
            // case 2048: test_collection_pgm<2048>(input_basename, output_basename); break;
            case 4095: test_collection_pgm<4095>(input_basename, output_basename); break;
            default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
        }
    else
        std::cerr << "ERROR: only support pgm index " << std::endl;

    return 0;
}
