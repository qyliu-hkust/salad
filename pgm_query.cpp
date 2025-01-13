#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h> // Linux
#include <sys/types.h>
#include "pgm_index_query.hpp"
// using namespace std;

std::string log_filename = "";
std::string decode_type = "";
std::string query_filename = "";
std::string query_type = "";
bool read_only = false;

template <uint64_t epsilon=64>
void test_collection_pgm(const std::string input_basename, const std::string output_filename) {

    typedef pgm_sequence::pgm_querier<uint32_t, epsilon> PGM_INDEX_QUERIER;
    PGM_INDEX_QUERIER index;
    index.test_query(output_filename, decode_type, query_filename, query_type);

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

    int mandatory = 8;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <index_filename> <epsilon> <read_only> <decode_type> <query_filename> <query_type> <log_filename>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_basename = argv[3];
    const std::string epsilon_str = argv[4];
    uint32_t epsilon = static_cast<uint32_t>(std::stoi(epsilon_str));
    const std::string read_only_str = argv[5];
    read_only = (read_only_str == "t");
    decode_type = argv[6];
    query_filename = argv[7];
    query_type = argv[8];
    log_filename = argv[9];

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
            // case 4096: test_collection_pgm<4096>(input_basename, output_basename); break;
            // case 8192: test_collection_pgm<8192>(input_basename, output_basename); break;
            // case 16384: test_collection_pgm<16384>(input_basename, output_basename); break;
            // case 32768: test_collection_pgm<32768>(input_basename, output_basename); break;
            // case 65536: test_collection_pgm<65536>(input_basename, output_basename); break;
            // case 131072: test_collection_pgm<131072>(input_basename, output_basename); break;
            // case 262144: test_collection_pgm<262144>(input_basename, output_basename); break;
            default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
        }
    else
        std::cerr << "ERROR: only support pgm index " << std::endl;
    return 0;
}
