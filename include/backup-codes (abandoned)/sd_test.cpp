#include <iostream>
#include <fstream>
#include <chrono>
#include "../../external/mm_file/include/mm_file/mm_file.hpp"
#include "sdsl/sd_vector.hpp"

using namespace std;
using namespace sdsl;
typedef uint32_t K;
size_t data_size;

std::vector<std::vector<K>> data_sequences;
std::vector<bit_vector> bvs;
std::vector<sd_vector<>> sd_vectors;

int main() {

    string input_basename = "/media/sirius/BA167427F361AA8B5/Data/CIFF/Output/rb04-4096/rb04-4096.docs";

    // vector<K> data;
    std::cerr << "Read File: " << input_basename << std::endl;
    mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
    K const* data = input.data();
    assert(data[0] == 1);
    std::cerr << "Universe Size: " << data[1] << std::endl;
    for (size_t i = 2; i < input.size();){
        size_t n = data[i];
        std::vector<K> sequence(data + i + 1, data + i + n);
        data_sequences.emplace_back(sequence);
        data_size += n;
        i += n + 1;
    }
    input.close();
    std::cerr << "Data Size: " << data_size << std::endl;

    bvs.resize(data_sequences.size());

    double count_size_in_bytes = 0;
    for (K i = 0; i < data_sequences.size(); ++i) {
        auto& bv = bvs[i];
        bv.resize(data_sequences[i].size() * 32);
        auto& data_sequence = data_sequences[i];
        for (K j = 0; j < data_sequence.size(); ++j) {
            K tmp = data_sequence[j];
            for (K k = 0; k < 32; ++k) {
                bv[j * 32 + k] = (tmp >> (31 - k)) & 1;
            }
        }
        count_size_in_bytes += size_in_bytes(bv);
    }

    cerr << "Total size: " << count_size_in_bytes << " in Bytes, "  << count_size_in_bytes / 1024.0 / 1024.0 / 1024.0 << "in GiB"<< endl;

    cerr << "Compress" << endl;
    sd_vectors.resize(data_sequences.size());
        count_size_in_bytes = 0;

    for (K i = 0; i < bvs.size(); ++i) {
        auto& bv = bvs[i];
        sd_vectors[i] = sd_vector<>(bv);
        count_size_in_bytes += size_in_bytes(sd_vectors[i]);
    }

    // 计算并输出压缩后的大小
    cerr << "Compress Rate: " << count_size_in_bytes / double(data_size) / double(sizeof(K))  << endl;
    cerr << "Total size: " << count_size_in_bytes << " in Bytes, "  << count_size_in_bytes / 1024.0 / 1024.0 / 1024.0 << "in GiB"<< endl;

    size_t error_point = 0;
    // 测试解压时间
    double avg_time = 0;
    size_t max_time = 0;
    size_t min_time = 0x7fffffff;
    size_t durations = 0;
    int ex_round = 5;
    for(K ex = 0; ex < ex_round + 1; ++ex) {
        auto start = chrono::high_resolution_clock::now();
        for (K i = 0; i < data_sequences.size(); ++i) {
            for (K j = 0; j < data_sequences[i].size(); ++j) {
                auto& sv = sd_vectors[i];
                uint32_t value = 0;

                // uint64_t rank_value = sv.get_int(); // rank操作
                for (K k = 0; k < 32; ++k) {
                    value |= (sv[j * 32 + k] << (31 - k));
                }
                // if (value != data_sequences[i][j] && ex == 0) {
                    // cerr << "Error: " << i << " " << j << " " << value << " " << data_sequences[i][j] << endl;
                    // error_point++;
                // }
            }
        }
        auto end = chrono::high_resolution_clock::now();
        durations = chrono::duration_cast<chrono::microseconds>(end - start).count();
        if (ex != 0) {
            avg_time += durations;
            if (durations > max_time) {
                max_time = durations;
            }
            if (durations < min_time) {
                min_time = durations;
            }
            cerr << "Round: " << ex << ", Time: " << durations << " microseconds" << endl;
        }
    }
    avg_time -= max_time;
    avg_time -= min_time;
    // 输出解压时间
    cerr << "Decode Time: " << avg_time / (ex_round - 2) << " microseconds " << avg_time / (ex_round - 2) * 1000.0 / data_size << endl;
    cerr << "Error Pointer: " << error_point << endl;

    return 0;
}