# pragma once

#include <vector>
#include <cstdint>
#include <codecfactory.h>
#include <deltautil.h>

class OptPFDCompressor {
public:
    size_t real_data_size;
    std::vector<uint32_t> compressed;
    size_t compressedSize;

    const std::string codec_name = "simdfastpfor256";
    // const std::string codec_name = "simdoptpfor";
    // const std::string codec_name = "simple16";

    // const std::string codec_name = "fastbinarypacking32";
    void compress_int(const std::vector<int> &input_int) {
        FastPForLib::CODECFactory factory;
        FastPForLib::IntegerCODEC &codec = *factory.getFromName(codec_name);

        std::vector<uint32_t> input(input_int.size());
        for (size_t i = 0; i < input_int.size(); i++) {
            input[i] = input_int[i] >= 0 ? input_int[i] : -input_int[i];
        }
        std::sort(input.begin(), input.end());
        real_data_size = input.size();
        compressed.resize(input.size() + 4096);
        size_t compressedSize = compressed.size();
        codec.encodeArray(input.data(), input.size(), compressed.data(), compressedSize);
        compressed.resize(compressedSize);
        compressed.shrink_to_fit();

        // compressionRatio(input, compressed);
    }

    void compress_uint(std::vector<uint32_t> input) {
        FastPForLib::CODECFactory factory;
        FastPForLib::IntegerCODEC &codec = *factory.getFromName(codec_name);
        // sort(input.begin(), input.end());
        real_data_size = input.size();
        compressed.resize(input.size() + 4096);
        compressedSize = compressed.size();
        codec.encodeArray(input.data(), input.size(), compressed.data(), compressedSize);
        compressed.resize(compressedSize);
        compressed.shrink_to_fit();

        // compressionRatio(input, compressed);
    }

    std::vector<uint32_t> decompress() {
        FastPForLib::CODECFactory factory;
        FastPForLib::IntegerCODEC &codec = *factory.getFromName(codec_name);
        std::vector<uint32_t> decompressed(real_data_size);
        size_t recoveredSize = decompressed.size();
        codec.decodeArray(compressed.data(), compressed.size(), decompressed.data(), recoveredSize);
        decompressed.resize(recoveredSize);
        return decompressed;
    }

    void compressionRatio(const std::vector<uint32_t> &input, const std::vector<uint32_t> &compressed) {
        std::cout << std::setprecision(3);
        std::cout << "You are using "
                  << 32.0 * static_cast<double>(compressed.size()) /
                         static_cast<double>(input.size())
                  << " bits per integer. " << std::endl;
        std::cout << "Compression ratio: "
                  << 100.0 * static_cast<double>(input.size()) /
                         static_cast<double>(compressed.size())
                  << " percentage" << std::endl;
    }

    uint64_t calculateSpaceUsage() {
        return compressedSize;
    }
};