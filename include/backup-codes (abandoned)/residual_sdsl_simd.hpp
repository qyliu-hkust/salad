#pragma once

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
using namespace sdsl;

namespace rsdsdslsimd{
template<typename K, int bits = 8>
class RSDSDSLSIMD {
protected:

    template <typename Index, typename Data>
    static void build(Index index, Data data, 
                        sdsl::int_vector<>& rsd_data,
                        sdsl::rrr_vector<63>& rsd_symbol) {

        sdsl::bit_vector rsd_symbol_item(data.size());  
        std::vector<K> approxKeyVecStore = index.seqSearchApproxKey(data.size());  
        

        for (size_t pos = 0; pos + 3 < data.size(); pos += 4){

            // Load data into AVX2 vectors
            __m256i true_key_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[pos]));
            __m256i approx_key_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&approxKeyVecStore[pos]));

            // Compute the absolute difference
            __m256i diff_vec = _mm256_sub_epi64(true_key_vec, approx_key_vec);

            // Compute the comparison mask
            __m256i cmp_mask = _mm256_cmpgt_epi64(true_key_vec, approx_key_vec);

            // Store the results, NO have -O3, sdsl::int_vector 默认是无符号, 需转换为有符号且bits+1
            rsd_data[pos + 0] =  (_mm256_extract_epi64(diff_vec, 0) < 0) ? -(_mm256_extract_epi64(diff_vec, 0)) : _mm256_extract_epi64(diff_vec, 0);
            rsd_data[pos + 1] =  (_mm256_extract_epi64(diff_vec, 1) < 0) ? -(_mm256_extract_epi64(diff_vec, 1)) : _mm256_extract_epi64(diff_vec, 1);
            rsd_data[pos + 2] =  (_mm256_extract_epi64(diff_vec, 2) < 0) ? -(_mm256_extract_epi64(diff_vec, 2)) : _mm256_extract_epi64(diff_vec, 2);
            rsd_data[pos + 3] =  (_mm256_extract_epi64(diff_vec, 3) < 0) ? -(_mm256_extract_epi64(diff_vec, 3)) : _mm256_extract_epi64(diff_vec, 3);

            rsd_symbol_item[pos + 0] = (_mm256_extract_epi64(cmp_mask, 0) != 0) ? 0 : 1;
            rsd_symbol_item[pos + 1] = (_mm256_extract_epi64(cmp_mask, 1) != 0) ? 0 : 1;
            rsd_symbol_item[pos + 2] = (_mm256_extract_epi64(cmp_mask, 2) != 0) ? 0 : 1;
            rsd_symbol_item[pos + 3] = (_mm256_extract_epi64(cmp_mask, 3) != 0) ? 0 : 1;
    
        }

        // util::bit_compress(rsd_data);   
        rsd_symbol = sdsl::rrr_vector<63>(rsd_symbol_item);
    }

public:

     sdsl::int_vector<> rsd_data;
     sdsl::rrr_vector<63> rsd_symbol;

    RSDSDSLSIMD() = default;

    template <typename Index, typename Data>
    RSDSDSLSIMD(const Index& index, const Data& data)
         :rsd_data(data.size(), 0 , bits), rsd_symbol(){
        build(index, data, rsd_data, rsd_symbol);
    }

    size_t rsd_data_bytes() const { return size_in_bytes(rsd_data); }
    size_t rsd_symbol_bytes() const { return size_in_bytes(rsd_symbol); }

    // 对rsd_data进行合格性检查
    size_t checkResidual(size_t epsilon) const { 
        size_t size = rsd_data.size();
        size_t count = 0;
        size_t i = 0;
        __m256i epsilon_vec = _mm256_set1_epi64x(epsilon + 1);

        for (; i + 3 < size; i += 4) {
            __m256i data_vec = _mm256_set_epi64x(rsd_data[i+3], rsd_data[i+2], rsd_data[i+1], rsd_data[i]);
            __m256i cmp_mask = _mm256_cmpgt_epi64(data_vec, epsilon_vec);   // dst[i+63:i] := ( a > b ) ? 1 : 0
            
            // 提取掩码，计算大于epsilon的数量
            int mask = _mm256_movemask_pd(_mm256_castsi256_pd(cmp_mask));   // if 1 dst = 1 , else dst = 0, mask 是一个数组
            count += __builtin_popcount(mask);   // 计算掩码中为 1 的位的数量，即大于 epsilon 的元素数量
        }

        // 处理剩余元素
        for (; i < size; ++i) {
            if (rsd_data[i] > epsilon + 1) {
                ++count;
            }
        }

        return count;
    }

};
}