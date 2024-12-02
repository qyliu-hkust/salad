#pragma once

#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
using namespace std;
using namespace sdsl;

namespace sdslrsd{
template<typename K, int bits = 8>
class SDSLRESIDUAL {
protected:

    template <typename Index, typename Data>
    static void build(Index index, Data data, 
                        sdsl::int_vector<>& rsd_data,
                        sdsl::rrr_vector<63>& rsd_symbol) {
        // sdsl::int_vector<> rsd_data_item(data.size(), 0, n);
        sdsl::bit_vector rsd_symbol_item(data.size());  
        for (size_t pos = 0; pos < data.size(); pos++){
            K true_key = data[pos];  
            K approximate_key = index.search(pos);  // 既然是依次搜索，每次调用index.search()计算approx_key太耗时。使用simd加速计算。
            rsd_data[pos] = (true_key > approximate_key ? true_key - approximate_key : approximate_key - true_key);
            rsd_symbol_item[pos] = (true_key > approximate_key ? 0 : 1);
        }
        util::bit_compress(rsd_data);   
        rsd_symbol = sdsl::rrr_vector<63>(rsd_symbol_item);
    }

public:

     sdsl::int_vector<> rsd_data;
     sdsl::rrr_vector<63> rsd_symbol;

    SDSLRESIDUAL() = default;

    template <typename Index, typename Data>
    SDSLRESIDUAL(const Index& index, const Data& data)
         :rsd_data(data.size(), 0 , bits), rsd_symbol(){
        build(index, data, rsd_data, rsd_symbol);
    }

    size_t rsd_data_bytes() const { return size_in_bytes(rsd_data); }
    size_t rsd_symbol_bytes() const { return size_in_bytes(rsd_symbol); }

    // 对rsd_data进行合格性检查
    size_t checkResidual(size_t epsilon) const { 
        size_t size = rsd_data.size();
        size_t count = 0;

        for (size_t i = 0; i < size; ++i) {
            if (rsd_data[i] > epsilon + 1) {
                ++count;
                std::cout << "residual[" << i << "]:" << rsd_data[i] <<std::endl;
            }
        }

        return count;
    }

};
}