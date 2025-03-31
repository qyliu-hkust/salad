#pragma once

#include <sdsl.hpp>
#include <iostream>
#include <string>
#include <vector>


namespace rsd{
template<typename K>
class RESIDUAL {
protected:

    template <typename Index, typename Data>
    static void build(Index index, Data data, 
                        sdsl::int_vector<> &rsd_data,
                        sdsl::bit_vector &rsd_symbol) {

        for (size_t i = 0; i < data.size(); i++){
            K true_key = data[i];  
            K approximate_key = index.search(i);
            rsd_data[i] = (true_key > approximate_key ? true_key - approximate_key : approximate_key - true_key);
            rsd_symbol[i] = (true_key > approximate_key ? 0 : 1);
        }

        sdsl::util::bit_compress(rsd_data);
    }

public:

    sdsl::int_vector<> rsd_data;  
    sdsl::bit_vector rsd_symbol;  
    

    RESIDUAL() = default;

    template <typename Index, typename Data>
    RESIDUAL(const Index &index, const Data &data)
         : rsd_data(data.size()),
          rsd_symbol(data.size()){
        build(index, data, rsd_data, rsd_symbol);
    }

    size_t data_bytes() const { return rsd_data.size() * rsd_data.width() / 8; }
    size_t symbol_bytes() const { return rsd_symbol.size() * rsd_symbol.width() / 8; }
};
}





















// #pragma once

// #include <sdsl.hpp>
// #include <iostream>
// #include <string>
// #include <vector>


// namespace rsd{
// template<typename K>
// class RESIDUAL {
// protected:

//     template <typename Index, typename Data>
//     static void build(Index index, Data data, 
//                         sdsl::int_vector<> &rsd_data,
//                         sdsl::int_vector<> &rsd_symbol) {

//         for (size_t i = 0; i < data.size(); i++){
//             K true_key = data[i];  
//             K approximate_key = index.search(i);
//             rsd_data[i] = (true_key > approximate_key ? true_key - approximate_key : approximate_key - true_key);
//             rsd_symbol[i] = (true_key > approximate_key ? 0 : 1);
//         }
//     }

// public:

//     sdsl::int_vector<> rsd_data;  
//     sdsl::int_vector<> rsd_symbol;  
    

//     RESIDUAL() = default;

//     template <typename Index, typename Data>
//     RESIDUAL(const Index &index, const Data &data)
//          : rsd_data(data.size()),
//           rsd_symbol(data.size()){
//         build(index, data, rsd_data, rsd_symbol);
//     }

//     size_t data_bytes() const { return rsd_data.size() * rsd_data.width() / 8; }
//     size_t symbol_bytes() const { return rsd_symbol.size() * rsd_symbol.width() / 8; }
// };
// }