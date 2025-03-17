#include <iostream>
#include <vector>
#include <list>
#include <functional>

// namespace pgm_sequence {
    template <typename KeyType, typename ValueType>
    class HashTable {
    public:

        HashTable() = default;

        HashTable(size_t size) : table(size) {}

        void insert(const KeyType& key, const ValueType& value) {
            size_t hashValue = hashFunction(key);
            while (hashValue > table.size()) hashValue -= table.size();
            for (auto& pair : table[hashValue]) {
                if (pair.first == key) {
                    pair.second = value;
                    return;
                }
            }
            table[hashValue].emplace_back(key, value);
        }

        bool find(const KeyType& key, ValueType& value) const {
            size_t hashValue = hashFunction(key);
            while (hashValue > table.size()) hashValue -= table.size();
            for (const auto& pair : table[hashValue]) {
                if (pair.first == key) {
                    value = pair.second;
                    return true;
                }
            }
            return false;
        }

        inline ValueType find_value (KeyType key) {
            size_t hashValue = hashFunction(key);
            while (hashValue > table.size()) hashValue -= table.size();
            for (const auto& pair : table[hashValue]) {
                if (pair.first == key) {
                    return pair.second;
                }
            }
            std::cerr << "Error: Key not found" << std::endl;
            return 0;
        }

        // inline ValueType operator [](const KeyType& key) {
        //     size_t hashValue = hashFunction(key);
        //     while (hashValue > table.size()) hashValue -= table.size();
        //     for (const auto& pair : table[hashValue]) {
        //         if (pair.first == key) {
        //             return pair.second;
        //         }
        //     }
        //     std::cerr << "Error: Key not found" << std::endl;
        //     return 0;
        // }

        void remove(const KeyType& key) {
            size_t hashValue = hashFunction(key) % table.size();
            table[hashValue].remove_if([&key](const std::pair<KeyType, ValueType>& pair) {
                return pair.first == key;
            });
        }

    private:
        std::vector<std::list<std::pair<KeyType, ValueType>>> table;
        std::hash<KeyType> hashFunction;
    };
// }
