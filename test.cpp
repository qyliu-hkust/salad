#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>

// 将7位数存入int_vector<64>
void store_7bit_numbers(const std::vector<uint8_t>& numbers, sdsl::int_vector<64>& vec) {
    // 计算总比特数和所需块数
    uint64_t total_bits = numbers.size() * 7;
    uint64_t required_blocks = (total_bits + 63) / 64; // 向上取整
    vec.resize(required_blocks);
    vec.width(64); // 确保每个元素为64位

    // 获取底层内存指针
    uint64_t* data = vec.data();
    uint64_t offset = 0; // 当前写入的比特偏移量

    // 遍历所有数，逐个写入
    for (auto num : numbers) {
        // 确保数值不超过7位（0~127）
        uint64_t value = num & 0x7F;
        // 使用SDSL函数写入7位
        sdsl::bits::write_int(data, value, offset, 7);
        offset += 7; // 更新偏移量
    }
}

// 从int_vector<64>中读取第index个7位数
uint8_t read_7bit_number(const sdsl::int_vector<64>& vec, uint64_t index) {
    uint64_t offset = index * 7; // 计算该数的起始比特位置
    return sdsl::bits::read_int(vec.data(), offset, 7); // 读取7位
}

template <typename T>
T* aligned_new(uint64_t num_elements) {
    void* ptr = std::aligned_alloc(64, num_elements * sizeof(T));
    if (!ptr) throw std::bad_alloc();
    return static_cast<T*>(ptr);
}

int main() {
    // 示例：存储10个7位数（0x00~0x7F）
    __m512i rerange_idx = _mm512_set_epi32(240, 224, 208, 192, 176, 160,  144,  128, 112, 96, 80,  64, 48, 32, 16, 0);
    // __m512i rerange_idx = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11,  12, 13, 14, 15);
    std::vector<int> a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<int> b (100);
    alignas(64) int *output = aligned_new<int>(16);
    __m512i result_v = _mm512_loadu_epi32(a.data());
    __m128i t0 = _mm512_extracti32x4_epi32(result_v, 0);
    __m128i t1 = _mm512_extracti32x4_epi32(result_v, 1);
    __m128i t2 = _mm512_extracti32x4_epi32(result_v, 2);
    __m128i t3 = _mm512_extracti32x4_epi32(result_v, 3);

    std::vector<int*> first_pointer(16);
    first_pointer[0] = b.data();

    *((int64_t*)(first_pointer[0])) = _mm_extract_epi64(t0, 0);  // a0, a1
    first_pointer[0] += 2;
    *((int64_t*)(first_pointer[0])) = _mm_extract_epi64(t0, 1);  // a0, a1

    for (int i = 0; i < 16; i++) {
        std::cout << b[i] << std::endl;
    }
    return 0;
}