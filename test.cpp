#include <iostream>
#include <vector>
#include <sys/mman.h>
#include <cstring>
#include <stdexcept>
#include <cerrno>
#include <string>
#include <unistd.h>

// 自定义分配器
template <typename T>
class HugePageAllocator {
    constexpr static size_t PAGE_SIZE = 2 * 1024 * 1024;
public:
    using value_type = T;

    // 构造函数和析构函数
    HugePageAllocator() = default;
    ~HugePageAllocator() = default;

    // 分配内存
    T* allocate(std::size_t n) {
        const size_t size = n * sizeof(T);
        void* ptr = mmap(nullptr, size,
                         PROT_READ | PROT_WRITE,
                         MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                         -1, 0);
        if (ptr == MAP_FAILED) {
            std::cerr << "Failed to allocate huge pages: " << std::strerror(errno) << std::endl;
            throw std::bad_alloc();
        }
        return static_cast<T*>(ptr);
    }

    // 释放内存
    void deallocate(T* p, std::size_t n) noexcept {
        if (munmap(p, PAGE_SIZE) == -1) {
            std::cerr << "Failed to unmap huge pages: " << std::strerror(errno) << std::endl;
        }
    }

    // 比较操作符（必须定义）
    bool operator == (const HugePageAllocator&) const { return true; }
    bool operator != (const HugePageAllocator&) const { return false; }
};

int main() {
    const size_t hugePageSize = 2 * 1024 * 1024; // 根据 Hugepagesize 设置
    const size_t n = 1000000; // 确保是大页的整数倍

    try {
        // 创建一个使用大页内存的 vector
        std::vector<int, HugePageAllocator<int>> current_value_vector;

        // 调整大小
        current_value_vector.resize(n);
        int *current_value = current_value_vector.data();

        // 测试写入和读取
        for (size_t i = 0; i < n; ++i) {
            current_value[i] = i;
        }

        // 验证数据
        for (size_t i = 0; i < 16; ++i) {
            std::cout << current_value[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Successfully used huge pages with std::vector!" << std::endl;
        while (true) {
            sleep(1);
        }

    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << std::endl;
    }

    return 0;
}