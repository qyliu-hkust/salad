#include <chrono>
#include <vector>
#include <iostream>
#include <memory>
#include <immintrin.h>
#include "sdsl/int_vector.hpp"
// #include <aligned_storage>

using namespace std;

// struct AlignedFloatArray {
//     alignas(32) float data[8];
// };
typedef int K;
const int align_val = 64;
int key_nums = 8;
int segments_size_simd = 2;
std::vector<float *> vec;
std::vector<float *> slope;
std::vector<float *> intercept;
std::vector<float *> signs;
std::vector<K *> covered;
std::vector<K *> first;
sdsl::int_vector<64> corrections;

template <typename T>
T* aligned_new(size_t num_elements) {
    void* ptr = std::aligned_alloc(align_val, num_elements * sizeof(T));
    if (!ptr) throw std::bad_alloc();
    return static_cast<T*>(ptr);
}

template <typename T>
void aligned_delete(T* ptr) {
    std::free(ptr);
}

// # pragma pack(4)
void add_float(int k) {
    // alignas(32) float* arr1 = new float[]{1,2,3,4,5,6,7,8};
    // alignas(32) float* arr1 = new float[]{5,6,7,8,9,10,11,12};
    // // cout << arr1 << endl;
    // vec.emplace_back(arr1);
    // alignas(32) float* arr2 = new float[]{1,2,3,4,5,6,7,8};
    // vec.emplace_back(arr2);
    // // AlignedFloatArray arr1 = AlignedFloatArray();
    K* covered_tmp = new K[key_nums];
    K* first_tmp = new K[key_nums];
    // std::vector<Segment> simd_segments(it, it + key_nums);
    float* slope_simd_tmp = aligned_new<float>(key_nums);
    float* intercept_simd_tmp = aligned_new<float>(key_nums);
    float* sings_simd_tmp = aligned_new<float>(key_nums);
    for (int i = 0;i < key_nums;i++) {
        covered_tmp[i] = i + 1;
        first_tmp[i] = i + 1;
        slope_simd_tmp[i] = i + 1;
        intercept_simd_tmp[i] = i + 1;
    }
    // covered[k] = covered_tmp;
    // first[k] = first_tmp;
    // slope[k] = slope_simd_tmp;
    // intercept[k] = intercept_simd_tmp;
    // signs[k] = sings_simd_tmp;
    covered.emplace_back(covered_tmp);
    first.emplace_back(first_tmp);
    slope.emplace_back(slope_simd_tmp);
    intercept.emplace_back(intercept_simd_tmp);

    // cout << & arr1 << endl;
    //
    // AlignedFloatArray arr2 = AlignedFloatArray();
    // for (int i = 0;i < 8;i++) {
    //     arr2.data[i] = i;
    // }
    // cout << & arr2 << endl;
    //
    // vec.emplace_back(arr1);
    // vec.emplace_back(arr2);
}
// # pragma pack()

int main() {
    // 添加一些浮点数数组到vector
    // slope.reserve(segments_size_simd);
    // intercept.reserve(segments_size_simd);
    // covered.reserve(segments_size_simd);
    // signs.reserve(segments_size_simd);
    // first.reserve(segments_size_simd);
    // 创建一个 sdsl::int_vector<64>，初始容量为 5，所有元素初始化为 -1
    sdsl::int_vector<32> corrections(1, 0);

    // 设置第一个元素为 -10
    corrections[0] = -1;

    // 输出第一个元素
    std::cout << "First element: " << corrections[0] << std::endl;

    return 0;

    for (int i = 0; i < segments_size_simd; i++) {
        add_float(i);
    }

    float* slopev = slope[0];
    float* interceptv = intercept[0];
    float* result = aligned_new<float>(key_nums);
    uint64_t avg_count = 0;

    int exp_num = 500000000;
    double *sloped = aligned_new<double>(key_nums);
    double *interceptd = aligned_new<double>(key_nums);
    double *resultd = aligned_new<double>(key_nums);


    for (int i = 0; i < exp_num; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        __m256 slope_v = _mm256_load_ps(slopev); // 加载斜率
        __m256 intercept_v = _mm256_load_ps(interceptv); // 加载斜率
        __m256 result_v = _mm256_add_ps(slope_v, intercept_v);
        _mm256_store_ps(result, result_v);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        avg_count += duration.count();
        if (i >= exp_num - 5) {
            cout << duration.count() << endl;
        }
    }
    cout << avg_count / exp_num << endl;
    cout << endl;
    avg_count = 0;
    for (int64_t i = 0; i < key_nums; i++) {
        slopev[i] = i * i;
        interceptv[i] = i * i;
    }


    for (int i = 0; i < exp_num; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        __m512 slove_v512 = _mm512_load_ps(slopev);
        __m512 intercept_v512 = _mm512_load_ps(interceptv);
        __m512 result_v512 = _mm512_add_ps(slove_v512, intercept_v512);
        _mm512_store_ps(result, result_v512);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        avg_count += duration.count();
        if (i >= exp_num - 5) {
            cout << duration.count() << endl;
        }
    }

    cout << avg_count / exp_num << endl;
    cout << endl;

    avg_count = 0;
    for (int64_t i = 0; i < key_nums; i++) {
        sloped[i] = i * i;
        interceptd[i] = i * i;
    }


    for (int i = 0; i < exp_num; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        __m512d slove_v512 = _mm512_load_pd(sloped);
        __m512d intercept_v512 = _mm512_load_pd(interceptd);
        __m512d result_v512 = _mm512_add_pd(slove_v512, intercept_v512);
        _mm512_store_pd(resultd, result_v512);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        avg_count += duration.count();
        if (i >= exp_num - 5) {
            cout << duration.count() << endl;
        }
    }

    cout << avg_count / exp_num << endl;
    cout << endl;


    int64_t* slopei = aligned_new<int64_t>(key_nums * 3);
    int64_t *intercepti = aligned_new<int64_t>(key_nums * 3);
    int64_t *resulti = aligned_new<int64_t>(key_nums * 3);
    for (int64_t i = 0; i < key_nums * 3; i++) {
        slopei[i] = i * i;
        intercepti[i] = i * i;
    }

    avg_count = 0;

    for (int i = 0; i < exp_num; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        __m512i slove_v512 = _mm512_load_epi64(slopei + key_nums);
        __m512i intercept_v512 = _mm512_load_epi64(intercepti + key_nums);
        __m512i result_v512 = _mm512_add_epi64(slove_v512, intercept_v512);
        _mm512_store_epi64(resulti, result_v512);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        avg_count += duration.count();
        if (i >= exp_num - 5) {
            cout << duration.count() << endl;
        }
    }

    cout << avg_count / exp_num << endl;


    // cout << endl;
    // for (int i = 0; i < key_nums; i++) {
    //     cout << result[i] << endl;
    // }

    // 由于使用了unique_ptr，内存会在vector析构时自动释放

    return 0;
}