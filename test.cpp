#include <immintrin.h>
#include <iostream>
#include <cstring>
#include<cstdio>
#include <vector>
#include <algorithm>
using namespace std;
alignas(64) static uint32_t reverseshuffle[]={ 15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };

// one shuffle can be replaced with a simple blend, as min/max are commutative
static int blendmasks[] = {
    0b11111111'00000000, // L1->L2
    0b11110000'11110000, // L2->L3
    0b11001100'11001100, // L3->L4
    0b10101010'10101010, // L4->L5
};
// changed shuffles for blend+shuflle per level
#define L 0
#define H 16
alignas(64) static uint32_t shuffles2[][16]={
    { 8|L, 9|L,10|L,11|L,12|L,13|L,14|L,15|L,  0|H, 1|H, 2|H, 3|H, 4|H, 5|H, 6|H, 7|H }, // L1->L2 for H
    { 4|L, 5|L, 6|L, 7|L, 0|H, 1|H, 2|H, 3|H, 12|L,13|L,14|L,15|L, 8|H, 9|H,10|H,11|H }, // L2->L3 for H
    { 2|L, 3|L, 0|H, 1|H, 6|L, 7|L, 4|H, 5|H, 10|L,11|L, 8|H, 9|H,14|L,15|L,12|H,13|H }, // L3->L4 for H
    { 1|L, 0|H, 3|L, 2|H, 5|L, 4|H, 7|L, 6|H,  9|L, 8|H,11|L,10|H,13|L,12|H,15|L,14|H }, // L4->L5 for H

    { 0|L, 0|H, 1|L, 1|H, 2|L, 2|H, 3|L, 3|H,  4|L, 4|H, 5|L, 5|H, 6|L, 6|H, 7|L, 7|H }, // output first 16 elements
    {15|H,15|L,14|H,14|L,13|H,13|L,12|H,12|L, 11|H,11|L,10|H,10|L, 9|H, 9|L, 8|H, 8|L }
    // output for second reversed for next iteration
};
#undef L
#undef H

inline size_t union_u32_normal(const uint32_t *a, const uint32_t *b, size_t a_size, size_t b_size, uint32_t *out) {
    const uint32_t *a_end = a + a_size;
    const uint32_t *b_end = b + b_size;
    uint32_t *out_end = out;
    while (a != a_end && b != b_end) {
        bool le = (*a <= *b);
        bool ge = (*a >= *b);
        *out_end = le ? *a : *b;
        a += le;
        b += ge;
        out_end++;
    }
    std::memcpy(out_end, a, (a_end - a) * sizeof(uint32_t));
    out_end += a_end - a;
    std::memcpy(out_end, b, (b_end - b) * sizeof(uint32_t));
    out_end += b_end - b;
    return out_end - out;
}

inline size_t union_u32_simd(const uint32_t *list1, const uint32_t *list2, size_t size1, size_t size2, uint32_t *result) {
    // ... 你的原始代码（保留不变，但添加调试输出） ...
    size_t count = 0;
    size_t i_a = 0, i_b = 0;
    const size_t st_a = ((size1-1)/16)*16;
    const size_t st_b = ((size2-1)/16)*16;
    uint32_t a_nextfirst, b_nextfirst;
    __m512i circlic_shuffle = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __m512i old = _mm512_set1_epi32(-1); // FIXME: 硬编码
    alignas(64) uint32_t maxtail[16];
    alignas(64) int32_t temp[16];

    if (i_a < st_a && i_b < st_b) {
        // load all the shuffles
        __m512i vL1L2 = _mm512_load_epi32(shuffles2[0]);
        __m512i vL2L3 = _mm512_load_epi32(shuffles2[1]);
        __m512i vL3L4 = _mm512_load_epi32(shuffles2[2]);
        __m512i vL4L5 = _mm512_load_epi32(shuffles2[3]);
        __m512i vL5Out_L = _mm512_load_epi32(shuffles2[4]);
        __m512i vL5Out_H = _mm512_load_epi32(shuffles2[5]);


        __mmask16 kL1L2 = blendmasks[0];
        __mmask16 kL2L3 = blendmasks[1];
        __mmask16 kL3L4 = blendmasks[2];
        __mmask16 kL4L5 = blendmasks[3];

        __m512i vreverse = _mm512_load_epi32(reverseshuffle);

        __m512i v_a = _mm512_loadu_epi32(list1);
        __m512i vb = _mm512_loadu_epi32(list2);
        __m512i v_b = _mm512_permutexvar_epi32(vreverse, vb);

        do {
	    // bitonic merge network
            cerr << "New Round" << endl;
            _mm512_store_epi32(temp, v_a);
            cerr << "v_a: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl << "v_b: ";
            _mm512_store_epi32(temp, v_b);
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl;
	    // level 1
	    __m512i min = _mm512_min_epi32(v_a, v_b); // find minimum of both vectors in each position
	    __m512i max = _mm512_max_epi32(v_a, v_b); // find maximum of both vectors in each position
	    __m512i L = _mm512_mask_blend_epi32(kL1L2, min, max); // blend minimum and maximum
	    __m512i H = _mm512_permutex2var_epi32(min, vL1L2, max);
            _mm512_store_epi32(temp, L);
            cerr << "Level 1 Lower: ";
            for (int i=0; i<16; i++) {
                    cerr << temp[i] << " ";
            }
            cerr << endl;
            _mm512_store_epi32(temp, H);
            cerr << "Higher: ";
            for (int i=0; i<16; i++) {
                    cerr << temp[i] << " ";
            }
            cerr << endl << endl;

	    // level 2
	    min = _mm512_min_epi32(L, H);
	    max = _mm512_max_epi32(L, H);
	    L = _mm512_mask_blend_epi32(kL2L3, min, max);
	    H = _mm512_permutex2var_epi32(min, vL2L3, max);
            _mm512_store_epi32(temp, L);
            cerr << "Level 2 Lower: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl;
            _mm512_store_epi32(temp, H);
            cerr << "Higher: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl << endl;
	    // level 3
	    min = _mm512_min_epi32(L, H);
	    max = _mm512_max_epi32(L, H);
	    L = _mm512_mask_blend_epi32(kL3L4, min, max);
	    H = _mm512_permutex2var_epi32(min, vL3L4, max);
            _mm512_store_epi32(temp, L);
            cerr << "Level 3 Lower: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl;
            _mm512_store_epi32(temp, H);
            cerr << "Higher: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl << endl;
	    // level 4
	    min = _mm512_min_epi32(L, H);
	    max = _mm512_max_epi32(L, H);
	    L = _mm512_mask_blend_epi32(kL4L5, min, max);
	    H = _mm512_permutex2var_epi32(min, vL4L5, max);
            _mm512_store_epi32(temp, L);
            cerr << "Level 4 Lower: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl;
            _mm512_store_epi32(temp, H);
            cerr << "Higher: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl << endl;
	    // level 5
	    min = _mm512_min_epi32(L, H);
	    max = _mm512_max_epi32(L, H);
	    L = _mm512_permutex2var_epi32(min, vL5Out_L, max);
	    H = _mm512_permutex2var_epi32(min, vL5Out_H, max);
            _mm512_store_epi32(temp, L);
            cerr << "Level 5 Lower: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl;
            _mm512_store_epi32(temp, H);
            cerr << "Higher: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl << endl;

	    __m512i recon = _mm512_mask_blend_epi32(0x7fff, old, L);
            _mm512_store_epi32(temp, recon);
            cerr << "Recon: ";
            for (int i=0; i<16; i++) {
                cerr << temp[i] << " ";
            }
            cerr << endl << endl;

	    __m512i dedup = _mm512_permutexvar_epi32(circlic_shuffle, recon);
	    __mmask16 kmask = _mm512_cmpneq_epi32_mask(dedup, L);
	    _mm512_mask_compressstoreu_epi32(&result[count], kmask, L);
            for (int i = count; i < count + _mm_popcnt_u32(kmask); i++) {
                cerr << result[i] << " ";
            }
	    count += _mm_popcnt_u32(kmask); // number of elements written
            cerr << endl << _mm_popcnt_u32(kmask) << endl;
	    // remember minimum for next iteration
	    old = L;

	    v_a = H;
	    // compare first element of the next block in both lists
	    a_nextfirst = list1[i_a+16];
	    b_nextfirst = list2[i_b+16];
	    // write minimum as above out to result
	    // keep maximum and do the same steps as above with next block
	    // next block from one list, which first element in new block is smaller
	    bool a_next = (a_nextfirst <= b_nextfirst);
	    i_a += a_next * 16;
	    i_b += !a_next * 16;
	    size_t index = a_next ? i_a: i_b;
	    const uint32_t *base = a_next ? list1: list2;
	    v_b = _mm512_loadu_epi32(&base[index]);
        } while (i_a < st_a && i_b < st_b);
        _mm512_store_epi32(maxtail, _mm512_permutexvar_epi32(vreverse, v_a));

        uint32_t endofblock = _mm_extract_epi32(_mm512_extracti32x4_epi32(old, 3), 3);

        size_t mti = 0;
        size_t mtsize = std::unique(maxtail, maxtail+16) - maxtail; // deduplicate tail
        if(a_nextfirst <= b_nextfirst){
            // endofblock needs to be considered too, for deduplication
            if(endofblock == std::min(maxtail[0],list1[i_a])) --count;
            // compare maxtail with list1
            while(mti < mtsize && i_a < size1){
                if(maxtail[mti] < list1[i_a]){
                    result[count++] = maxtail[mti];
                    mti++;
                }else if(maxtail[mti] > list1[i_a]){
                    result[count++] = list1[i_a];
                    i_a++;
                }else{
                    result[count++] = maxtail[mti];
                    mti++; i_a++;
                }
            }
            i_b += 16;
        }else{
            // endofblock needs to be considered too, for deduplication
            if(endofblock == std::min(maxtail[0],list2[i_b]))
                --count;
            // compare maxtail with list2
            while(mti < mtsize && i_b < size2){
                if(maxtail[mti] < list2[i_b]){
                    result[count++] = maxtail[mti];
                    mti++;
                }else if(maxtail[mti] > list2[i_b]){
                    result[count++] = list2[i_b];
                    i_b++;
                }else{
                    result[count++] = maxtail[mti];
                    mti++; i_b++;
                }
            }
            i_a += 16;
        }
        while(mti < mtsize){
            result[count++] = maxtail[mti++];
        }
    }

    count += union_u32_normal(list1+i_a, list2+i_b, size1-i_a, size2-i_b, result+count);

    return count;
}

int main() {
    // 测试用例：两个已排序数组
    uint32_t list1[200];
    uint32_t list2[200];
    for (size_t i=0; i<200; i++) {
        list1[i] = i * 2;
        list2[i] = i * 3;
    }
    size_t size1 = sizeof(list1)/sizeof(list1[0]);
    size_t size2 = sizeof(list2)/sizeof(list2[0]);

    uint32_t result[400];
    size_t count = union_u32_simd(list1, list2, size1, size2, result);

    std::cout << "\nCount: " << count << std::endl;

    std::cout << "Result: " << union_u32_normal(list1, list2, size1, size2, result) << std::endl;

    return 0;
}