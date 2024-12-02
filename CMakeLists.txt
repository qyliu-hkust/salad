cmake_minimum_required(VERSION 2.8)
project(PGM)
#add_compile_options(-mavx2)
add_compile_options(-mavx512f)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_STANDARD 17)
include_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/include)
link_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/lib)
#
#link_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/include/sdsl)
#
include_directories(${PROJECT_SOURCE_DIR}/include)
#add_subdirectory(${PROJECT_SOURCE_DIR}/external)
#add_subdirectory(${PROJECT_SOURCE_DIR}/include/sdsl)

#echo "${PROJECT_SOURCE_DIR}/include/sdsl/include"
#include_directories(${PROJECT_SOURCE_DIR}/include/sdsl/include/sdsl)
#find_library(SDSL_LIB sdsl)

#add_executable(pgm pgm.cpp)
#target_link_libraries(pgm sdsl)

add_executable(pgm_build pgm_build.cpp)
target_link_libraries(pgm_build sdsl)

add_executable(pgm_decode pgm_decode.cpp)
target_link_libraries(pgm_decode sdsl)

#add_executable(rsd_huffman rsd_huffman.cpp)
#target_link_libraries(rsd_huffman sdsl)

#add_executable(test test.cpp)
#target_link_libraries(test sdsl)
#target_link_libraries(pgm sdsl)
#./install.sh /mnt/home/xyzhu/libs/sdsl-lite
#./configure --prefix="/mnt/home/xyzhu/libs/python3.7/"
#alias python=/mnt/home/xyzhu/libs/python3.7/bin/python3
#alias pip=/mnt/home/xyzhu/libs/python3.7/bin/pip3