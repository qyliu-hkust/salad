cmake_minimum_required(VERSION 2.8)
project(PGM)
#add_compile_options(-mavx2)
add_compile_options(-mlzcnt)
#add_compile_options(-msse4.2)
add_compile_options(-mavx512f)
add_compile_options(-mavx512vl)  # Add this line
add_compile_options(-mavx512cd) # Add this line
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_STANDARD 17)

set(FASTPFOR_LIB_PATH ${PROJECT_SOURCE_DIR}/external/FastPFor/build)
set(FASTPFOR_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/external/FastPFor/headers)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external/Vp2SetOperation/include)

# sdsl-lite
include_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/include)
link_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/lib)

# FastPFor
include_directories(${FASTPFOR_INCLUDE_PATH})
link_directories(${FASTPFOR_LIB_PATH})

add_executable(pgm_build pgm_build.cpp)
target_link_libraries(pgm_build sdsl FastPFOR)

add_executable(pgm_decode pgm_decode.cpp)
target_link_libraries(pgm_decode sdsl)

add_executable(pgm_query pgm_query.cpp)
target_link_libraries(pgm_query sdsl)

#add_executable(pgm_residual_compressor pgm_residual_compressor.cpp)
#target_link_libraries(pgm_residual_compressor sdsl)
#
add_executable(test test.cpp)
target_link_libraries(test sdsl FastPFOR)