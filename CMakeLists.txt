 
cmake_minimum_required(VERSION 3.0)
project(MinimalistBT-CST)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O9")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")


include_directories(${sdsl_BINARY_DIR}/include)
include_directories(${sdsl_BINARY_DIR}/external/libdivsufsort/include)
include_directories(external/cst.jouni/libraries/TextIndexes/includes)
include_directories(external/cst.jouni/libraries/libcds2/includes)
include_directories(external/cst.jouni/libraries/tools/includes)
include_directories(external)
include_directories(includes)


add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(executables)