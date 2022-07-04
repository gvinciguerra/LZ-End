cmake_minimum_required(VERSION 3.12)
project(LZEnd)

set(CMAKE_CXX_STANDARD 17)

include_directories(lib/sdsl-lite/include include)

add_executable(example example.cpp)
