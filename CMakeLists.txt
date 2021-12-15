cmake_minimum_required(VERSION 3.19)
project(PHA_kHNDP)

set(CMAKE_CXX_STANDARD 11)

add_executable(PHA_kHNDP main.cpp src/graph.cpp src/functions.cpp)

include_directories("/usr/local/include")