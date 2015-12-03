cmake_minimum_required(VERSION 2.6)

project(myfdtd)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Wextra -Wall -std=c++11")

add_executable(myfdtd medium.cpp solver.cpp util.cpp main.cpp)

install(TARGETS myfdtd RUNTIME DESTINATION bin)
