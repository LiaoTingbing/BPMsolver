cmake_minimum_required(VERSION 3.30)
project(BPMsolver)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 CONFIG REQUIRED)
find_package(Armadillo CONFIG REQUIRED)
add_executable(BPMsolver main.cpp)

target_link_libraries(BPMsolver PRIVATE armadillo)
target_link_libraries(BPMsolver PRIVATE Eigen3::Eigen)
