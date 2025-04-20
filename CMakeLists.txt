cmake_minimum_required(VERSION 3.30)
project(BPMsolver)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_FLAGS "-lhdf5")

find_package(Eigen3 CONFIG REQUIRED)
find_package(Armadillo CONFIG REQUIRED)
find_package(hdf5 CONFIG REQUIRED)
find_package(superlu CONFIG REQUIRED)


add_executable(BPMsolver main.cpp
        function/commom.h
        function/Thomas.cpp
        function/Thomas.h
        MatrixOperator/dxdx/dxdxbase.cpp
        MatrixOperator/dxdx/dxdxbase.h
        MatrixOperator/dydy/dydybase.cpp
        MatrixOperator/dydy/dydybase.h
        MatrixOperator/dxdy/dxdybase.cpp
        MatrixOperator/dxdy/dxdybase.h
        MatrixOperator/dydx/dydxbase.cpp
        MatrixOperator/dydx/dydxbase.h
        MatrixOperator/test.h
        MatrixOperator/dxdx/dxdx.cpp
        MatrixOperator/dxdx/dxdx.h
        MatrixOperator/dxdy/dxdy.cpp
        MatrixOperator/dxdy/dxdy.h
        MatrixOperator/dydx/dydx.cpp
        MatrixOperator/dydx/dydx.h
        MatrixOperator/dydy/dydy.cpp
        MatrixOperator/dydy/dydy.h
        function/loadHDF5.cpp
        function/loadHDF5.h
        sim/BPM.cpp
        sim/BPM.h
        MatrixOperator/matrixHeader.h
        sim/bpmOneSetp.cpp
        sim/bpmOneSetp.h
)

target_link_libraries(BPMsolver PRIVATE armadillo)
target_link_libraries(BPMsolver PRIVATE Eigen3::Eigen)
target_link_libraries(BPMsolver PRIVATE hdf5::hdf5-shared hdf5::hdf5_hl-shared)
target_link_libraries(BPMsolver PRIVATE superlu::superlu)