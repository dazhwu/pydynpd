cmake_minimum_required(VERSION 3.9)
include(FetchContent)
project(ar_test)
# Set C++ 17 Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(mkl_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_CXX_FLAGS "-Wall -O3 -mavx2 -fopenmp")


FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v2.9.2
)

FetchContent_GetProperties(pybind11)
if (NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif ()


find_package(PkgConfig)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_subdirectory(pybind11)
pybind11_add_module(ar_test AR.cpp)


