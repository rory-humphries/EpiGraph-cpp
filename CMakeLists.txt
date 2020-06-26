cmake_minimum_required(VERSION 3.15)
project(EpiGraph)

set(CMAKE_CXX_STANDARD 14)

add_executable(EpiGraph edge.h sixrd_model.h graph.h distributions.h csv.h sixrd_ed_soa.cpp distributions.cpp spatial_utils.h meta_pop_network.h io.h travel_model.h)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()