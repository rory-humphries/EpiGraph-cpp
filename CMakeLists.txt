cmake_minimum_required(VERSION 3.13)
project(EpiGraph-cpp)
#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# Configure library
file(GLOB_RECURSE epigraph_src "EpiGraph/*.hpp" "EpiGraph/*.cpp")
add_library(EpiGraph SHARED ${epigraph_src})

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(EpiGraph Eigen3::Eigen)
#target_precompile_headers(EpiGraph PUBLIC EpiGraph/Core/NetEpiComp.hpp)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(EpiGraph OpenMP::OpenMP_CXX)
endif ()

target_include_directories(EpiGraph PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/libs ${CMAKE_CURRENT_LIST_DIR}/libs/spectra/include/)

# configure executables
add_subdirectory(simulations)



