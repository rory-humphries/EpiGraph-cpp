cmake_minimum_required(VERSION 3.13)
project(EpiGraph-cpp)
#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

# Configure library
file(GLOB_RECURSE epigraph_src "epigraph/*.hpp" "epigraph/*.cpp")
add_library(epigraph STATIC ${epigraph_src} epigraph/Core/Util.hpp)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(epigraph Eigen3::Eigen)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(epigraph OpenMP::OpenMP_CXX)
endif ()

target_include_directories(epigraph PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/libs)

# configure executables
add_subdirectory(simulations)
add_subdirectory(simulations/test_dir)




