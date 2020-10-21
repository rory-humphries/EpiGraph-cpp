cmake_minimum_required(VERSION 3.13)
project(EpiGraph-cpp)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

# Configure library

add_library(epigraph STATIC
        epigraph/sixrd_odes.hpp epigraph/distributions.cpp epigraph/distributions.hpp
        epigraph/eigen_util.hpp epigraph/io.h
        epigraph/random_matrix.hpp
        epigraph/sixrd_net_model.cpp epigraph/sixrd_net_model.hpp
        epigraph/sixrd_net_model_up.cpp epigraph/sixrd_net_model_up.hpp
        epigraph/spatial_util.cpp epigraph/spatial_util.hpp
        epigraph/travel_strategies.hpp)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(epigraph Eigen3::Eigen)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(epigraph OpenMP::OpenMP_CXX)
endif ()

target_include_directories(epigraph PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/libs)

# configure executables

add_executable(CountyDynLockDown simulations/dyn_county_lockdown.cpp)
add_executable(MultiPhase simulations/multi_phase_lockdown.cpp)


find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(CountyDynLockDown Eigen3::Eigen epigraph)
target_link_libraries(MultiPhase Eigen3::Eigen epigraph)


find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(CountyDynLockDown OpenMP::OpenMP_CXX)
    target_link_libraries(MultiPhase OpenMP::OpenMP_CXX)
endif ()


