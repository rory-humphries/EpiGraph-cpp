cmake_minimum_required(VERSION 3.15)
project(EpiGraph-cpp)

set(CMAKE_CXX_STANDARD 17)

include_directories("/usr/include/eigen3/")

add_executable(MultiPhaseLockDown
        epigraph/multi_phase_lockdown.cpp
        epigraph/sixrd_model.h epigraph/sixrd_model.cpp
        epigraph/distributions.h epigraph/distributions.cpp
        epigraph/spatial_util.h epigraph/spatial_util.cpp)

add_executable(DynLockDown
        epigraph/dyn_lockdown.cpp
        epigraph/sixrd_model.h epigraph/sixrd_model.cpp
        epigraph/distributions.h epigraph/distributions.cpp
        epigraph/spatial_util.h epigraph/spatial_util.cpp
        epigraph/io.h
        epigraph/random_matrix.h)

add_executable(CountyDynLockDown
        epigraph/dyn_county_lockdown.cpp
        epigraph/sixrd_model.h epigraph/sixrd_model.cpp
        epigraph/sixrd_odes.h
        epigraph/distributions.h epigraph/distributions.cpp
        epigraph/spatial_util.h epigraph/spatial_util.cpp
        epigraph/io.h
        epigraph/random_matrix.h)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(MultiPhaseLockDown Eigen3::Eigen)
target_link_libraries(DynLockDown Eigen3::Eigen)
target_link_libraries(CountyDynLockDown Eigen3::Eigen)


find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(MultiPhaseLockDown OpenMP::OpenMP_CXX)
    target_link_libraries(DynLockDown OpenMP::OpenMP_CXX)
    target_link_libraries(CountyDynLockDown OpenMP::OpenMP_CXX)

endif ()


