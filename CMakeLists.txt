cmake_minimum_required(VERSION 3.15)
project(EpiGraph)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories("/usr/include/eigen3/")

add_executable(MultiPhaseLockDown include/sixrd_model.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_util.h include/io.h include/travel_model.h src/multi_phase_lockdown.cpp include/eigen_util.h)
add_executable(DynLockDown include/sixrd_model.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_util.h include/io.h include/travel_model.h src/dyn_lockdown.cpp)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(MultiPhaseLockDown Eigen3::Eigen)
target_link_libraries(DynLockDown Eigen3::Eigen)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(MultiPhaseLockDown OpenMP::OpenMP_CXX)
    target_link_libraries(DynLockDown OpenMP::OpenMP_CXX)
endif ()


