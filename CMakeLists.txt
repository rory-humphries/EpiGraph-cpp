cmake_minimum_required(VERSION 3.15)
project(EpiGraph)

set(CMAKE_CXX_STANDARD 14)

include_directories("/usr/include/eigen3/")

add_executable(MultiPhase     include/sixrd_model.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_util.h include/io.h include/travel_model.h src/multi_phase.cpp include/eigen_util.h)
add_executable(VertexLockDown include/sixrd_model.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_util.h include/io.h include/travel_model.h src/vertex_dyn_lockdown.cpp)
add_executable(CountyLockDown include/sixrd_model.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_util.h include/io.h include/travel_model.h src/county_dyn_lockdown.cpp)

find_package (Eigen3 REQUIRED NO_MODULE)
target_link_libraries (MultiPhase Eigen3::Eigen)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MultiPhase OpenMP::OpenMP_CXX )
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

