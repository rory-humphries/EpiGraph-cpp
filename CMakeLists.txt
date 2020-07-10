cmake_minimum_required(VERSION 3.15)
project(EpiGraph)

set(CMAKE_CXX_STANDARD 14)

add_executable(MultiPhase     include/edge.h include/sixrd_model.h include/graph.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_utils.h include/meta_pop_network.h include/io.h include/travel_model.h src/multi_phase.cpp)
add_executable(VertexLockDown include/edge.h include/sixrd_model.h include/graph.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_utils.h include/meta_pop_network.h include/io.h include/travel_model.h src/vertex_dyn_lockdown.cpp)
add_executable(CountyLockDown include/edge.h include/sixrd_model.h include/graph.h include/distributions.h include/csv.h src/distributions.cpp include/spatial_utils.h include/meta_pop_network.h include/io.h include/travel_model.h src/county_dyn_lockdown.cpp)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()