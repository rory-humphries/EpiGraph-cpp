cmake_minimum_required(VERSION 3.15)
project(EpiGraph)

set(CMAKE_CXX_STANDARD 14)

add_executable(EpiGraph edge.h sixrd_model.h graph.h graph.cpp distributions.h csv.h sixrd_ed_soa.cpp distributions.cpp spatial_utils.h meta_pop_network.h network_builder.h io.h travel_model.h)