cmake_minimum_required(VERSION 3.15)
project(EpiGraph)

set(CMAKE_CXX_STANDARD 14)

add_executable(EpiGraph SIRX_ed.cpp edge.h spatial_graph.h distributions.h algos.h)