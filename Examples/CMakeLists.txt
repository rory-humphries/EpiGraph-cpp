link_libraries(EpiGraph)

add_executable(connected_component connected_component.cpp)
add_executable(degree degree.cpp)
add_executable(sir sir.cpp)
add_executable(si si.cpp)

#target_link_libraries(sir EpiGraph)
#target_link_libraries(si EpiGraph)
