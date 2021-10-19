add_executable(CountyLockDown CountyLockDown.cpp)
add_executable(CountryLockDown CountryLockDown.cpp)
add_executable(MinImpactLDZones MinImpactLDZones.cpp)
add_executable(CovidSimTest CovidSim_test.cpp)
add_executable(TwoNode Test.cpp)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(CountyLockDown EpiGraph)
target_link_libraries(CountryLockDown EpiGraph)
target_link_libraries(MinImpactLDZones EpiGraph)
target_link_libraries(CovidSimTest EpiGraph)
target_link_libraries(TwoNode EpiGraph)

