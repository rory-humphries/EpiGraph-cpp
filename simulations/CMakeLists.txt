add_executable(CountyLockDown CountyLockDown.cpp)
add_executable(CountryLockDown CountryLockDown.cpp)

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(CountyLockDown Eigen3::Eigen epigraph)
target_link_libraries(CountryLockDown Eigen3::Eigen epigraph)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(CountyLockDown OpenMP::OpenMP_CXX)
    target_link_libraries(CountryLockDown OpenMP::OpenMP_CXX)
endif ()
