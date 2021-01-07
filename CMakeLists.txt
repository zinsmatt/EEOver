cmake_minimum_required(VERSION 3.4...3.18)
project(ellipses_iou)

set(CMAKE_CXX_STANDARD 17)

include_directories(./)
link_directories(./)

add_subdirectory(pybind11)
pybind11_add_module(${PROJECT_NAME} SHARED ee.cpp solvers.cpp zsolve_quartic.cpp Roots3And4.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE
gsl 
gslcblas
)




