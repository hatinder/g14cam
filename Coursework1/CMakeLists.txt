cmake_minimum_required(VERSION 3.12)
project(Coursework1)

include_directories(../../Eigen)

set(CMAKE_CXX_STANDARD 11)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/spline_natural1.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/spline_natural2.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/spline_periodic.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(Coursework1 main.cpp CubicSpline.cpp CubicSpline.hpp CubicSplinePeriodic.cpp CubicSplinePeriodic.hpp RootFinding.cpp RootFinding.hpp)