cmake_minimum_required(VERSION 3.12)
project(Coursework1)

include_directories(../../Eigen)

set(CMAKE_CXX_STANDARD 11)

add_executable(Coursework1 main.cpp CubicSpline.cpp CubicSpline.hpp CubicSplinePeriodic.cpp CubicSplinePeriodic.hpp)