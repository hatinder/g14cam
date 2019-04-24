cmake_minimum_required(VERSION 3.13)
project(Code)
include_directories(../Eigen)

set(CMAKE_CXX_STANDARD 11)


add_executable(CW main.cpp GaussQuadrature.cpp GaussQuadrature.hpp IGaussQuadrature.hpp IUtility.hpp Utility.cpp Utility.hpp helper.hpp Problem2.cpp Problem2.hpp ODE.cpp ODE.hpp Problem3.cpp Problem3.hpp)