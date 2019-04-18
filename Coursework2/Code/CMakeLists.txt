cmake_minimum_required(VERSION 3.13)
project(Code)

set(CMAKE_CXX_STANDARD 11)


add_executable(CW main.cpp GaussQuadrature.cpp GaussQuadrature.hpp IGaussQuadrature.hpp IUtility.hpp Utility.cpp Utility.hpp helper.hpp)