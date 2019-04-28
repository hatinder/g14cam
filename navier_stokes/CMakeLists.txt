cmake_minimum_required(VERSION 3.13)
project(navier_stokes)

set(CMAKE_CXX_STANDARD 11)

add_executable(navier_stokes main.cpp ns2de.cpp ns2de.hpp)