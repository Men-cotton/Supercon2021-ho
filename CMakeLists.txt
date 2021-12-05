cmake_minimum_required(VERSION 3.20)
project(Supercon2021_ho)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(Supercon2021_ho
        replica.cpp
        sc21.h)
