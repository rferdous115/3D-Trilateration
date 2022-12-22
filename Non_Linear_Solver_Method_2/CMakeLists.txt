cmake_minimum_required(VERSION 3.23)
project(Non_Linear_Solver_Method_2)

set(EIGEN_DIR "/Users/rheinferdous/Downloads/eigen-3.4.0")
if(NOT EXISTS ${EIGEN_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
include_directories(${EIGEN_DIR})

set(CMAKE_CXX_STANDARD 23)

add_executable(Non_Linear_Solver_Method_2 main.cpp)
