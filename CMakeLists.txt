cmake_minimum_required(VERSION 3.25)
project(CompTask_V)

set(CMAKE_CXX_STANDARD 23)

include_directories(C:/MyFiles/ExtLib\ C++/boost_1_82_0)

add_executable(CompTask_V main.cpp RPN.cpp Integrator.cpp Integrator.cpp Task5.1.cpp Task5.2.cpp Task5.3.cpp InterpolateQF.cpp Matrix.cpp test.cpp invert_matrix.cpp)