cmake_minimum_required(VERSION 3.7)
project(RayTracer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -framework GLUT -framework OpenGL")

set(SOURCE_FILES main.cpp)
add_executable(RayTracer ${SOURCE_FILES})