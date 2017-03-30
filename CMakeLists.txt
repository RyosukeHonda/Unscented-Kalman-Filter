cmake_minimum_required(VERSION 3.6)
project(Unscented_Kalman_Filter)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/main.cpp
        src/tools.cpp
        src/ukf.cpp)

add_executable(Unscented_Kalman_Filter ${SOURCE_FILES})