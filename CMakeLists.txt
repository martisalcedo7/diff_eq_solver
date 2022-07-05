cmake_minimum_required(VERSION 3.10.2)

project(myApp
    VERSION 1.0
    DESCRIPTION "Runge-Kutta solver"
    LANGUAGES CXX)


add_executable(myApp src/main.cpp)

target_include_directories(myApp  
    PRIVATE ${CMAKE_SOURCE_DIR}/src 
    PRIVATE /usr/include/python3.6)

target_link_libraries(myApp
PRIVATE python3.6m)

target_compile_features(myApp PRIVATE cxx_std_17)

target_compile_options(myApp PRIVATE -Wall -Wextra -Wpedantic)

