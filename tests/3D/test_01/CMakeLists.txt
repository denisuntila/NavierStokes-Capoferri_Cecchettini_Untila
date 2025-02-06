cmake_minimum_required(VERSION 3.12.0)
project(test_01 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(main src/test_01.cpp ../../../src/NavierStokes.cpp)
deal_ii_setup_target(main)
