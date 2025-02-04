cmake_minimum_required(VERSION 3.12.0)
project(NavierStokes LANGUAGES CXX C)

include(./common/cmake-common.cmake)

add_executable(NavierStokes src/main.cpp src/NavierStokes3D.cpp)
deal_ii_setup_target(NavierStokes)

add_executable(postprocess src/postprocess.cpp src/NavierStokes3D.cpp)
deal_ii_setup_target(postprocess)





