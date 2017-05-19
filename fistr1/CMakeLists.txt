###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

project(Fistr C CXX Fortran)

include_directories(
  ${hecmw_SOURCE_DIR}/src/common
  ${hecmw_SOURCE_DIR}/src/hecmw
  ${hecmw_SOURCE_DIR}/src/visualizer
  ${hecmw_BINARY_DIR}/src
  ${Fistr_SOURCE_DIR}/src/common
  ${Fistr_SOURCE_DIR}/src/lib
  ${Fistr_BINARY_DIR}/src
)

add_subdirectory(src)

if(WITH_TOOLS)
  add_subdirectory(tools)
endif()
