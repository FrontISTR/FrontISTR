###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

project(Hecmw C CXX Fortran)

include_directories(
  ${Hecmw_SOURCE_DIR}/src/common
  ${Hecmw_SOURCE_DIR}/src/hecmw
  ${Hecmw_SOURCE_DIR}/src/visualizer
  ${Hecmw_BINARY_DIR}
)

add_subdirectory(src)

if(WITH_TOOLS)
  add_subdirectory(tools)
endif()
