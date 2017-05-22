###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

project(fistr C CXX Fortran)

set(fistr_SOURCES)
set(fistr_INCLUDE_DIRS)
set(fistr_DEFINITIONS)
set(fistr_LIBRARIES)

#
# Inherit definitions from fistr options
#
set(fistr_INCLUDE_DIRS ${FrontISTR_INCLUDE_DIRS})
set(fistr_DEFINITIONS ${FrontISTR_DEFINITIONS})
set(fistr_LIBRARIES ${FrontISTR_LIBRARIES})

#
# Build for libfistr.a
#
include(src/CMakeLists.txt)
add_library(fistr STATIC ${fistr_SOURCES})
target_include_directories(fistr PUBLIC ${fistr_INCLUDE_DIRS})
target_compile_definitions(fistr PUBLIC ${fistr_DEFINITIONS})

#
# Build for fistr1 (executable)
#

#include_directories(
#  ${hecmw_SOURCE_DIR}/src/common
#  ${hecmw_SOURCE_DIR}/src/hecmw
#  ${hecmw_SOURCE_DIR}/src/visualizer
#  ${hecmw_BINARY_DIR}/src
#  ${Fistr_SOURCE_DIR}/src/common
#  ${Fistr_SOURCE_DIR}/src/lib
#  ${Fistr_BINARY_DIR}/src
#)

#add_subdirectory(src)

#if(WITH_TOOLS)
#  add_subdirectory(tools)
#endif()
