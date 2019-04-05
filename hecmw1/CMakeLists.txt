###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

project(hecmw C CXX Fortran)

set(hecmw_SOURCES)
set(hecmw_INCLUDE_DIRS)
set(hecmw_DEFINITIONS)
set(hecmw_LIBRARIES)

#
# Inherit definitions from fistr options
#
set(hecmw_INCLUDE_DIRS ${FrontISTR_INCLUDE_DIRS})
set(hecmw_DEFINITIONS ${FrontISTR_DEFINITIONS})
set(hecmw_LIBRARIES ${FrontISTR_LIBRARIES})

#
# Build for libhecmw.a
#
include(src/CMakeLists.txt)

add_library(hecmw STATIC ${hecmw_SOURCES})
target_link_libraries(hecmw ${hecmw_LIBRARIES})
target_include_directories(hecmw PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(hecmw PUBLIC ${hecmw_DEFINITIONS})

#
# Build tools
#
if(WITH_TOOLS)
  add_subdirectory(tools)
endif()
