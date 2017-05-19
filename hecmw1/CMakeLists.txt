###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

project(hecmw C CXX Fortran)

set(hecmw_SOURCES)
set(hecmw_INCLUDE_DIR)
set(hecmw_DEFINITIONS)

#
# Inherit definitions from fistr options
#

#
# Build for libhecmw.a
#
include(src/CMakeLists.txt)

add_library(hecmw STATIC ${hecmw_SOURCES})
target_link_libraries(hecmw ${MPI_CXX_LIBRARIES})
target_include_directories(hecmw PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(hecmw PUBLIC ${hecmw_DEFINITIONS})

#
# Build tools
#
if(WITH_TOOLS)
  add_subdirectory(tools)
endif()
