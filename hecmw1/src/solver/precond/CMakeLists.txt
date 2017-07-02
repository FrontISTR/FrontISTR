###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

include(${CMAKE_CURRENT_LIST_DIR}/11/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/22/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/33/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/44/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/66/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/nn/CMakeLists.txt)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond.f90
)
