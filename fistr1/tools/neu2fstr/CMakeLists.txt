###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

include(${CMAKE_CURRENT_LIST_DIR}/HECD/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/NFD/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/converter/CMakeLists.txt)

list(APPEND neu2fstr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/neu2fstr.cpp
)

