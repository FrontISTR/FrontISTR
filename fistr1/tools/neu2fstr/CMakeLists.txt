###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

include(${CMAKE_CURRENT_LIST_DIR}/HECD/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/NFD/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/converter/CMakeLists.txt)

list(APPEND neu2fstr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/neu2fstr.cpp
)

