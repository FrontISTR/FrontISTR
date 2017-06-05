###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

list(APPEND neu2fstr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/CNFData.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_100.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_402.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_403.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_404.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_405.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_408.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_506.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_507.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFDB_601.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CNFMessage.cpp
)

list(APPEND neu2fstr_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)

