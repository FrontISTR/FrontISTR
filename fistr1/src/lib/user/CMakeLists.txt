###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/uelastic.f90
  ${CMAKE_CURRENT_LIST_DIR}/uhardening.f90
  ${CMAKE_CURRENT_LIST_DIR}/uload.f90
  ${CMAKE_CURRENT_LIST_DIR}/umat.f90
  ${CMAKE_CURRENT_LIST_DIR}/uyield.f90
)
