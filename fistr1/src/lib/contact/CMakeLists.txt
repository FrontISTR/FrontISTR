###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/bucket_search.f90
  ${CMAKE_CURRENT_LIST_DIR}/contact_lib.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_comm.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_def.F90
  ${CMAKE_CURRENT_LIST_DIR}/surf_ele.f90
)
