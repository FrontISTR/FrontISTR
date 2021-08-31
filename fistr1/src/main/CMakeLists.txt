###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND fistr1_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/fistr_main.f90
  ${CMAKE_CURRENT_LIST_DIR}/vhcall_wrapper.f90
  ${CMAKE_CURRENT_LIST_DIR}/main.c
)
