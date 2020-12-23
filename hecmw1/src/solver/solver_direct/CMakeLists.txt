###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ordering_qmd.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ordering_metis.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ordering.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_direct.f90
)
