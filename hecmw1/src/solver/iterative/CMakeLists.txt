###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_CG.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_BiCGSTAB.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GMRES.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GPBiCG.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_Iterative.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_Iterative_SXAT.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_Iterative_SXAT_entry.c
)

