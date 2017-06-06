###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_CG_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_BiCGSTAB_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_44.f90
)

