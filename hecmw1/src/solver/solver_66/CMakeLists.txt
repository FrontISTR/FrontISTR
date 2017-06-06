###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_CG_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_66.f90
)

