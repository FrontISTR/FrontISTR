###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/table_dyn.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_mat_ass_load.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_mat_ass_bc.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_mat_ass_bc_vl.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_mat_ass_bc_ac.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_mat_ass_couple.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_var_init.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_make_result.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_output.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_dynamic_nlexplicit.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_dynamic_nlimplicit.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_solve_dynamic.f90
)
