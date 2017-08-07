###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_lczeigenm.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_getamat.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_lanczos_util.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_matmult.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_mgs1.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_setMASS.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_tridiag.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_output.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_lanczos.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_solve_eigen.f90
)
