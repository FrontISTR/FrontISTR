###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_lanczos_util.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_setMASS.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_tridiag.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_output.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_EIG_lanczos.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_solve_eigen.f90
)
