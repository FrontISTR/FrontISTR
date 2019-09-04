###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/heat_init.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_io.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_get_amplitude.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_echo.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_bc_CFLUX.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_bc_DFLUX.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_bc_FILM.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_bc_FIXT.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_bc_RADIATE.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_boundary.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_capacity.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_mat_ass_conductivity.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_solve_main.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_solve_TRAN.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_solve_heat.f90
)
