###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_allocate.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_misc.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_ordering_CM.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_ordering_MC.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_reorder.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_con.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_ass.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_dump.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_local_matrix.f90
)

