###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_allocate.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_array_util.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_estimate_condition.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_local_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_ass.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_con.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_id.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_contact_lagrange.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_dump.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_misc.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_ordering_CM.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_ordering_MC.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_matrix_reorder.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_pair_array.f90
)

