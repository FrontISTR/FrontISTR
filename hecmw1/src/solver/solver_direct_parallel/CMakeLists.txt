###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/matrix2graph.c
  ${CMAKE_CURRENT_LIST_DIR}/matrix_repart.c
  ${CMAKE_CURRENT_LIST_DIR}/mem_util.c
  ${CMAKE_CURRENT_LIST_DIR}/separator_c2f_c.c
  ${CMAKE_CURRENT_LIST_DIR}/m_irjc_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_crs_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_child_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_matrix_partition_info.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_elap.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_direct_parallel.F90
)
