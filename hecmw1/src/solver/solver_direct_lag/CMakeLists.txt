###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/my_hecmw_util_f_lag.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_irjc_matrix_lag.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_child_matrix_lag.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_crs_matrix_lag.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_direct_serial_lag.f90
)
