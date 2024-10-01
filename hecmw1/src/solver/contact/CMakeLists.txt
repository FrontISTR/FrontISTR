###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_resid_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/sparse_matrix_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_MUMPS_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/set_arrays_DirectSolver.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_mkl.F90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_direct_serial_lag.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_contact_elim.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_contact.f90
)

