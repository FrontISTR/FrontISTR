###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/fstr_Ctrl_TimeInc.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_Cutback.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_Spring.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_Restart.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_NodalStress.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_StiffMatrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_Update.f90
  ${CMAKE_CURRENT_LIST_DIR}/readtemp.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_echo.f90
  ${CMAKE_CURRENT_LIST_DIR}/make_result.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ass_load.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_output.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_mat_con_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_Residual.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_AddContactStiff.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_AddBC.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_mat_resid_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/sparse_matrix_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_MUMPS_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/set_arrays_DirectSolver.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_mkl.F90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_direct_serial_lag.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_iter_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_solve_NonLinear.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_solve_NLGEOM.f90
)
