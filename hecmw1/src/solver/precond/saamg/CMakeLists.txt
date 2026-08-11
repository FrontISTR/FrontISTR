###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_util.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_param.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_comm.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_matrix.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_smoother.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_aggregate.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_prolongation.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_coarse.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_coarse_mumps.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_core.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_adapt.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG_verify.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAAMG.f90
)

target_include_directories(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
)
