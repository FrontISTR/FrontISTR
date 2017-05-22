###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_wrapper.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_estimate_condition.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_tuning_fx.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAINV_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_RIF_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_mat_id.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_ML_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_CG_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_BiCGSTAB_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GMRES_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GPBiCG_33.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_33.f90"
)

