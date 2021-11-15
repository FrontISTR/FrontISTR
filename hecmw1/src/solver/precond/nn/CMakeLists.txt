###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_nn_c.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_wrapper.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAINV_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_RIF_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_nn_f.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_ML_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_nn.f90
)

list(APPEND hecmw_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)

