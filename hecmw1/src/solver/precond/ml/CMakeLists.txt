###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_33_c.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_nn_c.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_wrapper.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_33_f.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_nn_f.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_ML_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_ML_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_ML.f90
)

list(APPEND hecmw_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)

