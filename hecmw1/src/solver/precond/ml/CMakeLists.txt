###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
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

target_include_directories(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
)

