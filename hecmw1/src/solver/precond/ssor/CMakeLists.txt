###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_11.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_22.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_33.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_44.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_66.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_nn.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR.f90
)

target_include_directories(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
)

