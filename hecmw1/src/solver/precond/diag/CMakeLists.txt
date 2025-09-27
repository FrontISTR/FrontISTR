###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_11.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_22.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG.f90
)

target_include_directories(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
)

