###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_66.f90
)

