###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_44.f90
)

