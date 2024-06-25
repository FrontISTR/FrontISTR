###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU.f90
)

list(APPEND hecmw_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)

