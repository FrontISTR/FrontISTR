###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_i.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_comm_f.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_comm_contact_f.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_misc.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_comm_group.f90
)

