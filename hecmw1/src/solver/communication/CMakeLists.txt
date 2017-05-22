###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_11.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_11i.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_22.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_22i.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_33.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_33i.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_44.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_44i.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_66.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_66i.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_mm.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_SR_mmi.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_comm_f.F90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_misc.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_comm_group.f90"
)

