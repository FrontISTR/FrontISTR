###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_11.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_22.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_tuning_fx.f90
)

