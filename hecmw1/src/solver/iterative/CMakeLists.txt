###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_CG.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_BiCGSTAB.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GMRES.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GPBiCG.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_Iterative.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GMRESR.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_GMRESREN.f90
)

