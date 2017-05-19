###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_22.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_CG_22.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_BLCG_22.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_22.f90"
)

