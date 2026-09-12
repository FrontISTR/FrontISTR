###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mf_graph.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mf_symbolic.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mf_kernel.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mf_dist.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mf_numeric.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_direct_mf.f90
)
