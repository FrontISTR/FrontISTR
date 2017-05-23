###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ClusterMKL_wrapper.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_direct_ClusterMKL.f90
)
