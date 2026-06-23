###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_MKL_wrapper.F90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_direct_MKL.f90
)
