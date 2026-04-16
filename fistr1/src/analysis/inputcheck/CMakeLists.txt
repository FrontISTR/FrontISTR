###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fstr_precheck.f90
  ${CMAKE_CURRENT_LIST_DIR}/precheck_LIB_2d.f90
  ${CMAKE_CURRENT_LIST_DIR}/precheck_LIB_3d.f90
  ${CMAKE_CURRENT_LIST_DIR}/precheck_LIB_shell.f90
)
