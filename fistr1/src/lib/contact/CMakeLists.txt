###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/bucket_search.f90
  ${CMAKE_CURRENT_LIST_DIR}/contact_lib.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_param.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_def.F90
  ${CMAKE_CURRENT_LIST_DIR}/surf_ele.f90
)
