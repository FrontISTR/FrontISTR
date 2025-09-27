###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/uelastic.f90
  ${CMAKE_CURRENT_LIST_DIR}/uhardening.f90
  ${CMAKE_CURRENT_LIST_DIR}/uload.f90
  ${CMAKE_CURRENT_LIST_DIR}/umat.f90
  ${CMAKE_CURRENT_LIST_DIR}/uyield.f90
)
