###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/sparse_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/sparse_matrix_hec.f90
)
