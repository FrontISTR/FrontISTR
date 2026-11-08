###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fstr_api_param.f90
)

target_include_directories(fistr PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
)
