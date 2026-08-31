###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_api_common.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_api_local_mesh.f90
)

target_include_directories(hecmw PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
)

