###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE 
  ${CMAKE_CURRENT_LIST_DIR}/hecmw2fstr_connect_conv.c
  ${CMAKE_CURRENT_LIST_DIR}/fstr_sort_index.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw2fstr_mesh_conv.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_setup_util.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_modifier.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_common.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_material.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_static.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_heat.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_eigen.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_ctrl_dynamic.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_get_prop.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_setup.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_debug_dump.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_precheck.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_rcap_io.F90
)

target_include_directories(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
)
