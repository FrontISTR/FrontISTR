###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_smoothing.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_geom.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_interference.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_elem_common.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_elem_slag.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_elem_alag.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_element.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_search.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_damping.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_assembly.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_mpc.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_output.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_mat_con_contact.f90
)
