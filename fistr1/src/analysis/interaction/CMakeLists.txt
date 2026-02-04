###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact_lib.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_mat_con_contact.f90
  ${CMAKE_CURRENT_LIST_DIR}/fstr_AddContactStiff.f90
)
