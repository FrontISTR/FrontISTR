###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

if(NOT WITH_PARACON)
  list(APPEND fistr1_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/fistr_main.f90
  )
else()
  list(APPEND fistr1_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/fistr_main_contact.f90
  )
endif()

