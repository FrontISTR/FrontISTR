###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/element.f90
  ${CMAKE_CURRENT_LIST_DIR}/hex20n.f90
  ${CMAKE_CURRENT_LIST_DIR}/hex8n.f90
  ${CMAKE_CURRENT_LIST_DIR}/line2n.f90
  ${CMAKE_CURRENT_LIST_DIR}/line3n.f90
  ${CMAKE_CURRENT_LIST_DIR}/prism15n.f90
  ${CMAKE_CURRENT_LIST_DIR}/prism6n.f90
  ${CMAKE_CURRENT_LIST_DIR}/quad4n.f90
  ${CMAKE_CURRENT_LIST_DIR}/quad8n.f90
  ${CMAKE_CURRENT_LIST_DIR}/quad9n.f90
  ${CMAKE_CURRENT_LIST_DIR}/quadrature.f90
  ${CMAKE_CURRENT_LIST_DIR}/tet10n.f90
  ${CMAKE_CURRENT_LIST_DIR}/tet4n.f90
  ${CMAKE_CURRENT_LIST_DIR}/tri3n.f90
  ${CMAKE_CURRENT_LIST_DIR}/tri6n.f90
)
