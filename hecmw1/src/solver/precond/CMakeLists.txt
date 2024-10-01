###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

include(${CMAKE_CURRENT_LIST_DIR}/bilu/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/diag/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/ml/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/rif/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/sainv/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/ssor/CMakeLists.txt)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond.f90
)
