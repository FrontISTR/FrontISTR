###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

include(${CMAKE_CURRENT_LIST_DIR}/communication/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/las/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/matrix/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/init/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/sparse_matrix/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/mumps/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/clustermkl/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/precond/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/iterative/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/solver_direct/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/solver_direct_parallel/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/solver_direct_lag/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/mpc/CMakeLists.txt)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver.f90
)
