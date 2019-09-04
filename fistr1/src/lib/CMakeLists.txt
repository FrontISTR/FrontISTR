###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

include(${CMAKE_CURRENT_LIST_DIR}/contact/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/element/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/physics/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/user/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/utilities/CMakeLists.txt)

list(APPEND fistr_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/GaussM.f90
  ${CMAKE_CURRENT_LIST_DIR}/dynamic_mass.f90
  ${CMAKE_CURRENT_LIST_DIR}/eigen_LIB.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB_CAPACITY.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB_CONDUCTIVITY.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB_DFLUX.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB_FILM.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB_NEUTRAL.f90
  ${CMAKE_CURRENT_LIST_DIR}/heat_LIB_RADIATE.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_common_struct.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_fstr.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_timepoint.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_out.f90
  ${CMAKE_CURRENT_LIST_DIR}/m_step.f90
  ${CMAKE_CURRENT_LIST_DIR}/precheck_LIB_2d.f90
  ${CMAKE_CURRENT_LIST_DIR}/precheck_LIB_3d.f90
  ${CMAKE_CURRENT_LIST_DIR}/precheck_LIB_shell.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_LINEQ.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_1d.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_2d.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_3d.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_3dIC.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_3d_vp.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_C3D8.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_Fbar.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_beam.f90
  ${CMAKE_CURRENT_LIST_DIR}/static_LIB_shell.f90
)
