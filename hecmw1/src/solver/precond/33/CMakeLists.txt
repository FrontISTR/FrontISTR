###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_wrapper_33.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_BILU_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_DIAG_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SSOR_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_SAINV_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_RIF_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_ML_helper_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_ML_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_precond_33.f90
)

