###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

include(${CMAKE_CURRENT_LIST_DIR}/analysis/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/common/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/lib/CMakeLists.txt)
#include(${CMAKE_CURRENT_LIST_DIR}/main/CMakeLists.txt)
#if(NOT WITH_PARACON)
#  set(Fistr_MAIN_SRCS main/fistr_main.f90)
#else()
#  set(Fistr_MAIN_SRCS main/fistr_main_contact.f90)
#endif()

#install(TARGETS fistr1 DESTINATION bin CONFIGURATIONS ${CMAKE_BUILD_TYPE})
