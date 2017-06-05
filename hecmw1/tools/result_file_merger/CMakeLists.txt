###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

list(APPEND rmerge_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/fstr_rmerge_util.c
  ${CMAKE_CURRENT_LIST_DIR}/fstr_rmerge.c
)

