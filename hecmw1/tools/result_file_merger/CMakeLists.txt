###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
add_executable(rmerge)

target_sources(rmerge
  PRIVATE
    fstr_rmerge_util.c
    fstr_rmerge.c
  PUBLIC
    fstr_rmerge_util.h
  )

target_include_directories(rmerge
  PUBLIC
  ../../src/common
  )

target_link_libraries(rmerge
  PUBLIC
  hecmw
  )

install(TARGETS rmerge
  DESTINATION bin
  )
