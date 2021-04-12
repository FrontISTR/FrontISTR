###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
add_executable(hec2rcap)
target_sources(hec2rcap
  PRIVATE
  hec2rcap.c
  )
target_link_libraries(hec2rcap
  PRIVATE
  hecmw
  )
target_include_directories(hec2rcap
  PUBLIC
  ../../src/hecmw
  )
install(TARGETS hec2rcap
  DESTINATION bin
  )
