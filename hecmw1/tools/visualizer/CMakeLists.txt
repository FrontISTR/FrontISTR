###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_executable(hecmw_vis1)
target_sources(hecmw_vis1
  PRIVATE
  visualizer.c
)
target_include_directories(hecmw_vis1
  PRIVATE
  ../../src/common
  ../../src/visualizer
  )
target_link_libraries(hecmw_vis1
  PRIVATE
  hecmw
  )
install(TARGETS hecmw_vis1
  DESTINATION bin
  )
