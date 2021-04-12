###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_executable(rconv)

target_sources(rconv
  PRIVATE
  hecmw_res_type_conv.c
)
target_include_directories(rconv
  PRIVATE
  ../../src/common
  )
target_link_libraries(rconv
  PRIVATE
  hecmw
  )
install(TARGETS rconv
  DESTINATION bin
  )
