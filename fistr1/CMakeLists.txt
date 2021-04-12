###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
add_executable(fistr1
  src/main/main.c
  src/main/fistr_main.f90
  )
target_link_libraries(fistr1
  PRIVATE
  hecmw
  fistr
  )
target_include_directories(fistr1
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/fistr1
  ${CMAKE_BINARY_DIR}/hecmw1
  )
set_target_properties(fistr1
  PROPERTIES
  LINKER_LANGUAGE C
  )

add_library(fistr)

add_subdirectory(src)

target_include_directories(fistr
  PRIVATE
    src/common
    ../hecmw1/src/common
    ${CMAKE_BINARY_DIR}/hecmw1
  )


#add_subdirectory(tools)
