###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_executable(neu2fstr)
target_sources(neu2fstr
  PRIVATE
    neu2fstr.cpp
)
target_include_directories(neu2fstr
  PRIVATE
  ${CMAKE_BINARY_DIR}/fistr1
  ${CMAKE_BINARY_DIR}/hecmw1
  )
target_link_libraries(neu2fstr
  PRIVATE
  hecmw
  fistr
  )

add_subdirectory(HECD)
add_subdirectory(NFD)
add_subdirectory(converter)
