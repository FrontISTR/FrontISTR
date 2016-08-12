cmake_minimum_required(VERSION 2.8.11)

project(Fistr C CXX Fortran)

include_directories(
  ${Hecmw_SOURCE_DIR}/src/common
  ${Hecmw_SOURCE_DIR}/src/hecmw
  ${Hecmw_SOURCE_DIR}/src/visualizer
  ${Hecmw_BINARY_DIR}/src
  ${Fistr_SOURCE_DIR}/src/common
  ${Fistr_SOURCE_DIR}/src/lib
  ${Fistr_BINARY_DIR}/src
)

add_subdirectory(src)

if(WITH_TOOLS)
  add_subdirectory(tools)
endif()
