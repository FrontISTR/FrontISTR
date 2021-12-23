cmake_minimum_required(VERSION 3.20)

find_package(Git QUIET)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD^1
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_HASH
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
string(TIMESTAMP BUILD_DATE)
file(READ "VERSION" fistr_ver)
string(REGEX MATCH "[0-9](.[0-9](.[0-9])?)?" fistr_ver ${fistr_ver})
message("===== FrontISTR v" ${fistr_ver} " (" ${GIT_HASH} ") =====")

#cmake_policy(SET CMP0048 NEW)
project(FrontISTR
  VERSION ${fistr_ver}
  LANGUAGES Fortran C CXX
  DESCRIPTION "Open-Source Large-Scale Parallel FEM Program for Nonlinear Structural Analysis:FrontISTR"
  HOMEPAGE_URL "https://www.frontistr.com/"
  )

find_package(MPI)
find_package(OpenMP)
configure_file(
  "${PROJECT_SOURCE_DIR}/FrontISTRConfig.h.in"
  "${PROJECT_BINARY_DIR}/FrontISTRConfig.h"
  )

add_subdirectory(hecmw1)
add_subdirectory(fistr1)

