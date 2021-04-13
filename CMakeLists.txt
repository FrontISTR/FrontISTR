###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

###############################################################
# CMakeLists.txt for FrontISTR
###############################################################
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
project(FrontISTR
	VERSION ${fistr_ver}
	DESCRIPTION "Open-Source Large-Scale Parallel FEM Program for Nonlinear Structural Analysis:FrontISTR"
	HOMEPAGE_URL "https://www.frontistr.com/"
	LANGUAGES Fortran C CXX
  )

include(FeatureSummary)
include(CMakeDependentOption)
include(FortranCInterface)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

###############################################################
# options
###############################################################
find_package(MPI COMPONENTS C CXX Fortran)
find_package(OpenMP COMPONENTS C CXX Fortran)
find_package(LAPACK)
find_package(METIS)
option(with_tools "compile tools" ON)
option(with_mpi "for parallel environment with MPI" ${MPI_FOUND})
option(with_openmp "for multi-(core|processor) environment with OpenMP" ${OPENMP_FOUND})
option(with_metis "use METIS" ${METIS_FOUND})
feature_summary(WHAT ALL)

configure_file(
  "${PROJECT_SOURCE_DIR}/FrontISTRConfig.h.in"
  "${PROJECT_BINARY_DIR}/FrontISTRConfig.h"
  )
add_subdirectory(hecmw1)
#add_subdirectory(fistr1)

