###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

###############################################################
# CMakeLists.txt for FrontISTR
###############################################################
cmake_minimum_required(VERSION 2.8.11)

###############################################################
# Check building mode whether out-of-source builds
###############################################################
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "
*************** ERROR ***************
 Only support for out-of-source building.
 You may need to remove CMakeCache.txt file and CMakeFiles/ directory.
 Please make new directory such as build/ and run cmake in it again.
[example]
 rm CMakeCache.txt; rm -rf CMakeFiles/
 mkdir build/; cd build/
 cmake ..
 make")
endif()

project(FrontISTR Fortran C CXX)
###############################################################
# Includes standard macros
###############################################################
include(FeatureSummary)
include(CMakeDependentOption)
include(CMakeDetermineFortranCompiler)

###############################################################
# Add module path for FindPackage.cmake
###############################################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

###############################################################
# Initialize variable of definitions, link libraries, includes
###############################################################
set(FrontISTR_DEFINITIONS)
set(FrontISTR_LIBRARIES)
set(FrontISTR_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} ${PROJECT_BINARY_DIR})

###############################################################
# Build type (RELEASE or DEBUG)
###############################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Build type (RELEASE or DEBUG)" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS RELEASE DEBUG)
endif()
if(CMAKE_BUILD_TYPE MATCHES "DEBUG")
  list(APPEND FrontISTR_DEFINITIONS "DEBUG")
  # Enable debug options
  find_package(Backtrace)
  if(Backtrace_FOUND)
    list(APPEND FrontISTR_LIBRARIES ${Backtrace_LIBRARY})
    list(APPEND FrontISTR_INCLUDE_DIR ${Backtrace_INCLUDE_DIR})
    set(CMAKE_VERBOSE_MAKEFILE ON)
    # GNU Compiler
    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
      set(CMAKE_C_FLAGS_DEBUG
        "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra")
      set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
      set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -Wextra -fbounds-check -fbacktrace")
      if(${DEBUG_EXTRA})
        set(CMAKE_C_FLAGS_DEBUG
          "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address")
        set(CMAKE_CXX_FLAGS_DEBUG
          "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
        set(CMAKE_Fortran_FLAGS_DEBUG
          "${CMAKE_Fortran_FLAGS_DEBUG} -fsanitize=address")
      endif()
    # Intel Compiler
    elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
      set(CMAKE_C_FLAGS_DEBUG
        "${CMAKE_C_FLAGS_DEBUG} -Wall -Wbrief -Wcheck")
      set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wbrief -Wcheck")
      set(CMAKE_Fortran_FLAGS_DEBUG
        #"${CMAKE_Fortran_FLAGS_DEBUG} -warn all -check all -traceback")
        "${CMAKE_Fortran_FLAGS_DEBUG} -warn nounused -check uninit,bound -traceback")
    endif()
  endif()
elseif (CMAKE_BUILD_TYPE MATCHES "RELEASE")
  list(APPEND FrontISTR_DEFINITIONS "NDEBUG")
endif()

###############################################################
# Version number
###############################################################
#include(GetGitRevisionDescription)
#get_git_head_revision(GIT_REFSPEC GIT_SHA1)
find_package(Git QUIET)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD OUTPUT_VARIABLE GIT_SHA1 OUTPUT_STRIP_TRAILING_WHITESPACE)

file(READ "VERSION" ver)
string(REGEX MATCH "[0-9](.[0-9](.[0-9])?)?" vertxt ${ver})
string(REGEX MATCHALL "[0-9]" verlist "${vertxt}")
list(LENGTH verlist verlength)
set(VERSION_MAJOR 0)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
list(GET verlist 0 VERSION_MAJOR)
if(${verlength} GREATER 1)
  list(GET verlist 1 VERSION_MINOR)
endif()
if(${verlength} GREATER 2)
  list(GET verlist 2 VERSION_PATCH)
endif()
#execute_process(COMMAND date "+%Y-%m-%dT%H:%M:%S%z" OUTPUT_VARIABLE BUILD_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
string(TIMESTAMP BUILD_DATE "%Y-%m-%dT%H:%M:%S%z")
set(GIT_HASH \"${GIT_SHA1}\" CACHE STRING "GIT hash number of FrontISTR")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH BUILD_DATE GIT_HASH)

#############################################
# Architecture and compiler dependent options
#############################################
if(WIN32)
  set(WINDOWS TRUE)
endif()

###############################################################
# Options
###############################################################
# -DWINDOWS          : build on Windows
# -DWITH_TOOLS       : compile tools
# -DWITH_MPI         : for parallel environment with MPI
# -DWITH_OPENMP      : for multi-(core|processor) environment
# -DWITH_REFINER     : compile with REVOCAP_Refiner
# -DWITH_REVOCAP     : compile with REVOCAP_Coupler
# -DWITH_METIS       : compile with METIS graph partitioning package
#   -DMETIS_VER_4    : uses METIS-4
## -DWITH_PARMETIS    : compile with ParMETIS
##   -DPARMETIS_VER_3 : uses PARMETIS-3
# -DWITH_MKL         : compile with MKL PARDISO
# -DWITH_MUMPS       : compile with MUMPS direct solver
# -DWITH_LAPACK      : compile with LAPACK
# -DWITH_ML          : compile with Trilinos multigrid preconditioning package
# -DWITH_DOC         : generate API documents
# -DOLD_RES_FORMAT   : use old style result file( comment and global output disabled)

###################
# Find external packages
###################
option(WITH_MPI "for parallel environment with MPI" ON)
if(WITH_MPI)
  find_package(MPI)
  set(WITH_MPI "${MPI_FOUND}")
  #find_package(Parmetis)
  find_package(Scalapack)
endif()
find_package(OpenMP)
find_package(MKL)
find_package(LAPACK)
find_package(Metis)
find_package(Mumps)
find_package(Refiner)
find_package(Revocap)
find_package(Trilinos)
find_package(Doxygen)
option(WITH_TOOLS "compile tools" ON)
option(WITH_DOC "Generate API documents." ${WITH_DOC})
option(WINDOWS "build on windows" ${WINDOWS})
option(WITH_OPENMP "for multi-(core|processor) environment" ${OPENMP_FOUND})
option(WITH_MKL "compile with MKL PARDISO" ${MKL_FOUND})
option(WITH_LAPACK "for estimating number of condition" ${LAPACK_FOUND})
option(WITH_METIS "compile with METIS" ${METIS_FOUND})
option(METIS_VER_4 "uses METIS-4" OFF)
option(WITH_MUMPS "compile with MUMPS direct solver" ${MUMPS_FOUND})
#option(WITH_PARMETIS "compile with ParMETIS" ${PARMETIS_FOUND})
#option(PARMETIS_VER_3 "uses METIS-3" OFF)
option(WITH_REFINER "compile with REVOCAP_Refiner" ${REFINER_FOUND})
option(WITH_REVOCAP "compile with REVOCAP_Coupler" ${REVOCAP_FOUND})
option(WITH_ML "compile with Trilinos ML" ${Trilinos_FOUND})
option(OLD_RES_FORMAT "use old style result file" ${OLD_RES_FORMAT})
option(DEBUG_EXTRA "enable memory leak check (only for gcc4.9 <)" OFF)
feature_summary(WHAT ALL)

# Provisional solution FOR GCC Fortran 10 compiling error
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()

###################
# -DWITH_OPENMP
###################
if(WITH_OPENMP AND OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
else()
  set(WITH_OPENMP OFF)
endif()

###################
# -DWITH_MKL
###################
if(WITH_MKL)
  set(BLA_VENDOR "Intel10_64lp")
  list(APPEND FrontISTR_INCLUDE_DIRS ${MKL_INCLUDE_PATH})
  list(APPEND FrontISTR_LIBRARIES ${MKL_LIBRARIES})
  list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_MKL")
  set(WITH_LAPACK ON)
else()
  set(WITH_MKL OFF)
endif()

###################
# -DWITH_REFINER
###################
if(WITH_REFINER AND REFINER_FOUND)
  list(APPEND FrontISTR_INCLUDE_DIRS ${REFINER_INCLUDE_PATH})
  list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_REFINER")
  list(APPEND FrontISTR_LIBRARIES ${REFINER_LIBRARIES})
else()
  set(WITH_REFINER OFF)
  set(HECMW_WITH_REFINER OFF)
endif()

###################
# -DWITH_REVOCAP
###################
if(WITH_REVOCAP AND REVOCAP_FOUND)
  list(APPEND FrontISTR_INCLUDE_DIRS ${REVOCAP_INCLUDE_PATH})
else()
  set(WITH_REVOCAP OFF)
endif()

###################
# -DWITH_ML
###################
if(WITH_ML)
  if(Trilinos_FOUND)
    if(Trilinos_VERSION VERSION_LESS 14.0.0)
      list(FIND Trilinos_TPL_LIST MPI trilinos_mpi_index)
      if((trilinos_mpi_index EQUAL -1 AND WITH_MPI) OR (trilinos_mpi_index GREATER -1 AND NOT WITH_MPI))
        message(WARNING "Trilinos found but turned off due to incompatible MPI setting")
        set(WITH_ML OFF)
      else()
        list(APPEND FrontISTR_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
        list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_ML")
        list(APPEND FrontISTR_LIBRARIES ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
      endif()
    else()
      if((NOT ML_ENABLE_MPI AND WITH_MPI) OR (ML_ENABLE_MPI AND NOT WITH_MPI))
        message(WARNING "Trilinos found but turned off due to incompatible MPI setting")
        set(WITH_ML OFF)
      else()
        list(APPEND FrontISTR_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS})
        list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_ML")
        list(APPEND FrontISTR_LIBRARIES ${Trilinos_LIBRARIES})
      endif()
    endif()
  else()
    message(WARNING "Trilinos not found")
    set(WITH_ML OFF)
  endif()
endif()
mark_as_advanced(Trilinos_DIR)

###################
# -DWITH_MUMPS
###################
if(WITH_MUMPS)
  if(MUMPS_FOUND)
    if(WITH_MPI)
      if((SCALAPACK_FOUND AND LAPACK_FOUND) OR SCALAPACK_MKL)
        list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_MUMPS")
        list(APPEND FrontISTR_INCLUDE_DIRS ${MUMPS_INCLUDE_PATH})
        list(APPEND FrontISTR_LIBRARIES ${MUMPS_LIBRARIES} ${SCALAPACK_LIBRARIES})
        mark_as_advanced(scalapack_DIR)
      else()
        message(WARNING "MUMPS could not be enabled. Please also install scalapack and lapack.")
        set(WITH_MUMPS OFF)
      endif()
    else()
      list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_MUMPS")
      list(APPEND FrontISTR_INCLUDE_DIRS ${MUMPS_INCLUDE_PATH})
      list(APPEND FrontISTR_LIBRARIES ${MUMPS_LIBRARIES})
    endif()
  else()
    message(WARNING "MUMPS not found")
    set(WITH_MUMPS OFF)
  endif()
endif()

###################
# -DWITH_PARMETIS
###################
#if(WITH_PARMETIS AND PARMETIS_FOUND)
#  list(APPEND FrontISTR_INCLUDE_DIRS ${PARMETIS_INCLUDE_PATH})
#  list(APPEND FrontISTR_LIBRARIES ${PARMETIS_LIBRARIES})
#  if(NOT PARMETIS_VER_3)
#  else()
#  endif()
#endif()

###################
# -DWITH_METIS
# -DMETIS_VER_4
###################
if(WITH_METIS AND METIS_FOUND)
  list(APPEND FrontISTR_INCLUDE_DIRS ${METIS_INCLUDE_PATH})
  list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_METIS")
  list(APPEND FrontISTR_DEFINITIONS "HECMW_PART_WITH_METIS")
  list(APPEND FrontISTR_LIBRARIES ${METIS_LIBRARIES})
  if(NOT METIS_VER_4)
    list(APPEND FrontISTR_DEFINITIONS "HECMW_METIS_VER=5")
  else()
    list(APPEND FrontISTR_DEFINITIONS "HECMW_METIS_VER=4")
  endif()
endif()

###################
# -DWITH_LAPACK
###################
if(WITH_LAPACK)
  if(LAPACK_FOUND AND BLAS_FOUND)
    list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_LAPACK")
    list(APPEND FrontISTR_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  else()
    set(WITH_LAPACK OFF)
  endif()
endif()


###################
# -DOLD_RES_FORMAT
###################
if(OLD_RES_FORMAT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DOLD_RES_FORMAT")
endif()

###################
# -DWITH_MPI
###################
if(WITH_MPI AND MPI_Fortran_FOUND AND MPI_CXX_FOUND AND MPI_C_FOUND)
  list(APPEND FrontISTR_INCLUDE_DIRS
    ${MPI_C_INCLUDE_PATH}
    ${MPI_CXX_INCLUDE_PATH}
    ${MPI_Fortran_INCLUDE_PATH})

  list(APPEND FrontISTR_LIBRARIES
    ${MPI_C_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${MPI_Fortran_LIBRARIES})

else()
  set(WITH_MPI OFF)
  set(HECMW_SERIAL ON)
  list(APPEND FrontISTR_DEFINITIONS "HECMW_SERIAL")
endif()
mark_as_advanced(MPI_EXTRA_LIBRARY MPI_LIBRARY)

#############################################
# Architecture and compiler dependent options
#############################################
###################
# -DWINDOWS
###################
if(WINDOWS)
  list(APPEND FrontISTR_DEFINITIONS "_WINDOWS")
  list(APPEND FrontISTR_LIBRARIES "ws2_32")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()


###################
# -DWITH_DOC
###################
if(WITH_DOC)
  if(DOXYGEN_FOUND)
    message(STATUS "Generate API documents. Please type \"make doc\"")
    file(COPY ${CMAKE_SOURCE_DIR}/doc/fistr.css DESTINATION ${CMAKE_BINARY_DIR}/doc/)
    file(COPY ${CMAKE_SOURCE_DIR}/doc/fistr_logo.png DESTINATION ${CMAKE_BINARY_DIR}/doc/)

    set(doxy_main_page "${CMAKE_SOURCE_DIR}/README.md")
    set(doxyfile_in ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in)
    set(doxyfile ${CMAKE_BINARY_DIR}/doc/Doxyfile)
    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Generate API documents."
      VERBATIM
    )
  else()
    message(WARNING "doxygen not found.")
  endif()
endif()

###############################################################################
# Print STATUS
###############################################################################
message(STATUS
  "Build FrontISTR for ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION}")
message(STATUS "BUILD Type\t\t: ${CMAKE_BUILD_TYPE}")
message(STATUS "C compiler\t\t: ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler\t\t: ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran compiler\t: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Installation path\t: ${CMAKE_INSTALL_PREFIX}")

if(CMAKE_BUILD_TYPE MATCHES "DEBUG")
  message(STATUS "Defines\t\t: ${FrontISTR_DEFINITIONS}")
  message(STATUS "Include paths\t: ${FrontISTR_INCLUDE_DIRS}")
  message(STATUS "Library search path\t: ${CMAKE_PREFIX_PATH}")
  message(STATUS "Linking libraries\t: ${FrontISTR_LIBRARIES}")
endif()

###############################################################################
# Make configuration header
###############################################################################
configure_file(
  "${PROJECT_SOURCE_DIR}/FrontISTRConfig.h.in"
  "${PROJECT_BINARY_DIR}/FrontISTRConfig.h"
)

###############################################################################
# Build hecmw1, fistr1 and fistr
###############################################################################
# Build targets
# hecmw (libhecmw.a)
# fistr (libfistr.a)
# fistr1 (fistr1)
#
# Optional : rmerge, hecmw_part, hec2rcap, rconv, hecmw_vis1, neu2fstr
#
add_subdirectory(hecmw1)
add_subdirectory(fistr1)

###############################################################################
# Tests
###############################################################################
enable_testing()
add_subdirectory(tests)

###############################################################################
# Package generator (deb, rpm, exe etc.)
###############################################################################
include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Large-scale structural analysis program with FEM ")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VENDOR "FrontISTR Commons")
set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
set(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
set(CPACK_PACKAGE_CONTACT "${CPACK_PACKAGE_VENDOR}")
if(WIN32 AND NOT UNIX)
  # Windows
else()
  # Unix
  set(CPACK_BINARY_DEB ON)
  set(CPACK_BINARY_RPM ON)
  set(CPACK_GENERATOR "TBZ2;DEB;RPM")
  set(CPACK_STRIP_FILES ON)
endif()
include(CPack)



