###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

###############################################################################
# CMakeLists.txt for FrontISTR
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

project(FrontISTR Fortran C CXX)

###############################################################################
# Includes standard macros
###############################################################################
include(FeatureSummary)
include(CMakeDependentOption)

###############################################################################
# Add module path for FindPackage.cmake
###############################################################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

###############################################################################
# Add Fortran module path
###############################################################################
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

###############################################################################
# Initialize variable of definitions, link libraries, includes
###############################################################################
set(FrontISTR_DEFINITIONS)
set(FrontISTR_LIBRARIES)
set(FrontISTR_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY})

###############################################################################
# Options
###############################################################################
# -DWINDOWS          : build on Windows
# -DWITH_TOOLS       : compile tools
# -DWITH_MPI         : for parallel environment with MPI
# -DWITH_OPENMP      : for multi-(core|processor) environment
# -DWITH_REFINER     : compile with REVOCAP_Refiner
# -DWITH_REVOCAP     : compile with REVOCAP_Coupler
# -DWITH_PARACON     : for parallel contact (unimplemented)
# -DWITH_METIS       : compile with METIS graph partitioning package
#   -DMETIS_VER_4    : uses METIS-4
# -DWITH_PARMETIS    : compile with ParMETIS (unimplemented)
#   -DPARMETIS_VER_3 : uses PARMETIS-3
# -DWITH_MKL         : compile with MKL PARDISO (unimplemented)
# -DWITH_MUMPS       : compile with MUMPS direct solver
# -DWITH_LAPACK      : compile with LAPACK
# -DWITH_ML          : compile with Trilinos multigrid preconditioning package
# -DBUILD_DOC        : generate API documents

###############################################################################
# Archtecture and compiler dependent options
###############################################################################
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(BLA_VENDOR "Intel")
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Windows")
  #
  # MinGW specific options
  #
  if(CMAKE_GENERATOR MATCHES "(MinGW|MSYS)")
    list(APPEND FrontISTR_DEFINITIONS "_WINDOWS")
    set(CMAKE_Fortran_FLAGS "-fno-range-check")
    set(WINDOWS ON)
  endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(CMAKE_C_COMPILER "gcc")
  set(CMAKE_CXX_COMPILER "g++")
  set(CMAKE_Fortran_COMPILER "gfortran")
else()
endif()

###############################################################################
# Build type (RELEASE or DEBUG)
###############################################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Build type (RELEASE or DEBUG)")
endif()
if(CMAKE_BUILD_TYPE MATCHES "DEBUG")
  list(APPEND FrontISTR_DEFINITIONS "DEBUG")
elseif (CMAKE_BUILD_TYPE MATCHES "RELEASE")
  list(APPEND FrontISTR_DEFINITIONS "NDEBUG")
endif()

###############################################################################
# Find external libraries and set WITH_XX=ON if found libraries
###############################################################################
# Find external packages
find_package(LAPACK)
find_package(Metis)
find_package(MPI)
find_package(Mumps)
find_package(OpenMP)
find_package(Parmetis)
find_package(Refiner)
find_package(Revocap)
find_package(Scalapack)
find_package(Trilinos)
find_package(Doxygen)
feature_summary(WHAT ALL)

option(WINDOWS "build on windows" ${WINDOWS})
option(WITH_TOOLS "compile tools" ON)
option(WITH_MPI "for parallel environment with MPI" ${MPI_C_FOUND})
option(WITH_OPENMP "for multi-(core|processor) environment" ${OPENMP_FOUND})
option(WITH_REFINER "compile with REVOCAP_Refiner" ${REFINER_FOUND})
option(WITH_REVOCAP "compile with REVOCAP_Coupler" ${REVOCAP_FOUND})
option(WITH_PARACON "for parallel contact" OFF)
option(WITH_METIS "compile with METIS" ${METIS_FOUND})
option(METIS_VER_4 "uses METIS-4" OFF)
option(WITH_PARMETIS "compile with ParMETIS" ${PARMETIS_FOUND})
option(PARMETIS_VER_3 "uses METIS-3" OFF)
option(WITH_MUMPS "compile with MUMPS direct solver" ${MUMPS_FOUND})
option(WITH_ML "compile with Trilinos ML" ${Trilinos_FOUND})
option(BUILD_DOC "Generate API documents." ${DOXYGEN_FOUND})

##############
# -DWINDOWS
##############
if(WINDOWS)
  list(APPEND FrontISTR_DEFINITIONS "_WINDOWS")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()
# -DWITH_MPI
if(WITH_MPI)
  list(APPEND FrontISTR_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
  list(APPEND FrontISTR_LIBRARIES ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
else()
  set(WITH_MPI OFF)
  list(APPEND FrontISTR_DEFINITIONS "HECMW_SERIAL")
endif()


# -DWITH_OPENMP
if(WITH_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
else()
  set(WITH_OPENMP OFF)
endif()

# -DWITH_REFINER
if(WITH_REFINER)
  list(APPEND FrontISTR_INCLUDE_DIRS ${REFINER_INCLUDE_PATH})
  list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_REFINER")
  list(APPEND FrontISTR_LIBRARIES ${REFINER_LIBRARIES})
else()
  set(WITH_REFINER OFF)
  set(HECMW_WITH_REFINER OFF)
  list(REMOVE_ITEM FrontISTR_DEFINITIONS "HECMW_WITH_REFINER")
endif()

# -DWITH_REVOCAP
if(WITH_REVOCAP)
  list(APPEND FrontISTR_INCLUDE_DIRS ${REVOCAP_INCLUDE_PATH})
else()
  set(WITH_REVOCAP OFF)
endif()

# -DWITH_PARACON
if(WITH_PARACON)
  list(APPEND FrontISTR_DEFINITIONS "PARA_CONTACT")
else()
  set(WITH_PARACON OFF)
  list(REMOVE_ITEM FrontISTR_DEFINITIONS "PARA_CONTACT")
endif()

# -DWITH_METIS, -DMETIS_VER_4
if(WITH_METIS)
  if(METIS_FOUND)
    list(APPEND FrontISTR_INCLUDE_DIRS ${METIS_INCLUDE_PATH})
    list(APPEND FrontISTR_DEFINITIONS "HECMW_WITH_METIS")
    list(APPEND FrontISTR_DEFINITIONS "HECMW_PART_WITH_METIS")
    list(APPEND FrontISTR_LIBRARIES ${METIS_LIBRARIES})
    if(NOT METIS_VER_4)
      list(APPEND FrontISTR_DEFINITIONS "HECMW_METIS_VER=5")
    else()
      list(APPEND FrontISTR_DEFINITIONS "HECMW_METIS_VER=4")
    endif()
  else()
    list(REMOVE_ITEM FrontISTR_DEFINITIONS "HECMW_WITH_METIS")
    list(REMOVE_ITEM FrontISTR_DEFINITIONS "HECMW_PART_WITH_METIS")
  endif()
endif()

# -DWITH_PARMETIS
if(WITH_PARMETIS)
  if(PARMETIS_FOUND)
    list(APPEND FrontISTR_INCLUDE_DIRS ${PARMETIS_INCLUDE_PATH})
    list(APPEND FrontISTR_LIBRARIES ${PARMETIS_LIBRARIES})
  endif()
endif()

# -DWITH_MKL (unimplemented)
# W.I.P.
option(WITH_MKL "compile with MKL PARDISO" OFF)
if(WITH_MKL)
  get_filename_component(ICC icc PROGRAM)
  set(CMAKE_C_COMPILER ${ICC} CACHE FILEPATH "Intel C Compiler" FORCE)
  get_filename_component(ICPC icpc PROGRAM)
  set(CMAKE_CXX_COMPILER ${ICPC} CACHE FILEPATH "Intel C++ Compiler" FORCE)
  get_filename_component(IFORT ifort PROGRAM)
  set(CMAKE_Fortran_COMPILER ${IFORT} CACHE FILEPATH "Intel Fortran Compiler" FORCE)
  set(BLA_VENDOR "Intel")
  find_package(LAPACK)
  find_package(Scalapack)
  add_definitions(-DWITH_MKL)
endif()
#set(WITH_LAPACK ON)

# -DWITH_MUMPS
if(WITH_MUMPS)
  if(MUMPS_FOUND)
    list(APPEND FrontISTR_DEFINITIONS "WITH_MUMPS")
    list(APPEND FrontISTR_INCLUDE_DIRS ${MUMPS_INCLUDE_PATH})
    list(APPEND FrontISTR_LIBRARIES ${MUMPS_LIBRARIES})
    mark_as_advanced(scalapack_DIR)
  else()
    message(FATAL_ERROR "MUMPS could not found")
  endif()
else()
  set(WITH_MUMPS OFF)
endif()

# -DWITH_LAPACK
if(WITH_LAPACK)
  if(LAPACK_FOUND AND BLAS_FOUND)
    add_definitions(-DWITH_LAPACK)
  else()
    set(BLAS_LIBRARIES "BLAS_LIBRARIES-NOTFOUND" CACHE FILEPATH "BLAS Library")
    set(LAPACK_LIBRARIES "LAPACK_LIBRARIES-NOTFOUND" CACHE FILEPATH "LAPACK Library")
    message(WARNING "Please specify -DBLAS_LIBRARIES and -DLAPACK_LIBRARIES manually.")
    set(WITH_LAPACK OFF)
    remove_definitions(-DWITH_LAPACK)
  endif()
endif()

# -DWITH_ML
if(WITH_ML)
  if(Trilinos_FOUND)
    list(APPEND FrontISTR_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
    list(APPEND FrontISTR_DEFINITIONS HECMW_WITH_ML)
    list(APPEND FrontISTR_LIBRARIES ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  else()
    #message(WARNING "Please turn on WITH_MPI")
    message(WARNING "Trilinos not found")
    set(WITH_ML OFF)
    list(REMOVE_ITEM FrontISTR_DEFINITIONS HECMW_WITH_ML)
  endif()
endif()

# -DBUILD_DOC
if(BUILD_DOC)
  if(DOXYGEN_FOUND)
    message(STATUS "Generate API documents. Please type \"make doc\"")
    file(COPY ${CMAKE_SOURCE_DIR}/doc/fistr.css DESTINATION ${CMAKE_BINARY_DIR}/doc/)
    file(COPY ${CMAKE_SOURCE_DIR}/doc/fistr_logo.png DESTINATION ${CMAKE_BINARY_DIR}/doc/)

    set(doxy_main_page "${CMAKE_SOURCE_DIR}/INSTALL.md")
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
message(STATUS "C compiler is ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER}")
message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER}")
message(STATUS "Installation path is ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_PREFIX_PATH})
  message(STATUS "Library search path is ${CMAKE_PREFIX_PATH}")
else()
  message(STATUS "Library search path is $ENV{HOME}/local /usr/local /usr etc.")
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
add_subdirectory(examples)

###############################################################################
# Package generator (deb, rpm, exe etc.)
###############################################################################
include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Large-scale structural analysis program with FEM ")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VENDOR "FrontISTR Forum, The University of Tokyo")
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
