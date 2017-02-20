###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
message(STATUS ${CMAKE_SYSTEM_NAME}${CMAKE_SYSTEM_PROCESSOR}${CMAKE_SYSTEM_VERSION})
###############################################################################
# CMakeLists.txt for FrontISTR
#
# FrontISTR
# ├── CMakeLists.txt
# ├── cmake
# │   ├── FindMetis.cmake
# │   ├── FindMumps.cmake
# │   ├── FindParMetis.cmake
# │   ├── FindRefiner.cmake
# │   ├── FindRevocap.cmake
# │   └── FindScalapack.cmake
# ├── examples
# │   ├── CMakeLists.txt
# │   └── test_FrontISTR.rb
# ├── fistr
# │   └── CMakeLists.txt
# ├── fistr1
# │   ├── CMakeLists.txt
# │   ├── src
# │   │   └── CMakeLists.txt
# │   └── tools
# │       └── CMakeLists.txt
# └── hecmw1
#      ├── CMakeLists.txt
#      ├── src
#      │   └── CMakeLists.txt
#      └── tools
#         └── CMakeLists.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)
include(CMakeDependentOption)

###############################################################################
# Build type (Release or Debug)
###############################################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Build type (RELEASE or DEBUG)")
endif()

if(CMAKE_BUILD_TYPE MATCHES "RELEASE")
  add_definitions(-DNDEBUG)
  set(NDEBUG ON)
elseif(CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions(-DDEBUG)
  set(DEBUG ON)
endif()

###############################################################################
# Project name
###############################################################################
project(FrontISTR Fortran C CXX)

###############################################################################
# Version number
###############################################################################
set(VERSION_MAJOR 5 CACHE STRING "Major version number of FrontISTR")
set(VERSION_MINOR 0 CACHE STRING "Minor version number of FrontISTR")
set(VERSION_PATCH 0 CACHE STRING "Patch version number of FrontISTR")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)


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

# includes standard macros
include(FeatureSummary)
include(CMakeDependentOption)

# add module path for FindPackage.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

# Find external packages
find_package(MPI)
find_package(OpenMP)
find_package(Refiner)
find_package(Revocap)
find_package(Metis)
find_package(ParMetis)
find_package(Mumps)
find_package(scalapack)
find_package(LAPACK)
find_package(Trilinos)
find_package(Doxygen)
feature_summary(WHAT ALL)
message(STATUS "Installation path is ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_PREFIX_PATH})
  message(STATUS "Library search path is ${CMAKE_PREFIX_PATH}")
else()
  message(STATUS "Library search path is $ENV{HOME}/local /usr/local /usr etc.")
endif()

# Options
###############################################################################
# Windows (MinGW) specific setting
###############################################################################
if(CMAKE_GENERATOR MATCHES "(MinGW|MSYS)")
  set(WINDOWS ON)
endif()
option(WINDOWS "build on windows" ${WINDOWS})
option(WITH_TOOLS "compile tools" ON)
option(WITH_MPI "for parallel environment with MPI" ${MPI_Fortran_FOUND})
option(WITH_OPENMP "for multi-(core|processor) environment" ${OPENMP_FOUND})
option(WITH_REFINER "compile with REVOCAP_Refiner" ${REFINER_FOUND})
option(WITH_REVOCAP "compile with REVOCAP_Coupler" ${REVOCAP_FOUND})
option(WITH_PARACON "for parallel contact" OFF)
option(WITH_METIS "compile with METIS" ${METIS_FOUND})
option(METIS_VER_4 "uses METIS-4" OFF)
option(WITH_PARMETIS "compile with ParMETIS" ${PARMETIS_FOUND})
option(PARMETIS_VER_3 "uses METIS-3" OFF)
option(WITH_MUMPS "compile with MUMPS direct solver" ${MUMPS_FOUND})
option(WITH_ML "compile with Trilinos multigrid preconditioning package" ${Trilinos_FOUND})
option(BUILD_DOC "Generate API documents." ${DOXYGEN_FOUND})

# -DWINDOWS
if(WINDOWS)
  add_definitions(-D_WINDOWS)
  set(CMAKE_Fortran_FLAGS "-fno-range-check")
  set(_WINDOWS ON)
endif()

# -DWITH_MPI
if(WITH_MPI)
  if(MPI_Fortran_FOUND AND MPI_CXX_FOUND AND MPI_C_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH})
    set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortra_LINK_FLAGS})
    mark_as_advanced(MPI_EXTRA_LIBRARY MPI_LIBRARY)
  else()
    set(WITH_MPI OFF)
  endif()
else()
  add_definitions(-DHECMW_SERIAL)
  set(HECMW_SERIAL ON)
endif()

# -DWITH_OPENMP
if(WITH_OPENMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
  else()
    set(WITH_OPENMP OFF)
  endif()
endif()

# -DWITH_REFINER
if(WITH_REFINER)
  if(REFINER_FOUND)
    include_directories(${REFINER_INCLUDE_PATH})
    add_definitions(-DHECMW_WITH_REFINER)
    set(HECMW_WITH_REFINER ON)
  else()
    set(WITH_REFINER OFF)
    set(HECMW_WITH_REFINER OFF)
    remove_definitions(-DHECMW_WITH_REFINER)
  endif()
endif()

# -DWITH_REVOCAP
if(WITH_REVOCAP)
  if(REVOCAP_FOUND)
    include_directories(${REVOCAP_INCLUDE_PATH})
  else()
    set(WITH_REVOCAP OFF)
  endif()
endif()

# -DWITH_PARACON (unimplemented)
if(WITH_PARACON)
  add_definitions(-DPARA_CONTACT)
  set(PARA_CONTACT ON)
else()
  set(WITH_PARACON OFF)
  set(PARA_CONTACT OFF)
  remove_definitions(-DPARA_CONTACT)
endif()

# -DWITH_METIS, -DMETIS_VER_4
if(WITH_METIS)
  if(METIS_FOUND)
    include_directories(${METIS_INCLUDE_PATH})
    add_definitions(-DHECMW_WITH_METIS)
    add_definitions(-DHECMW_PART_WITH_METIS)
    set(HECMW_WITH_METIS ON)
    set(HECMW_PART_WITH_METIS ON)
    if(NOT METIS_VER_4)
      add_definitions(-DHECMW_METIS_VER=5)
      set(HECMW_METIS_VER 5)
    else()
      add_definitions(-DHECMW_METIS_VER=4)
      set(HECMW_METIS_VER 4)
    endif()
  else()
    remove_definitions(-DHECMW_WITH_METIS)
    remove_definitions(-DHECMW_PART_WITH_METIS)
    set(HECMW_WITH_METIS OFF)
    set(HECMW_PART_WITH_METIS OFF)
  endif()
endif()

# -DWITH_PARMETIS (unimplemented)
if(WITH_PARMETIS)
  if(PARMETIS_FOUND)
    include_directories(${PARMETIS_INCLUDE_PATH})
  endif()
endif()

# -DWITH_MKL (unimplemented)
# W.I.P.
#option(WITH_MKL "compile with MKL PARDISO" OFF)
#set(BLA_VENDOR "Intel")
#set(WITH_LAPACK ON)

# -DWITH_MUMPS
if(WITH_MUMPS)
  if(MUMPS_FOUND)
    add_definitions(-DWITH_MUMPS)
    include_directories(${MUMPS_INCLUDE_PATH})
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
    include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
    set(HECMW_WITH_ML ON)
    add_definitions(-DHECMW_WITH_ML)
    mark_as_advanced(Trilinos_DIR)
  else()
    #message(WARNING "Please turn on WITH_MPI")
    message(WARNING "Trilinos not found")
    set(WITH_ML OFF)
    remove_definitions(-DHECMW_WITH_ML)
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
add_subdirectory(fistr)

# libfistr.a depends on libhecmw.a
add_dependencies(fistr hecmw)
# fistr1 depends on libfistr.a and libhecmw.a
add_dependencies(fistr1 fistr hecmw)

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
