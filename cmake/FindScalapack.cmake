###############################################################################
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# SCALAPACK_FOUND
# SCALAPACK_LIBRARIES
#
if(SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND TRUE)
  RETURN()
endif()

set(SCALAPACK_LIB
  "${SCALAPACK_ROOT}"
  "${SCALAPACK_ROOT}/lib"
  "$ENV{SCALAPACK_ROOT}"
  "$ENV{SCALAPACK_ROOT}/lib"
  "${CMAKE_LIBRARY_PREFIX}"
  "${CMAKE_SOURCE_DIR}/scalapack"
  "${CMAKE_SOURCE_DIR}/scalapack/lib")

find_library(SCALAPACK_LIBRARIES
  NAMES
  "scalapack" "scalapack-mpi"
  HINTS
  ${SCALAPACKLIB})

if(SCALAPACK_LIBRARIES)
  message(STATUS "Found SCALAPACK")
  set(SCALAPACK_FOUND TRUE)
endif()
