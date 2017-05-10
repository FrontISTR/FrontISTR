###############################################################################
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

# please use scalapack-config.cmake
# Exists in ${SCALAPACK_INSTALL_DIR}/lib/cmake/scalapack-2.0.2/
#
# Variables:
#
# SCALAPACK_FOUND
# SCALAPACK_LIBRARIES
#
if(SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND TRUE)
  RETURN()
endif()

if(WITH_MKL)
  message(STATUS "LAPACK Vendor is ${BLA_VENDOR}")
  find_library(_MKL_SCALAPACK
    NAMES mkl_scalapack_lp64
    HINTS $ENV{MKLROOT}/lib/inel64
  )
  find_library(_MKL_BLACS_INTELMPI
    NAMES mkl_blacs_intelmpi_lp64
    HINTS $ENV{MKLROOT}/lib/intel64
  )
  set(SCALAPACK_LIBRARIES ${_MKL_SCALAPACK} ${_MKL_BLACS_INTELMPI} FORCE)
else()
  find_library(SCALAPACK_LIBRARIES
    NAMES scalapack
    HINTS ${CMAKE_SOURCE_DIR}/../scalapack-2.0.2/build/lib
    $ENV{HOME}/metis-5.0.1/build/lib
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/metis/lib
    /usr/local/lib
    /usr/lib
  )
endif()

if(SCALAPACK_LIBRARIES)
  message(STATUS "Found SCALAPACK")
  set(SCALAPACK_FOUND TRUE)
endif()
