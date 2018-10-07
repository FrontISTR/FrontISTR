###############################################################################
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

#
# Variables:
#
# SCALAPACK_FOUND      TRUE if FindScalapack found scalapack
# SCALAPACK_MKL        TRUE if FindScalapack found Intel MKL scalapack
# SCALAPACK_LIBRARIES  scalapack libraries
#
# env SCALAPACK_ROOT   Set SCALAPACK_ROOT environment variable,
#                      where scalapack are.
#    ex. export SCALAPACK_ROOT=/home/someone/somewhere/scalapack-2.0.2
#
if(SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND TRUE)
  RETURN()
endif()

include(FindPackageHandleStandardArgs)

if(EXISTS $ENV{MKLROOT})
  find_file(_MKL_SCALAPACK_INCLUDE
    NAMES "mkl_scalapack.h"
    HINTS "$ENV{MKLROOT}/include"
    NO_SYSTEM_ENVIRONMENT_PATH)

  if(_MKL_SCALAPACK_INCLUDE)
    set(SCALAPACK_INCLUDE_PATH ${_MKL_SCALAPACK_INCLUDE})
    find_library(_MKL_SCALAPACK_LP64
      NAMES mkl_scalapack_lp64
      HINTS $ENV{MKLROOT}/lib/intel64)

    find_library(_MKL_INTEL_LP64
      NAMES mkl_intel_lp64
      HINTS $ENV{MKLROOT}/lib/intel64)

    find_library(_MKL_INTEL_THREAD
      NAMES mkl_intel_thread
    HINTS $ENV{MKLROOT}/lib/intel64)

    find_library(_MKL_GNU_THREAD
      NAMES mkl_gnu_thread
      HINTS $ENV{MKLROOT}/lib/intel64)

    find_library(_MKL_CORE
      NAMES mkl_core
      HINTS $ENV{MKLROOT}/lib/intel64)

    find_library(_MKL_BLACS_INTELMPI_LP64
      NAMES mkl_blacs_intelmpi_lp64
      HINTS $ENV{MKLROOT}/lib/intel64)

    find_library(_MKL_BLACS_OPENMPI_LP64
      NAMES mkl_blacs_openmpi_lp64
      HINTS $ENV{MKLROOT}/lib/intel64)

    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
      set(SCALAPACK_LIBRARIES
        ${_MKL_SCALAPACK_LP64}
        ${_MKL_INTEL_LP64}
        ${_MKL_INTEL_THREAD}
        ${_MKL_CORE}
        ${_MKL_BLACS_INTELMPI_LP64}
        iomp5
        pthread
        m
        dl
        CACHE STRING "MKL ScaLAPACK for Intel")

    elseif(CMAKE_C_COMPILER_ID MATCHES "GNU")
      set(SCALAPACK_LIBRARIES
        ${_MKL_SCALAPACK_LP64}
        ${_MKL_INTEL_LP64}
        ${_MKL_GNU_THREAD}
        ${_MKL_CORE}
        ${_MKL_BLACS_INTELMPI_LP64}
        gomp
        pthread
        m
        dl
        CACHE STRING "MKL ScaLAPACK for GCC")
    endif()
    set(SCALAPACK_MKL ON)
    find_package_handle_standard_args(Scalapack
      DEFAULT_MSG SCALAPACK_LIBRARIES SCALAPACK_INCLUDE_PATH)
  endif()
else()
  find_library(SCALAPACK_LIBRARIES
    NAMES scalapack scalapack-openmpi
    HINTS $ENV{SCALAPACK_ROOT}/build/lib
    ${CMAKE_SOURCE_DIR}/../scalapack-2.0.2/build/lib
    $ENV{HOME}/scalapack-2.0.2/build/lib
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/scalapack/lib
    /usr/local/lib
    /usr/lib
    /usr/lib/x86_64-linux-gnu)
  unset(WITH_MKL)
  find_package_handle_standard_args(Scalapack
    DEFAULT_MSG SCALAPACK_LIBRARIES)
endif()

if(SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND TRUE)
endif()

mark_as_advanced(_MKL_SCALAPACK_LP64 _MKL_INTEL_LP64
  _MKL_GNU_THREAD _MKL_INTEL_THREAD _MKL_CORE
  _MKL_BLACS_INTELMPI_LP64 _MKL_SCALAPACK_INCLUDE _MKL_BLACS_OPENMPI_LP64)

