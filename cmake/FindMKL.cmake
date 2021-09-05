###############################################################################
# Copyright (c) 2020 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# MKL_FOUND         TRUE if FindMKL found MKL
# MKL_INCLUDE_PATH  Inclue path of MKL
# MKL_LIBRARIES     MKL libraries
#
# env MKL_ROOT      Set MKL_ROOT environment variable,
#
if(MKL_LIBRARIES)
  set(MKL_FOUND TRUE)
  RETURN()
endif()

find_path(MKL_INCLUDE_PATH
  NAMES mkl.h
  HINTS $ENV{MKLROOT}/include
  $ENV{HOME}/local/include
  $ENV{HOME}/.local/include
  ${CMAKE_INCLUDE_PATH}
  /opt/intel/mkl/include
  /usr/local/include/mkl
  /usr/include/mkl
  /usr/local/include
  /usr/include
)
find_library(_MKL_INTEL_LP64          NAMES mkl_intel_lp64           HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
find_library(_MKL_INTEL_THREAD        NAMES mkl_intel_thread         HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
find_library(_MKL_GNU_THREAD          NAMES mkl_gnu_thread           HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
find_library(_MKL_CORE                NAMES mkl_core                 HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(MKL_LIBRARIES
    ${_MKL_INTEL_LP64}
    ${_MKL_INTEL_THREAD}
    ${_MKL_CORE}
    iomp5
    pthread
    m
    dl
    CACHE STRING "MKL for Intel")

elseif(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(MKL_LIBRARIES
    ${_MKL_INTEL_LP64}
    ${_MKL_GNU_THREAD}
    ${_MKL_CORE}
    gomp
    pthread
    m
    dl
    CACHE STRING "MKL for GCC")
endif()
if(MKL_INCLUDE_PATH AND MKL_LIBRARIES)
  set(MKL_FOUND TRUE)
endif()

mark_as_advanced(MKL_INCLUDE_PATH MKL_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_PATH)

