###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# env BLOPEX_ROOT    Set BLOPEX_ROOT environment variable, where metis are.
# ex. export BLOPEX_ROOT=/home/someone/somewhere/BLOPEX

if(BLOPEX_LIBRARIES)
  set(BLOPEX_FOUND TRUE)
  RETURN()
endif()

set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} $ENV{BLOPEX_ROOT}/blopex_abstract/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} $ENV{BLOPEX_ROOT}/blopex_fortran)

find_library(BLOPEX_LIB
  NAMES BLOPEX
  HINTS ${LIB_SEARCH_PATH}
)

find_library(BLOPEX_DRIVER_LIB
  NAMES BLOPEX_driver
  HINTS ${LIB_SEARCH_PATH}
)

find_path(BLOPEX_INCLUDE_PATH
  NAMES lobpcg.h
  HINTS $ENV{BLOPEX_ROOT}/blopex_abstract/include
)

if(BLOPEX_DRIVER_LIB AND BLOPEX_LIB AND BLOPEX_INCLUDE_PATH)
  set(BLOPEX_LIBRARIES ${BLOPEX_DRIVER_LIB} ${BLOPEX_LIB})
  set(BLOPEX_FOUND ON)
endif()
mark_as_advanced(BLOPEX_INCLUDE_PATH BLOPEX_LIBRARIES BLOPEX_DRIVER_LIB BLOPEX_LIB)
