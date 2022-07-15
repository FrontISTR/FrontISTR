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

find_library(BLOPEX_LIB
  NAMES BLOPEX
  HINTS ${LIB_SEARCH_PATH}
)

find_library(BLOPEX_DRIVER_LIB
  NAMES BLOPEX_driver
  HINTS ${LIB_SEARCH_PATH}
)

find_path(BLOPEX_INCLUDE_1_PATH
  NAMES lobpcg.h
  HINTS $ENV{BLOPEX_ROOT}/blopex_abstract/include
)

find_path(BLOPEX_INCLUDE_2_PATH
  NAMES multi_vector.h
  HINTS $ENV{BLOPEX_ROOT}/blopex_serial_double/multivector
)

find_path(BLOPEX_INCLUDE_3_PATH
  NAMES pcg_multi.h
  HINTS $ENV{BLOPEX_ROOT}/blopex_serial_double/pcg_multi
)

if(BLOPEX_LIB AND BLOPEX_DRIVER_LIB AND BLOPEX_INCLUDE_1_PATH AND BLOPEX_INCLUDE_2_PATH AND BLOPEX_INCLUDE_3_PATH)
  list(APPEND BLOPEX_INCLUDE_PATH ${BLOPEX_INCLUDE_1_PATH} ${BLOPEX_INCLUDE_2_PATH} ${BLOPEX_INCLUDE_3_PATH})
  set(BLOPEX_LIBRARIES ${BLOPEX_DRIVER_LIB} ${BLOPEX_LIB})
  set(BLOPEX_FOUND ON)
endif()

mark_as_advanced(BLOPEX_LIBRARIES BLOPEX_INCLUDE_PATH)
