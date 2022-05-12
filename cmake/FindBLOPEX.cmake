###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

if(BLOPEX_LIBRARIES)
  set(BLOPEX_FOUND TRUE)
  RETURN()
endif()

set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} /Users/morita/git/blopex/blopex_abstract/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} /Users/morita/git/blopex/blopex_fortran)

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
  HINTS /Users/morita/git/blopex/blopex_abstract/include
)

if(BLOPEX_DRIVER_LIB AND BLOPEX_LIB AND BLOPEX_INCLUDE_PATH)
  set(BLOPEX_LIBRARIES ${BLOPEX_DRIVER_LIB} ${BLOPEX_LIB})
  set(BLOPEX_FOUND ON)
endif()
mark_as_advanced(BLOPEX_INCLUDE_PATH BLOPEX_LIBRARIES BLOPEX_DRIVER_LIB BLOPEX_LIB)
