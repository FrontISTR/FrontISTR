###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# PARMETIS_FOUND         TRUE if FindParmetis found parmetis
# PARMETIS_INCLUDE_PATH  Include path of parmetis
# PARMETIS_LIBRARIES     parmetis libraries
#
# PARMETIS_VER_3         Set TRUE, if use parmetis-3
# env PARMETIS_ROOT      Set PARMETIS_ROOT envioenment variable,
#                        where parmetis are.
#    ex. export PARMETIS_ROOT=/home/someone/somewhere/parmetis-4.0.3
#
if(PARMETIS_LIBRARIES)
  set(PARMETIS_FOUND TRUE)
  RETURN()
endif()

if(NOT PARMETIS_VER_3)
  find_path(PARMETIS_INCLUDE_PATH
    NAMES parmetis.h
    HINTS ${PARMETIS_ROOT}/include
    ${CMAKE_SOURCE_DIR}/../parmetis-4.0.3/include
    $ENV{HOME}/parmetis-4.0.3/include
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/parmetis/include
    /usr/local/include
    /usr/parmetis/include
    /usr/include
  )
  find_library(PARMETIS_LIBRARIES
    NAMES parmetis
    HINTS $ENV{PARMETIS_ROOT}/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libparmetis
    ${CMAKE_SOURCE_DIR}/../parmetis-4.0.3/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libparmetis
    $ENV{HOME}/parmetis-4.0.3/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libparmetis
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/parmetis/lib
    /usr/local/lib
    /usr/parmetis/lib
    /usr/lib
  )
else()
  find_path(PARMETIS_INCLUDE_PATH
    NAMES parmetis.h
    HINTS $ENV{PARMETIS_ROOT}/Lib
    ${CMAKE_SOURCE_DIR}/../ParMetis-3.2.0/Lib
    $ENV{HOME}/ParMetis-3.2.0/
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/parmetis/include
    /usr/local/include
    /usr/parmetis/include
    /usr/include
  )
  find_library(PARMETIS_LIBRARIES
    NAMES parmetis
    HINTS $ENV{PARMETIS_ROOT}
    ${CMAKE_SOURCE_DIR}/../ParMetis-3.2.0
    $ENV{HOME}/ParMetis-3.2.0
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/parmetis/lib
    /usr/local/lib
    /usr/parmetis/lib
    /usr/lib
  )
endif()

if(PARMETIS_INCLUDE_PATH AND PARMETIS_LIBRARIES)
  set(PARMETIS_FOUND TRUE)
endif()

mark_as_advanced(PARMETIS_INCLUDE_PATH PARMETIS_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Parmetis
  DEFAULT_MSG PARMETIS_LIBRARIES PARMETIS_INCLUDE_PATH)
