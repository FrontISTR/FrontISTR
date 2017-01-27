###############################################################################
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# PARMETIS_FOUND
# PARMETIS_INCLUDE_PATH
# PARMETIS_LIBRARIES
# PARMETIS_VER_4
#
if(PARMETIS_LIBRARIES)
  set(PARMETIS_FOUND TRUE)
  RETURN()
endif()

if(NOT PARMETIS_VER_3)
  find_path(PARMETIS_INCLUDE_PATH
    NAMES parmetis.h
    HINTS ${CMAKE_SOURCE_DIR}/../parmetis-4.0.3/include
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
    HINTS ${CMAKE_SOURCE_DIR}/../parmetis-4.0.3/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libparmetis
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
    PATH ${CMAKE_SOURCE_DIR}/../ParMetis-3.2.0/Lib
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
    PATH ${CMAKE_SOURCE_DIR}/../ParMetis-3.2.0
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
  message(STATUS "Found parmetis")
  set(PARMETIS_FOUND TRUE)
endif()

mark_as_advanced(PARMETIS_INCLUDE_PATH PARMETIS_LIBRARIES)
