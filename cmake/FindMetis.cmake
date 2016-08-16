###############################################################################
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# METIS_FOUND
# METIS_INCLUDE_PATH
# METIS_LIBRARIES
# METIS_VER_4
#
if(NOT METIS_VER_4)
  find_path(METIS_INCLUDE_PATH
    NAMES metis.h
    PATHS ${CMAKE_SOURCE_DIR}/../metis-5.1.0/build/include
    $ENV{HOME}/metis-5.0.1/build/include
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    /usr/local/metis/include
    /usr/local/include
    /usr/metis/include
    /usr/include
  )
  find_library(METIS_LIBRARIES
    NAMES metis
    PATHS ${CMAKE_SOURCE_DIR}/../metis-5.1.0/build/libmetis
    $ENV{HOME}/metis-5.0.1/build/libmetis
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    /usr/local/metis/lib
    /usr/local/lib
    /usr/metis/lib
    /usr/lib
  )
else()
  find_path(METIS_INCLUDE_PATH
    NAMES metis.h
    PATH ${CMAKE_SOURCE_DIR}/../metis-4.0.3/Lib
    $ENV{HOME}/metis-4.0.3/
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    /usr/local/metis/include
    /usr/local/include
    /usr/metis/include
    /usr/include
  )
  find_library(METIS_LIBRARIES
    NAMES metis
    PATH ${CMAKE_SOURCE_DIR}/../metis-4.0.3
    $ENV{HOME}/metis-4.0.3
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    /usr/local/metis/lib
    /usr/local/lib
    /usr/metis/lib
    /usr/lib
  )
endif()

if(METIS_INCLUDE_PATH AND METIS_LIBRARIES)
  set(METIS_FOUND TRUE)
endif()
