###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# METIS_FOUND         TRUE if FindMetis found metis
# METIS_INCLUDE_PATH  Include path of metis
# METIS_LIBRARIES     metis libraries
#
# METIS_VER_4         Set TRUE, if use metis-4
# env METIS_ROOT      Set METIS_ROOT environment variable,
#                     where metis are.
#    ex. export METIS_ROOT=/home/someone/somewhere/metis-5.1.0
#
if(METIS_LIBRARIES)
  set(METIS_FOUND TRUE)
  RETURN()
endif()

if(NOT METIS_VER_4)
  find_path(METIS_INCLUDE_PATH
    NAMES metis.h
    HINTS $ENV{METIS_ROOT}/include
    ${CMAKE_SOURCE_DIR}/../metis-5.1.0/include
    $ENV{HOME}/metis-5.1.0/include
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/metis/include
    /usr/local/include
    /usr/metis/include
    /usr/include
  )
  find_library(METIS_LIBRARIES
    NAMES metis
    HINTS $ENV{METIS_ROOT}/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
    ${CMAKE_SOURCE_DIR}/../metis-5.1.0/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
    $ENV{HOME}/metis-5.1.0/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/metis/lib
    /usr/local/lib
    /usr/metis/lib
    /usr/lib
  )
else()
  find_path(METIS_INCLUDE_PATH
    NAMES metis.h
    HINTS $ENV{METIS_ROOT}/Lib
    ${CMAKE_SOURCE_DIR}/../metis-4.0.3/Lib
    $ENV{HOME}/metis-4.0.3/
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/metis/include
    /usr/local/include
    /usr/metis/include
    /usr/include
  )
  find_library(METIS_LIBRARIES
    NAMES metis
    HINTS $ENV{METIS_ROOT}
    ${CMAKE_SOURCE_DIR}/../metis-4.0.3
    $ENV{HOME}/metis-4.0.3
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/metis/lib
    /usr/local/lib
    /usr/metis/lib
    /usr/lib
  )
endif()

if(METIS_INCLUDE_PATH AND METIS_LIBRARIES)
  set(METIS_FOUND TRUE)
endif()

mark_as_advanced(METIS_INCLUDE_PATH METIS_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Metis
  DEFAULT_MSG METIS_LIBRARIES METIS_INCLUDE_PATH)

