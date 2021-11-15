###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# REFINER_FOUND        TRUE if FindRefiner found REVOCAP_Refiner
# REFINER_INCLUDE_PATH Include path of REVOCAP_Refiner
# REFINER_LIBRARIES    REVOCAP_Refiner libraries
#
# env REFINER_ROOT     Set REFINER_ROOT environment variable,
#                      where REVOCAP_Refiner are.
#    ex. export REFINER_ROOT=/home/someone/somewhere/REVOCAP_Refiner-1.1.04
#
if(REFINER_LIBRARIES)
  set(REFINER_FOUND TRUE)
  RETURN()
endif()

string(TOLOWER ${CMAKE_SYSTEM_NAME} REFINER_SYSTEM_NAME)
find_path(REFINER_INCLUDE_PATH
  NAMES rcapRefiner.h
  HINTS $ENV{REFINER_ROOT}/Refiner
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Refiner-1.1.04/Refiner
    $ENV{HOME}/REVOCAP_Refiner-1.1.04/Refiner
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    NO_DEFAULT_PATH
)
find_library(REFINER_LIBRARIES
  NAMES RcapRefiner
  HINTS $ENV{REFINER_ROOT}/lib/${CMAKE_SYSTEM_PROCESSOR}-${REFINER_SYSTEM_NAME}
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Refiner-1.1.04/lib/${CMAKE_SYSTEM_PROCESSOR}-${REFINER_SYSTEM_NAME}
    $ENV{HOME}/REVOCAP_Refiner-1.1.04/lib/${CMAKE_SYSTEM_PROCESSOR}-${REFINER_SYSTEM_NAME}
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/lib
    /usr/lib
    NO_DEFAULT_PATH
)

if(REFINER_INCLUDE_PATH AND REFINER_LIBRARIES)
  set(REFINER_FOUND TRUE)
  set(HECMW_WITH_REFINER)
endif()

mark_as_advanced(REFINER_INCLUDE_PATH REFINER_LIBRARIES REFINER_SYSTEM_NAME)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(REVOCAP_Refiner
  DEFAULT_MSG REFINER_LIBRARIES REFINER_INCLUDE_PATH)

