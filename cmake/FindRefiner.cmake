###############################################################################
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# env REVINER_ROOT
# REFINER_FOUND
# REFINER_INCLUDE_PATH
# REFINER_LIBRARIES
#
find_path(REFINER_INCLUDE_PATH
  NAMES rcapRefiner.h
  PATHS $ENV{REFINER_ROOT}/Refiner
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Refiner-1.1.04/Refiner
    $ENV{HOME}/REVOCAP_Refiner-1.1.04/Refiner
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    /usr/local/include
    /usr/include
    NO_DEFAULT_PATH
)
find_library(REFINER_LIBRARIES
  NAMES RcapRefiner
  PATHS $ENV{REFINER_ROOT}/lib
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Refiner-1.1.04/lib/x86_64-linux
    $ENV{HOME}/REVOCAP_Refiner-1.1.04/lib/x86_64-linux
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    /usr/local/lib
    /usr/lib
	NO_DEFAULT_PATH
)

if(REFINER_INCLUDE_PATH AND REFINER_LIBRARIES)
set(REFINER_FOUND TRUE)
  set(HECMW_WITH_REFINER)
endif()
