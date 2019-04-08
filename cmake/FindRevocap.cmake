###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# REVOCAP_FOUND         TRUE if FindRevocap found REVOCAP_Coupler
# REVOCAP_INCLUDE_PATH  Include path of REVOCAP_Coupler
# REVOCAP_LIBRARIES     REVOCAP_Coupler libraries
#
# env REVOCAP_ROOT      Set REVOCAP_ROOT envionment variable,
#                       where REVOCAP_Coupler are.
#    ex. REVOCAP_ROOT=/home/someone/somewhere/REVOCAP_Coupler-2.1
#
if(REVOCAP_LIBRARIES)
  set(REVOCAP_FOUND TRUE)
  RETURN()
endif()

find_path(REVOCAP_INCLUDE_PATH
  NAMES rcap.h
  HINTS $ENV{REVOCAP_ROOT}/librcap
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    NO_DEFAULT_PATH
)

find_library(REVOCAP_RCAP_LIBRARY
  NAMES rcap
  HINTS $ENV{REVOCAP_ROOT}/librcap
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/lib
    /usr/lib
    NO_DEFAULT_PATH
)

find_library(REVOCAP_RCAPF_LIBRARY
  NAMES rcapf
  HINTS $ENV{REVOCAP_ROOT}/librcap
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/lib
    /usr/lib
)

if(REVOCAP_INCLUDE_PATH AND REVOCAP_RCAP_LIBRARY AND REVOCAP_RCAPF_LIBRARY)
  set(REVOCAP_LIBRARIES ${REVOCAP_LIBRARIES} ${REVOCAP_RCAP_LIBRARY} ${REVOCAP_RCAPF_LIBRARY})
  set(REVOCAP_FOUND TRUE)
endif()

mark_as_advanced(REVOCAP_INCLUDE_PATH REVOCAP_RCAP_LIBRARY REVOCAP_RCAPF_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(REVOCAP_Coupler
  DEFAULT_MSG REVOCAP_LIBRARIES REVOCAP_INCLUDE_PATH)

