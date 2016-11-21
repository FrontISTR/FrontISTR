#
# env REVOCAP_ROOT
#
# REVOCAP_FOUND
# REVOCAP_INCLUDE_PATH
# REVOCAP_LIBRARIES
#
find_path(REVOCAP_INCLUDE_PATH
  NAMES rcap.h
  PATHS $ENV{REVOCAP_ROOT}/librcap
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    /usr/local/include
    /usr/include
    NO_DEFAULT_PATH
)

find_library(REVOCAP_RCAP_LIBRARY
  NAMES rcap
  PATHS $ENV{REVOCAP_ROOT}/librcap
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    /usr/local/lib
    /usr/lib
    NO_DEFAULT_PATH
)

find_library(REVOCAP_RCAPF_LIBRARY
  NAMES rcapf
  PATHS $ENV{REVOCAP_ROOT}/librcap
    ${CMAKE_SOURCE_DIR}/../REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/REVOCAP_Coupler-2.1/librcap
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    /usr/local/lib
    /usr/lib
)

if(REVOCAP_INCLUDE_PATH AND REVOCAP_RCAP_LIBRARY AND REVOCAP_RCAPF_LIBRARY)
  set(REVOCAP_LIBRARIES ${REVOCAP_LIBRARIES}
    ${REVOCAP_RCAP_LIBRARY} ${REVOCAP_RCAPF_LIBRARY}
  )
  set(REVOCAP_FOUND TRUE)
else()
  set(WITH_REVOCAP OFF)
endif()
