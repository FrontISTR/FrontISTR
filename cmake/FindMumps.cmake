###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# MUMPS_FOUND         TRUE if FindMumps found mumps
# MUMPS_INCLUDE_PATH  Include path of mumps
# MUMPS_LIBRARIES     mumps libraries
#
# env MUMPS_ROOT      Set MUMPS_ROOT environment variable,
#                     where mumps are.
#    ex. export MUMPS_ROOT=/home/someone/somewhere/MUMPS_5.1.1
#
if(MUMPS_LIBRARIES)
  set(MUMPS_FOUND TRUE)
  RETURN()
endif()

if($ENV{MUMPS_ROOT})
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} $ENV{MUMPS_ROOT}/lib)
endif()
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.3.3/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.3.1/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.3.0/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.1/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.1/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.0/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.2/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.1/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.0/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.0.2/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.0.1/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} $ENV{HOME}/local/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} $ENV{HOME}/.local/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_LIBRARY_PATH})
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} /usr/local/lib)
set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} /usr/lib)

if(NOT WITH_MPI)
  if($ENV{MUMPS_ROOT})
    set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} $ENV{MUMPS_ROOT}/libseq)
  endif()
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.3.3/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.3.1/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.3.0/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.1/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.1/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.0/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.2/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.1/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.0/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.0.2/libseq)
  set(LIB_SEARCH_PATH ${LIB_SEARCH_PATH} ${CMAKE_SOURCE_DIR}/../MUMPS_5.0.1/libseq)
endif()

find_library(MUMPS_D_LIB
  NAMES dmumps
  HINTS ${LIB_SEARCH_PATH}
)
find_library(MUMPS_COMMON_LIB
  NAMES mumps_common
  HINTS ${LIB_SEARCH_PATH}
)
find_library(MUMPS_PORD_LIB
  NAMES pord
  HINTS ${LIB_SEARCH_PATH}
)

if(NOT WITH_MPI)
  find_library(MUMPS_MPISEQ_LIB
    NAMES mpiseq
    HINTS ${LIB_SEARCH_PATH}
  )
endif()

find_path(MUMPS_INCLUDE_PATH
  NAMES mumps_compat.h
  HINTS $ENV{MUMPS_ROOT}/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.1/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.2.0/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.2/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.1/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.1.0/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.0.2/include
    ${CMAKE_SOURCE_DIR}/../MUMPS_5.0.1/include
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/include
    /usr/include
)

if(WITH_MPI)
  if(MUMPS_D_LIB AND MUMPS_COMMON_LIB AND MUMPS_PORD_LIB AND MUMPS_INCLUDE_PATH)
    set(MUMPS_LIBRARIES ${MUMPS_D_LIB} ${MUMPS_COMMON_LIB} ${MUMPS_PORD_LIB})
    set(MUMPS_FOUND ON)
  endif()
  mark_as_advanced(MUMPS_INCLUDE_PATH MUMPS_LIBRARIES MUMPS_D_LIB MUMPS_COMMON_LIB MUMPS_PORD_LIB)
else()
  if(MUMPS_D_LIB AND MUMPS_COMMON_LIB AND MUMPS_PORD_LIB AND MUMPS_MPISEQ_LIB AND MUMPS_INCLUDE_PATH)
    set(MUMPS_LIBRARIES ${MUMPS_D_LIB} ${MUMPS_COMMON_LIB} ${MUMPS_PORD_LIB} ${MUMPS_MPISEQ_LIB})
    set(MUMPS_FOUND ON)
  endif()
  mark_as_advanced(MUMPS_INCLUDE_PATH MUMPS_LIBRARIES MUMPS_D_LIB MUMPS_COMMON_LIB MUMPS_PORD_LIB MUMPS_MPISEQ_LIB)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mumps
  DEFAULT_MSG MUMPS_LIBRARIES MUMPS_INCLUDE_PATH)
