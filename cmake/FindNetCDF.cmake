###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# NETCDF_FOUND         TRUE if FindNetCDF found netcdf
# NETCDF_INCLUDE_PATH  Include path of netcdf
# NETCDF_LIBRARIES     netcdf libraries
#
# env NETCDF_ROOT      Set NETCDF_ROOT environment variable,
#                      where netcdf is installed.
#    ex. export NETCDF_ROOT=/home/someone/somewhere/netcdf-4.9.2
#
if(NETCDF_LIBRARIES)
  set(NETCDF_FOUND TRUE)
  RETURN()
endif()

find_path(NETCDF_INCLUDE_PATH
  NAMES netcdf.h
  HINTS $ENV{NETCDF_ROOT}/include
    $ENV{NETCDF_DIR}/include
    $ENV{HOME}/local/include
    $ENV{HOME}/.local/include
    ${CMAKE_INCLUDE_PATH}
    /usr/local/netcdf/include
    /usr/local/include
    /usr/include
)

find_library(NETCDF_LIBRARIES
  NAMES netcdf
  HINTS $ENV{NETCDF_ROOT}/lib
    $ENV{NETCDF_ROOT}/lib64
    $ENV{NETCDF_DIR}/lib
    $ENV{NETCDF_DIR}/lib64
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    ${CMAKE_LIBRARY_PATH}
    /usr/local/netcdf/lib
    /usr/local/lib
    /usr/local/lib64
    /usr/lib
    /usr/lib64
)

if(NETCDF_INCLUDE_PATH AND NETCDF_LIBRARIES)
  set(NETCDF_FOUND TRUE)
endif()

mark_as_advanced(NETCDF_INCLUDE_PATH NETCDF_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDE_PATH)
