###############################################################################
# Copyright (c) 2020 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

# Variables:
#
# MKL_FOUND         TRUE if FindMKL found MKL
# MKL_INCLUDE_PATH  Include path of MKL
# MKL_LIBRARIES     MKL libraries
#
# env MKL_ROOT      Set MKL_ROOT environment variable,
#
if(MKL_LIBRARIES)
  set(MKL_FOUND TRUE)
  RETURN()
endif()

find_path(MKL_INCLUDE_PATH
  NAMES mkl.h
  HINTS $ENV{MKLROOT}/include
  $ENV{HOME}/local/include
  $ENV{HOME}/.local/include
  ${CMAKE_INCLUDE_PATH}
  /opt/intel/mkl/include
  /usr/local/include/mkl
  /usr/include/mkl
  /usr/local/include
  /usr/include
)
find_library(_MKL_INTEL_LP64          NAMES mkl_intel_lp64           HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
find_library(_MKL_INTEL_THREAD        NAMES mkl_intel_thread         HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
find_library(_MKL_GNU_THREAD          NAMES mkl_gnu_thread           HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
find_library(_MKL_CORE                NAMES mkl_core                 HINTS $ENV{MKLROOT}/lib/intel64 $ENV{HOME}/local/lib $ENV{HOME}/.local/lib /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib )
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(MKL_LIBRARIES
    ${_MKL_INTEL_LP64}
    ${_MKL_INTEL_THREAD}
    ${_MKL_CORE}
    iomp5
    pthread
    m
    dl
    CACHE STRING "MKL for Intel")

elseif(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(MKL_LIBRARIES
    ${_MKL_INTEL_LP64}
    ${_MKL_GNU_THREAD}
    ${_MKL_CORE}
    gomp
    pthread
    m
    dl
    CACHE STRING "MKL for GCC")
endif()
# Cluster PARDISO needs the BLACS wrapper matching the MPI implementation.
# MKLBLACS_LIBRARY can be supplied explicitly when MPI cannot be identified.
if(WITH_MPI AND MKL_LIBRARIES)
  set(_mkl_mpi_version "${MPI_C_LIBRARY_VERSION_STRING};${MPI_Fortran_LIBRARY_VERSION_STRING}")
  if(_mkl_mpi_version MATCHES "Open MPI")
    set(_mkl_blacs_name mkl_blacs_openmpi_lp64)
  elseif(_mkl_mpi_version MATCHES "Intel.*MPI|MPICH")
    set(_mkl_blacs_name mkl_blacs_intelmpi_lp64)
  endif()
  if(_mkl_blacs_name)
    get_filename_component(_mkl_library_dir "${_MKL_CORE}" DIRECTORY)
    find_library(MKLBLACS_LIBRARY NAMES ${_mkl_blacs_name}
      HINTS "${_mkl_library_dir}" $ENV{MKLROOT}/lib/intel64 $ENV{MKLROOT}/lib
      /opt/intel/mkl/lib/intel64 /usr/lib/x86_64-linux-gnu)
  endif()
  if(MKLBLACS_LIBRARY)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      # Shared wrappers must survive --as-needed; static wrappers also need
      # extraction and export so MKL can find MKLMPI_Get_wrappers via dlsym.
      list(APPEND MKL_LIBRARIES "-Wl,--push-state,--no-as-needed")
      if(MKLBLACS_LIBRARY MATCHES "\\.a$")
        list(APPEND MKL_LIBRARIES "-Wl,--undefined=MKLMPI_Get_wrappers" "-Wl,--export-dynamic")
      endif()
      list(APPEND MKL_LIBRARIES ${MKLBLACS_LIBRARY} "-Wl,--pop-state")
    else()
      list(APPEND MKL_LIBRARIES ${MKLBLACS_LIBRARY})
    endif()
    set(MKL_LIBRARIES "${MKL_LIBRARIES}" CACHE STRING "MKL libraries including MPI BLACS" FORCE)
  endif()
  mark_as_advanced(MKLBLACS_LIBRARY)
  unset(_mkl_mpi_version)
  unset(_mkl_blacs_name)
  unset(_mkl_library_dir)
endif()

if(MKL_INCLUDE_PATH AND MKL_LIBRARIES)
  set(MKL_FOUND TRUE)
endif()

mark_as_advanced(MKL_INCLUDE_PATH MKL_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_PATH)
