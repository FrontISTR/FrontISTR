# Cluster PARDISO requires the BLACS library matching the actual MPI runtime.
# Set MKLBLACS_LIBRARY explicitly for an unrecognized MPI implementation or
# when cross-compiling (MPI_Get_library_version cannot be run on the host).
if(NOT MKLBLACS_LIBRARY)
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
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKLBLACS DEFAULT_MSG MKLBLACS_LIBRARY)
mark_as_advanced(MKLBLACS_LIBRARY)
unset(_mkl_mpi_version)
unset(_mkl_blacs_name)
unset(_mkl_library_dir)
