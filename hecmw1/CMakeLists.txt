add_library(hecmw)
target_include_directories(hecmw
  PRIVATE
  src/common
  src/hecmw
  src/visualizer
  )
target_link_libraries(hecmw
  PRIVATE
  OpenMP::OpenMP_C
  OpenMP::OpenMP_CXX
  OpenMP::OpenMP_Fortran
  MPI::MPI_C
  MPI::MPI_CXX
  MPI::MPI_Fortran
  )
add_subdirectory(src)

# tools
add_subdirectory(tools)
