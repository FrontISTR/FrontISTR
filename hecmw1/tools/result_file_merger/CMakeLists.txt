add_executable(rmerge)

target_sources(rmerge
  PRIVATE
  fstr_rmerge_util.c
  fstr_rmerge.c
)

target_link_libraries(rmerge
  PRIVATE
  hecmw
  OpenMP::OpenMP_C
  OpenMP::OpenMP_CXX
  OpenMP::OpenMP_Fortran
  MPI::MPI_C
  MPI::MPI_CXX
  MPI::MPI_Fortran
  )

target_include_directories(rmerge
  PRIVATE
  ${CMAKE_SOURCE_DIR}/hecmw1/src/common
  )

