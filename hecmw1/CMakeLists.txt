###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_library(hecmw)

add_subdirectory(src)

add_subdirectory(tools)

target_link_libraries(hecmw
  PRIVATE
  MPI::MPI_C
  MPI::MPI_CXX
  MPI::MPI_Fortran
  OpenMP::OpenMP_C
  OpenMP::OpenMP_CXX
  OpenMP::OpenMP_Fortran
  #  METIS::METIS
  )
target_include_directories(hecmw
  PRIVATE
  src/visualizer
  src/common
  )
target_compile_definitions(hecmw
  PUBLIC 
  )
