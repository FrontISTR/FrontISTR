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

#install(TARGETS hecmw
#  EXPORT hecmw-export
#  LIBRARY       DESTINATION lib
#  INCLUDES      DESTINATION include
#  )
#install(EXPORT hecmw-export
#  FILE hecmw-config.cmake
#  DESTINATION share/cmake/hecmw/
#  EXPORT_LINK_INTERFACE_LIBRARIES
#  )
#export(TARGETS hecmw-export
#  FILE hecmw-config.cmake
#  EXPORT_LINK_INTERFACE_LIBRARIES
#  )
