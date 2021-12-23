
add_subdirectory(11)
add_subdirectory(22)
add_subdirectory(33)
add_subdirectory(44)
add_subdirectory(66)
add_subdirectory(nn)

target_sources(hecmw
  PRIVATE
  hecmw_precond.f90
)
