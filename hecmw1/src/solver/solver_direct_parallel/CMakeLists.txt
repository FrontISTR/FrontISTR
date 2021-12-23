
target_sources(hecmw
  PRIVATE
  matrix2graph.c
  matrix_repart.c
  mem_util.c
  separator_c2f_c.c
  m_irjc_matrix.f90
  m_crs_matrix.f90
  m_child_matrix.f90
  m_matrix_partition_info.f90
  m_elap.F90
  hecmw_solver_direct_parallel.F90
)
