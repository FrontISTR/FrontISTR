cmake_minimum_required(VERSION 2.8.11)

if(WITH_METIS)
	set(EXTRA_LIBS ${EXTRA_LIBS} ${METIS_LIBRARIES})
endif()

if(WITH_MPI)
	set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()


# hec2rcap
add_executable(hec2rcap
	hec2rcap/hec2rcap.c
)
target_link_libraries(hec2rcap hecmw ${EXTRA_LIBS})

# hecmw_part1
add_executable(hecmw_part1
  partitioner/hecmw_partlex.c
  partitioner/hecmw_part_log.c
  partitioner/hecmw_mesh_hash_sort.c
  partitioner/hecmw_mesh_edge_info.c
  partitioner/hecmw_part_get_control.c
  partitioner/hecmw_init_for_partition.c
  partitioner/hecmw_graph.c
  partitioner/hecmw_partition.c
  partitioner/hecmw_partitioner.c
)
target_link_libraries(hecmw_part1 hecmw ${EXTRA_LIBS})

# rmerge
add_executable(rmerge
	result_file_merger/fstr_rmerge.c
	result_file_merger/fstr_rmerge_util.c
)
target_link_libraries(rmerge hecmw ${EXTRA_LIBS})

# rconv
add_executable(rconv
  result_type_converter/hecmw_res_type_conv.c
)
target_link_libraries(rconv hecmw ${EXTRA_LIBS})

# hecmw_vis1
add_executable(hecmw_vis1
  visualizer/visualizer.c
)
target_link_libraries(hecmw_vis1 hecmw ${EXTRA_LIBS})
