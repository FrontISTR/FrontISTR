###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_executable(hecmw_part1)

target_sources(hecmw_part1
  PRIVATE
    hecmw_partlex.c
    hecmw_part_log.c
    hecmw_mesh_hash_sort.c
    hecmw_mesh_edge_info.c
    hecmw_part_get_control.c
    hecmw_init_for_partition.c
    hecmw_graph.c
    hecmw_partition.c
    hecmw_partitioner.c
  PUBLIC
    hecmw_partlex.h
    hecmw_part_log.h
    hecmw_mesh_hash_sort.h
    hecmw_mesh_edge_info.h
    hecmw_part_get_control.h
    hecmw_init_for_partition.h
    hecmw_graph.h
    hecmw_partition.h
  )

target_link_libraries(hecmw_part1
  PRIVATE
  hecmw
  )

target_include_directories(hecmw_part1
  PRIVATE
    ../../src/common
  )

install(TARGETS hecmw_part1
  DESTINATION bin
  )
