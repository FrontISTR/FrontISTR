###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw_part1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_partlex.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_part_log.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mesh_hash_sort.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_mesh_edge_info.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_part_get_control.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_init_for_partition.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_graph.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_partition.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_partitioner.c
)
