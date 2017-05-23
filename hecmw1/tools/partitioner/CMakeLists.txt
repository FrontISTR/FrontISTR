###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_part1_SOURCES
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
