###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPAND hecmw_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_transfer_result_c.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_transfer_data_c.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_transfer_data_c2f.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dist_result_copy_c2f.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_data_transfer.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dynamic_load_balancing.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_mesh2graph.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_read_control.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_migrate.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_elem_sr_adapt.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_comm_util.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/output_ucd_c.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_mem_util.c"
  "${CMAKE_CURRENT_SOURCE_DIR}/hecmw_dlb_output_result.c"
)
