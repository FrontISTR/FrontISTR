###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_int_sr.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_int_sr_rev.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_item_sr.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_real_sr.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_stack_sr.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_active.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_adjemb.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_bc_pointer.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_cell_comm_table.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_edge_comm_table.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_edge_info.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_error_exit.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_extemb.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_get_new_cell_info.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_grid_smooth.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_allocate.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_new_cell_341.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_new_cell_351.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_new_cell.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_new_node.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_new_mesh.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_repro_comm_table.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_init.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_adapt_proc.f90
)
