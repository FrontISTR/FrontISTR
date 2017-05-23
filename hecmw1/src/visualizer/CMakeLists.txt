###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_surface_main.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_comm_util.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_mem_util.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_connectivity_build.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_surface_compute.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_case_table.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_calc_attr.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_intersection_find.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_patch_const.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_tetra_intersect.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_combine.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_rendering.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_define_parameters.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_font_texture.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_color_mapping.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_color_composite_sf.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_generate_histogram_sf.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_subimage_composite_sf.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_endian.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_fstr_output_femap.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_pvr_main.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_comm_util.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_mem_util.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_new_refine.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_ray_trace.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_define_parameters.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_font_texture.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_color_mapping.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_ucd_trans.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_voxel_gen.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_read_voxel.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_mesh_extent.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_color_composite_vr.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_generate_histogram_vr.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_subimage_composite_vr.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_vis_read_control.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_visualizer.c
)

list(APPEND hecmw_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)

