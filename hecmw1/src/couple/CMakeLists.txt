###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

list(APPEND hecmw_SOURCES
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_comm.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_copy_c2f.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_copy_f2c.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_control.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_info.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_info_if.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_get_mesh.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_boundary_info.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_bounding_box.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_background_cell.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_mapped_point.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_judge.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_inter_iftable.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_intra_iftable.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_weight.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_s2n_average.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_s2n_with_area.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_s2n_dist_node.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_s2n_dist_surf.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_n2s_average.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_n2s_with_area.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_interpolate_info.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_init.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_startup.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_finalize.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_init_if.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_startup_if.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_if.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_finalize_if.c"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_define_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_struct_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_copy_c2f_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_copy_f2c_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_info_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_get_mesh_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_init_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_startup_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_cleanup_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_f.f90"
  "${CMAKE_CURRENT_LIST_DIR}/hecmw_couple_finalize_f.f90"
)

list(APPEND hecmw_INCLUDE_DIRS
  "${CMAKE_CURRENT_LIST_DIR}"
)
