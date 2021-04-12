###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(hecmw
	PRIVATE
    hecmw_transfer_result_c.f90
    hecmw_transfer_data_c.f90
    hecmw_transfer_data_c2f.f90
    hecmw_dist_result_copy_c2f.f90
    hecmw_dlb_data_transfer.c
    hecmw_dynamic_load_balancing.c
    hecmw_dlb_mesh2graph.c
    hecmw_dlb_read_control.c
    hecmw_dlb_migrate.c
    hecmw_dlb_elem_sr_adapt.c
    hecmw_dlb_comm_util.c
    output_ucd_c.c
    hecmw_dlb_mem_util.c
    hecmw_dlb_output_result.c
  PUBLIC
    hecmw_dlb_comm_util.h
    hecmw_repart.h
)
