/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_DIST_PRINT_INCLUDED
#define HECMW_DIST_PRINT_INCLUDED

#include <stdio.h>
#include "hecmw_struct.h"


extern void HECMW_dist_print_flags(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_header(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_gridfile(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_files(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_zero_temp(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_node(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_elem(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_pe(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_adapt(const struct hecmwST_local_mesh *mesh, FILE *fp);


extern void HECMW_dist_print_section(const struct hecmwST_section *sect, FILE *fp);


extern void HECMW_dist_print_material(const struct hecmwST_material *material, FILE *fp);


extern void HECMW_dist_print_mpc(const struct hecmwST_mpc *mpc, FILE *fp);


extern void HECMW_dist_print_amp(const struct hecmwST_amplitude *amp, FILE *fp);


extern void HECMW_dist_print_ngrp(const struct hecmwST_node_grp *ngrp, FILE *fp);


extern void HECMW_dist_print_egrp(const struct hecmwST_elem_grp *egrp, FILE *fp);


extern void HECMW_dist_print_sgrp(const struct hecmwST_surf_grp *sgrp, FILE *fp);


extern void HECMW_dist_print(const struct hecmwST_local_mesh *mesh, FILE *fp);

#endif

