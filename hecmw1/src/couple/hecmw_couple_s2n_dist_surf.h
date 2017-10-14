/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_S2N_DIST_SURF
#define INC_HECMW_COUPLE_S2N_DIST_SURF

#include "hecmw_struct.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_inter_iftable.h"

extern struct hecmw_couple_weight_list *HECMW_couple_s2n_dist_surf(
    const struct hecmwST_local_mesh *mesh_src,
    const struct hecmwST_local_mesh *mesh_dst,
    const struct hecmw_couple_comm *comm_src,
    const struct hecmw_couple_comm *comm_dst,
    const struct hecmw_couple_comm *intercomm,
    const struct hecmw_couple_boundary *boundary_src,
    const struct hecmw_couple_boundary *boundary_dst,
    const struct hecmw_couple_mapped_point *mapped_point,
    const struct hecmw_couple_inter_iftable *inter_tbl);

#endif /* INC_HECMW_COUPLE_S2N_DIST_SURF */
