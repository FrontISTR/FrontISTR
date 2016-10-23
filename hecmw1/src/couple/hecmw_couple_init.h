/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_INIT
#define INC_HECMW_COUPLE_INIT

#include "hecmw_struct.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_intra_iftable.h"
#include "hecmw_couple_inter_iftable.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_weight.h"

struct hecmw_couple_info {
	int unit_specifier_src;
	int unit_specifier_dst;
	struct hecmw_couple_comm *comm_src;
	struct hecmw_couple_comm *comm_dst;
	struct hecmw_couple_comm *intercomm;
	struct hecmw_couple_boundary *boundary_src;
	struct hecmw_couple_boundary *boundary_dst;
	struct hecmw_couple_intra_iftable *intra_tbl_src;
	struct hecmw_couple_intra_iftable *intra_tbl_dst;
	struct hecmw_couple_inter_iftable *inter_tbl;
	struct hecmw_couple_mapped_point *mapped_point;
	struct hecmw_couple_weight_list *ip_list_pre;
	struct hecmw_couple_weight_list *ip_list_main;
	struct hecmw_couple_weight_list *ip_list_post;
};

extern struct hecmw_couple_info *
HECMW_couple_get_info(const char *boundary_id);

extern int
HECMW_couple_init(const char *couple_name,
		struct hecmwST_local_mesh *unit1_mesh, struct hecmwST_local_mesh *unit2_mesh);

extern void
HECMW_couple_free_init(const char *boundary_id);

#endif	/* INC_HECMW_COUPLE_INIT */
