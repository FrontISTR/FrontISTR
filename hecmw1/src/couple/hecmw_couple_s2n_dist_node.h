/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#ifndef INC_HECMW_COUPLE_S2N_DIST_NODE
#define INC_HECMW_COUPLE_S2N_DIST_NODE

#include "hecmw_struct.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_inter_iftable.h"

extern struct hecmw_couple_weight_list *
HECMW_couple_s2n_dist_node(
		const struct hecmwST_local_mesh *mesh_src, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_boundary *boundary_dst,
		const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_inter_iftable *inter_tbl);

#endif	/* INC_HECMW_COUPLE_S2N_DIST_NODE */
