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




#ifndef INC_HECMW_COUPLE_INTERPOLATE_INFO
#define INC_HECMW_COUPLE_INTERPOLATE_INFO

#include "hecmw_struct.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_inter_iftable.h"

extern struct hecmw_couple_weight_list *
HECMW_couple_make_pre_ip_list(
		const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_intra_iftable *intra_tbl_src);

extern struct hecmw_couple_weight_list *
HECMW_couple_make_main_ip_list(
		const struct hecmwST_local_mesh *mesh_src, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_boundary *boundary_dst,
		const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		const struct hecmw_couple_inter_iftable *inter_tbl);

extern struct hecmw_couple_weight_list *
HECMW_couple_make_post_ip_list(
		const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst,
		const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_intra_iftable *intra_tbl_dst);

#endif	/* INC_HECMW_COUPLE_INTERPOLATE_INFO */
