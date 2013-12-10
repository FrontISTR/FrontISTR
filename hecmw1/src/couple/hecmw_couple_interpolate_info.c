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




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "hecmw_struct.h"
#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_inter_iftable.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_s2n_average.h"
#include "hecmw_couple_s2n_with_area.h"
#include "hecmw_couple_s2n_dist_node.h"
#include "hecmw_couple_s2n_dist_surf.h"
#include "hecmw_couple_n2s_average.h"
#include "hecmw_couple_n2s_with_area.h"
#include "hecmw_couple_interpolate_info.h"



extern struct hecmw_couple_weight_list *
HECMW_couple_make_pre_ip_list(
		const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_intra_iftable *intra_tbl_src)
{
	struct hecmw_couple_weight_list *weight_list = NULL, *p;

	if(mesh_src == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_pre_ip_list(): 'mesh_src' is NULL");
		return NULL;
	}
	if(boundary_src == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_pre_ip_list(): 'boundary_src' is NULL");
		return NULL;
	}
	if(comm_src == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_pre_ip_list(): 'comm_src' is NULL");
		return NULL;
	}

	if((weight_list = HECMW_couple_alloc_weight_list()) == NULL) goto error;
	p = weight_list;

	if(boundary_src->data_type == HECMW_COUPLE_NODE_GROUP) {			/* node group		*/
		p->next = NULL;

	} else if(boundary_src->data_type == HECMW_COUPLE_SURFACE_GROUP) {	/* surface group	*/
		p->next = HECMW_couple_s2n_average(mesh_src, boundary_src);
/*		p->next = HECMW_couple_s2n_with_area(mesh_src, boundary_src);	*/
		if(p->next == NULL) goto error;

#if 0
	} else if(boundary_src->data_type == HECMW_COUPLE_ELEMENT_GROUP) {	/* element group	*/
		p->next = HECMW_couple_e2n_by_average(mesh_src, boundary_src);
		if(p->next == NULL) goto error;
#endif

	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_DATATYPE, "");
		goto error;
	}

	return weight_list;

error:
	return NULL;
}


extern struct hecmw_couple_weight_list *
HECMW_couple_make_main_ip_list(
		const struct hecmwST_local_mesh *mesh_src, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_boundary *boundary_dst,
		const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		const struct hecmw_couple_inter_iftable *inter_tbl)
{
	struct hecmw_couple_weight_list *weight_list = NULL, *p;

	if(comm_src == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_main_ip_list(): 'comm_src' is NULL");
		goto error;
	}
	if(comm_dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_main_ip_list(): 'comm_dst' is NULL");
		goto error;
	}
	if(intercomm == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_main_ip_list(): 'intercomm' is NULL");
		goto error;
	}
	if(inter_tbl == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_main_ip_list(): 'inter_tbl' is NULL");
		goto error;
	}
	if(comm_src->is_member) {
		if(mesh_src == NULL) {
			HECMW_set_error(HECMWCPL_E_INVALID_ARG,
					"HECMW_couple_make_main_ip_list(): 'mesh_src' is NULL");
			goto error;
		}
		if(boundary_src == NULL) {
			HECMW_set_error(HECMWCPL_E_INVALID_ARG,
					"HECMW_couple_make_main_ip_list(): 'boundary_src' is NULL");
			goto error;
		}
	}
	if(comm_dst->is_member) {
		if(mesh_dst == NULL) {
			HECMW_set_error(HECMWCPL_E_INVALID_ARG,
					"HECMW_couple_make_main_ip_list(): 'mesh_dst' is NULL");
			goto error;
		}
		if(boundary_dst == NULL) {
			HECMW_set_error(HECMWCPL_E_INVALID_ARG,
					"HECMW_couple_make_main_ip_list(): 'boundary_dst' is NULL");
			goto error;
		}
		if(mapped_point == NULL) {
			HECMW_set_error(HECMWCPL_E_INVALID_ARG,
					"HECMW_couple_make_main_ip_list(): 'mapped_point' is NULL");
			goto error;
		}
	}

	if((weight_list = HECMW_couple_alloc_weight_list()) == NULL) goto error;
	p = weight_list;

	if(boundary_src->data_type == HECMW_COUPLE_NODE_GROUP) {
		p->next = HECMW_couple_s2n_dist_node(mesh_src, mesh_dst, comm_src, comm_dst, intercomm,
				boundary_src, boundary_dst, mapped_point, inter_tbl);
		if(p->next == NULL) goto error;

	} else if(boundary_src->data_type == HECMW_COUPLE_SURFACE_GROUP) {
		p->next = HECMW_couple_s2n_dist_surf(mesh_src, mesh_dst, comm_src, comm_dst, intercomm,
				boundary_src, boundary_dst, mapped_point, inter_tbl);
		if(p->next == NULL) goto error;

#if 0
	} else if(boundary_src->data_type == HECMW_COUPLE_ELEMENT_GROUP) {
#endif

	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_DATATYPE, "");
		goto error;
	}

	return weight_list;

error:
	return NULL;
}


extern struct hecmw_couple_weight_list *
HECMW_couple_make_post_ip_list(
		const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst,
		const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_intra_iftable *intra_tbl_dst)
{
	struct hecmw_couple_weight_list *weight_list = NULL, *p;

	if(mesh_dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_post_ip_list(): 'mesh_dst' is NULL");
		return NULL;
	}
	if(boundary_dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_post_ip_list(): 'boundary_dst' is NULL");
		return NULL;
	}
	if(mapped_point == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_post_ip_list(): 'mapped_point' is NULL");
		return NULL;
	}
	if(comm_dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_post_ip_list(): 'comm_dst' is NULL");
		return NULL;
	}

	if((weight_list = HECMW_couple_alloc_weight_list()) == NULL) goto error;
	p = weight_list;

	if(boundary_dst->data_type == HECMW_COUPLE_NODE_GROUP) {			/* node group		*/
		p->next = NULL;

	} else if(boundary_dst->data_type == HECMW_COUPLE_SURFACE_GROUP) {	/* surface group	*/
		p->next = HECMW_couple_n2s_average(mesh_dst, boundary_dst, comm_dst, intra_tbl_dst);
/*		p->next = HECMW_couple_n2s_with_area(mesh_dst, boundary_dst);	*/
		if(p->next == NULL) goto error;

#if 0
	} else if(boundary_dst->data_type == HECMW_COUPLE_ELEMENT_GROUP) {	/* element group */
#endif

	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_DATATYPE, "");
		goto error;
	}

	return weight_list;

error:
	return NULL;
}
