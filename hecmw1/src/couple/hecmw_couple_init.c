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
#include <assert.h>
#include <errno.h>

#include "hecmw_msgno.h"
#include "hecmw_config.h"
#include "hecmw_struct.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_control.h"
#include "hecmw_couple_info.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_bounding_box.h"
#include "hecmw_couple_background_cell.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_intra_iftable.h"
#include "hecmw_couple_inter_iftable.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_interpolate_info.h"
#include "hecmw_couple_init.h"



struct couple_info {
	char *boundary_id;						
	struct hecmw_couple_info *couple_info;	
	struct couple_info *next;				
} couple_list = {
	NULL,	/* boundary_id	*/
	NULL,	/* couple_info	*/
	NULL,	/* next			*/
};


/*================================================================================================*/

static void
free_couple_info(struct hecmw_couple_info *p)
{
	if(p == NULL) return;

	HECMW_couple_free_comm(p->comm_src);
	HECMW_couple_free_comm(p->comm_dst);
	HECMW_couple_free_comm(p->intercomm);
	HECMW_couple_free_boundary_info(p->boundary_src);
	HECMW_couple_free_boundary_info(p->boundary_dst);
	HECMW_couple_free_intra_iftable(p->intra_tbl_src);
	HECMW_couple_free_intra_iftable(p->intra_tbl_dst);
	HECMW_couple_free_inter_iftable(p->inter_tbl);
	HECMW_couple_free_mapped_point(p->mapped_point);
	HECMW_couple_free_weight_list(p->ip_list_pre);
	HECMW_couple_free_weight_list(p->ip_list_main);
	HECMW_couple_free_weight_list(p->ip_list_post);

	p->comm_src      = NULL;
	p->comm_dst      = NULL;
	p->intercomm     = NULL;
	p->boundary_src  = NULL;
	p->boundary_dst  = NULL;
	p->intra_tbl_src = NULL;
	p->intra_tbl_dst = NULL;
	p->inter_tbl     = NULL;
	p->mapped_point  = NULL;
	p->ip_list_pre   = NULL;
	p->ip_list_main  = NULL;
	p->ip_list_post  = NULL;
	HECMW_free(p);
	p = NULL;
}



static void
free_couple_list(struct couple_info *p)
{
	if(p == NULL) return;

	free_couple_info(p->couple_info);
	HECMW_free(p);
}



static struct hecmw_couple_info *
alloc_couple_info(void)
{
	struct hecmw_couple_info *p = NULL;

	p = (struct hecmw_couple_info *)HECMW_malloc(sizeof(struct hecmw_couple_info));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	p->unit_specifier_src = HECMW_COUPLE_UNIT_UNDEF;
	p->unit_specifier_dst = HECMW_COUPLE_UNIT_UNDEF;
	p->comm_src      = NULL;
	p->comm_dst      = NULL;
	p->intercomm     = NULL;
	p->boundary_src  = NULL;
	p->boundary_dst  = NULL;
	p->intra_tbl_src = NULL;
	p->intra_tbl_dst = NULL;
	p->inter_tbl     = NULL;
	p->mapped_point  = NULL;
	p->ip_list_pre   = NULL;
	p->ip_list_main  = NULL;
	p->ip_list_post  = NULL;

	return p;
}



static struct couple_info *
alloc_couple_list(void)
{
	struct couple_info *p = NULL;

	p = (struct couple_info *)HECMW_malloc(sizeof(struct couple_info));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	p->boundary_id = NULL;
	p->couple_info = NULL;
	p->next        = NULL;

	if((p->couple_info = alloc_couple_info()) == NULL) {
		HECMW_free(p);
		return NULL;
	}

	return p;
}



static struct couple_info *
cmp_couple_list(const char *boundary_id)
{
	struct couple_info *p;

	for(p=couple_list.next; p; p=p->next) {
		if(strcmp(p->boundary_id, boundary_id) == 0) return p;
	}

	return NULL;
}



static void
del_couple_list(const char *boundary_id)
{
	struct couple_info *p, *q;

	for(p=couple_list.next, q=&couple_list; p; p=p->next) {
		if(strcmp(p->boundary_id, boundary_id) == 0) {
			q->next = p->next;
			free_couple_list(p);
			break;
		}
	}
}



static int
add_couple_list(const char *boundary_id, int unit_specifier_src, int unit_specifier_dst,
		struct hecmw_couple_comm *comm_src,
		struct hecmw_couple_comm *comm_dst,
		struct hecmw_couple_comm *intercomm,
		struct hecmw_couple_boundary *boundary_src,
		struct hecmw_couple_boundary *boundary_dst,
		struct hecmw_couple_intra_iftable *intra_tbl_src,
		struct hecmw_couple_intra_iftable *intra_tbl_dst,
		struct hecmw_couple_inter_iftable *inter_tbl,
		struct hecmw_couple_mapped_point *mapped_point,
		struct hecmw_couple_weight_list *ip_list_pre,
		struct hecmw_couple_weight_list *ip_list_main,
		struct hecmw_couple_weight_list *ip_list_post)
{
	struct couple_info *p = NULL;

	if((p = cmp_couple_list(boundary_id)) != NULL) {
		HECMW_set_error(HECMWCPL_E, "");	/*@@*/
		return -1;
	}

	if((p = alloc_couple_list()) == NULL) return -1;

	p->boundary_id = HECMW_strdup(boundary_id);
	if(p->boundary_id == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	p->couple_info->unit_specifier_src = unit_specifier_src;
	p->couple_info->unit_specifier_dst = unit_specifier_dst;

	p->couple_info->comm_src      = comm_src;
	p->couple_info->comm_dst      = comm_dst;
	p->couple_info->intercomm     = intercomm;
	p->couple_info->boundary_src  = boundary_src;
	p->couple_info->boundary_dst  = boundary_dst;
	p->couple_info->intra_tbl_src = intra_tbl_src;
	p->couple_info->intra_tbl_dst = intra_tbl_dst;
	p->couple_info->inter_tbl     = inter_tbl;
	p->couple_info->mapped_point  = mapped_point;
	p->couple_info->ip_list_pre   = ip_list_pre;
	p->couple_info->ip_list_main  = ip_list_main;
	p->couple_info->ip_list_post  = ip_list_post;

	p->next          = couple_list.next;
	couple_list.next = p;

	return 0;

error:
	return -1;
}



extern struct hecmw_couple_info *
HECMW_couple_get_info(const char *boundary_id)
{
	struct couple_info *p;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple_get_info(): 'boundary_id' is NULL");
		return NULL;
	}

	if((p = cmp_couple_list(boundary_id)) == NULL) {
		HECMW_set_error(HECMWCPL_E, "");	/*@@*/
		return NULL;
	}

	return p->couple_info;
}


/*================================================================================================*/

extern void
HECMW_couple_free_init(const char *boundary_id)
{
	del_couple_list(boundary_id);
}



extern int
HECMW_couple_init(const char *boundary_id,
		struct hecmwST_local_mesh *mesh_unit1, struct hecmwST_local_mesh *mesh_unit2)
{
	struct hecmwST_local_mesh *mesh_src, *mesh_dst;
	struct hecmw_couple_boundary *boundary_src = NULL, *boundary_dst = NULL;
	struct hecmw_couple_bounding_box *bbox_src = NULL, *bbox_dst = NULL;
	struct hecmw_couple_background_cell *bgcell_src = NULL, *bgcell_dst = NULL;
	struct hecmw_couple_comm *comm_src = NULL, *comm_dst = NULL, *intercomm = NULL;
	struct hecmw_couple_intra_iftable *intra_tbl_src = NULL, *intra_tbl_dst = NULL;
	struct hecmw_couple_inter_iftable *inter_tbl = NULL;
	struct hecmw_couple_mapped_point *mapped_point_dst = NULL;
	struct hecmw_couple_weight_list *ip_list_pre = NULL;
	struct hecmw_couple_weight_list *ip_list_main = NULL;
	struct hecmw_couple_weight_list *ip_list_post = NULL;
	char src_unit_id[HECMW_NAME_LEN+1], dst_unit_id[HECMW_NAME_LEN+1];
	int is_unit1_memb, is_unit2_memb, is_src_memb, is_dst_memb;
	int unit_specifier_src, unit_specifier_dst;
	int direction;

	/* check argument */
	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple_init(): 'boundary_id' is NULL");
		return HECMW_ERROR;
	}

	is_unit1_memb = HECMW_couple_is_unit_member(boundary_id, HECMW_COUPLE_UNIT1);
	if(is_unit1_memb < 0) return HECMW_ERROR;
	is_unit2_memb = HECMW_couple_is_unit_member(boundary_id, HECMW_COUPLE_UNIT2);
	if(is_unit2_memb < 0) return HECMW_ERROR;

	if(is_unit1_memb && mesh_unit1 == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple_init(): 'mesh_unit1' is NULL");
		return HECMW_ERROR;
	}
	if(is_unit2_memb && mesh_unit2 == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple_init(): 'mesh_unit2' is NULL");
		return HECMW_ERROR;
	}

	/* set each couple unit */
	HECMW_couple_ctrl_get_direction(boundary_id, &direction);
	if(direction == HECMW_COUPLE_UNIT1_TO_UNIT2) {			/* UNIT1 -> UNIT2	*/
		unit_specifier_src = HECMW_COUPLE_UNIT1;
		unit_specifier_dst = HECMW_COUPLE_UNIT2;
		is_src_memb = is_unit1_memb;
		is_dst_memb = is_unit2_memb;
		mesh_src = mesh_unit1;
		mesh_dst = mesh_unit2;

	} else if(direction == HECMW_COUPLE_UNIT2_TO_UNIT1) {	/* UNIT2 -> UNIT1	*/
		unit_specifier_src = HECMW_COUPLE_UNIT2;
		unit_specifier_dst = HECMW_COUPLE_UNIT1;
		is_src_memb = is_unit2_memb;
		is_dst_memb = is_unit1_memb;
		mesh_src = mesh_unit2;
		mesh_dst = mesh_unit1;

	} else {												/* error			*/
		HECMW_set_error(HECMWCPL_E_INVALID_DIRECTION, "");
		goto error;
	}

	if(HECMW_couple_get_unit_id(boundary_id, unit_specifier_src,
				src_unit_id, HECMW_NAME_LEN+1) == NULL) goto error;
	if(HECMW_couple_get_unit_id(boundary_id, unit_specifier_dst,
				dst_unit_id, HECMW_NAME_LEN+1) == NULL) goto error;
	if((comm_src = HECMW_couple_get_intracomm(boundary_id, unit_specifier_src)) == NULL) goto error;
	if((comm_dst = HECMW_couple_get_intracomm(boundary_id, unit_specifier_dst)) == NULL) goto error;
	if((intercomm = HECMW_couple_get_intercomm(boundary_id)) == NULL) goto error;

	if(is_src_memb) {
		boundary_src = HECMW_couple_set_boundary_info(boundary_id, unit_specifier_src, mesh_src);
		if(boundary_src == NULL) goto error;

		bbox_src = HECMW_couple_set_bounding_box(boundary_id, mesh_src, boundary_src);
		if(bbox_src == NULL) goto error;

		bgcell_src = HECMW_couple_set_background_cell(boundary_id, mesh_src, bbox_src, boundary_src);
		if(bgcell_src == NULL) goto error;

		intra_tbl_src = HECMW_couple_make_intra_iftable(mesh_src, boundary_src, comm_src);
		if(intra_tbl_src == NULL) goto error;
	}

	if(is_dst_memb) {
		boundary_dst = HECMW_couple_set_boundary_info(boundary_id, unit_specifier_dst, mesh_dst);
		if(boundary_dst == NULL) goto error;

		bbox_dst = HECMW_couple_set_bounding_box(boundary_id, mesh_dst, boundary_dst);
		if(bbox_dst == NULL) goto error;

		bgcell_dst = HECMW_couple_set_background_cell(boundary_id, mesh_dst, bbox_dst, boundary_dst);
		if(bgcell_dst == NULL) goto error;

		mapped_point_dst = HECMW_couple_set_mapped_point(boundary_id, mesh_dst, boundary_dst);
		if(mapped_point_dst == NULL) goto error;

		intra_tbl_dst = HECMW_couple_make_intra_iftable(mesh_dst, boundary_dst, comm_dst);
		if(intra_tbl_dst == NULL) goto error;
	}

	/* make interface table for inter-communication */
	inter_tbl = HECMW_couple_set_map_data(mesh_src, mesh_dst, comm_src, comm_dst, intercomm,
			boundary_src, bbox_src, bbox_dst, bgcell_src, mapped_point_dst);
	if(inter_tbl == NULL) goto error;

	/* make weight list */
	if(comm_src->is_member) {
		ip_list_pre = HECMW_couple_make_pre_ip_list(mesh_src, boundary_src, comm_src, intra_tbl_src);
		if(ip_list_pre == NULL) goto error;
	}

	ip_list_main = HECMW_couple_make_main_ip_list(mesh_src, mesh_dst, boundary_src,
			boundary_dst, mapped_point_dst, comm_src, comm_dst, intercomm, inter_tbl);
	if(ip_list_main == NULL) goto error;

	if(comm_dst->is_member) {
		ip_list_post = HECMW_couple_make_post_ip_list(mesh_dst, boundary_dst,
				mapped_point_dst, comm_dst, intra_tbl_dst);
		if(ip_list_post == NULL) goto error;
	}

	if(add_couple_list(boundary_id, unit_specifier_src, unit_specifier_dst,
				comm_src, comm_dst, intercomm, boundary_src, boundary_dst,
				intra_tbl_src, intra_tbl_dst, inter_tbl, mapped_point_dst,
				ip_list_pre, ip_list_main, ip_list_post)) goto error;

	HECMW_couple_free_bounding_box(bbox_src);
	HECMW_couple_free_bounding_box(bbox_dst);
	HECMW_couple_free_background_cell(bgcell_src);
	HECMW_couple_free_background_cell(bgcell_dst);
	return HECMW_SUCCESS;

error:
	HECMW_couple_free_bounding_box(bbox_src);
	HECMW_couple_free_bounding_box(bbox_dst);
	HECMW_couple_free_background_cell(bgcell_src);
	HECMW_couple_free_background_cell(bgcell_dst);
	return HECMW_ERROR;
}

