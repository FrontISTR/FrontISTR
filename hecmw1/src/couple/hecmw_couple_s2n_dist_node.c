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

#include "hecmw_common_define.h"
#include "hecmw_struct.h"
#include "hecmw_malloc.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_mapped_point.h"
#include "hecmw_couple_inter_iftable.h"
#include "hecmw_couple_s2n_dist_node.h"


#define FRAC_1_2 (0.5)


#define FRAC_1_3 (0.33333333333333333)


#define FRAC_1_4 (0.25)


#define EPS_ZERO (1.0E-24)


/*================================================================================================*/

struct link_list {
	int id;						
	double weight;				
	struct link_list *next;		
};


struct hecmw_couple_vertex {
	double x;	
	double y;	
	double z;	
};


struct hecmw_couple_vector {
	double x;	
	double y;	
	double z;	
};


/*================================================================================================*/

static void
free_link_list(struct link_list *r)
{
	struct link_list *p, *q;

	for(p=r; p; p=q) {
		q = p->next;
		HECMW_free(p);
	}
}



static int
intercomm_d2s_coord(const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, double **coord)
{
	int *sendbuf_index = NULL, *recvbuf_index = NULL;
	double *sendbuf = NULL;
	int size, rtc, i;

	if(comm_dst->is_member) {
		sendbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_import+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<=inter_tbl->n_neighbor_pe_import; i++) {
			sendbuf_index[i] = 3 * inter_tbl->import_index[i];
		}

		size = sizeof(double) * (inter_tbl->import_index[inter_tbl->n_neighbor_pe_import] * 3 + 1);
		sendbuf = (double *)HECMW_malloc(size);
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<inter_tbl->import_index[inter_tbl->n_neighbor_pe_import]; i++) {
			sendbuf[3*i]   = mapped_point->coord[3*(inter_tbl->import_item[i])];
			sendbuf[3*i+1] = mapped_point->coord[3*(inter_tbl->import_item[i])+1];
			sendbuf[3*i+2] = mapped_point->coord[3*(inter_tbl->import_item[i])+2];
		}
	}

	if(comm_src->is_member) {
		recvbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_export+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<=inter_tbl->n_neighbor_pe_export; i++) {
			recvbuf_index[i] = 3 * inter_tbl->export_index[i];
		}

		size = sizeof(double) * (inter_tbl->export_index[inter_tbl->n_neighbor_pe_export] * 3 + 1);
		*coord = (double *)HECMW_malloc(size);
		if(*coord == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	rtc = HECMW_couple_inter_send_recv(
			inter_tbl->n_neighbor_pe_import, inter_tbl->neighbor_pe_import, sendbuf_index, sendbuf,
			inter_tbl->n_neighbor_pe_export, inter_tbl->neighbor_pe_export, recvbuf_index, *coord,
			HECMW_DOUBLE, intercomm->comm);
	if(rtc != 0) goto error;

	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recvbuf_index);
	return 0;

error:
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recvbuf_index);
	return -1;
}



static int
s2n_dist_node_tet1(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src, int id,
		const struct hecmw_couple_vertex *coord_dst,
		struct link_list *weight_list)
{
	struct link_list *p;
	struct hecmw_couple_vertex coord[3], gravity;
	double d, r_d_surf, r_d_node[3], r_d_sum = 0.0;
	int node_id[3], node, n, i;

	for(n=0, i=boundary_src->elem_node_index[id]; i<boundary_src->elem_node_index[id+1]; i++) {
		node_id[n] = boundary_src->elem_node_item[i];
		node       = boundary_src->node->item[node_id[n]];
		coord[n].x = mesh_src->node[3*(node-1)];
		coord[n].y = mesh_src->node[3*(node-1)+1];
		coord[n].z = mesh_src->node[3*(node-1)+2];
		n++;
	}

	for(i=0; i<3; i++) {
		d = sqrt((coord[i].x-coord_dst->x)*(coord[i].x-coord_dst->x) +
				(coord[i].y-coord_dst->y)*(coord[i].y-coord_dst->y) +
				(coord[i].z-coord_dst->z)*(coord[i].z-coord_dst->z));
		r_d_node[i] = 1.0 / (d+EPS_ZERO);
		r_d_sum    += r_d_node[i];
	}

	for(i=0; i<3; i++) {
		p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
		if(p == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		p->id     = node_id[i];
		p->weight = r_d_node[i] / r_d_sum;
		p->next   = weight_list->next;
		weight_list->next = p;
	}

	return 0;
}



static int
s2n_dist_node_hex1(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src, int id,
		const struct hecmw_couple_vertex *coord_dst,
		struct link_list *weight_list)
{
	struct link_list *p;
	struct hecmw_couple_vertex coord[4], gravity;
	double d, r_d_surf, r_d_node[4], r_d_sum = 0.0;
	int node_id[4], node, n, i;

	for(n=0, i=boundary_src->elem_node_index[id]; i<boundary_src->elem_node_index[id+1]; i++) {
		node_id[n] = boundary_src->elem_node_item[i];
		node       = boundary_src->node->item[node_id[n]];
		coord[n].x = mesh_src->node[3*(node-1)];
		coord[n].y = mesh_src->node[3*(node-1)+1];
		coord[n].z = mesh_src->node[3*(node-1)+2];
		n++;
	}

	for(i=0; i<4; i++) {
		d = sqrt((coord[i].x-coord_dst->x)*(coord[i].x-coord_dst->x) +
				(coord[i].y-coord_dst->y)*(coord[i].y-coord_dst->y) +
				(coord[i].z-coord_dst->z)*(coord[i].z-coord_dst->z));
		r_d_node[i] = 1.0 / (d+EPS_ZERO);
		r_d_sum += r_d_node[i];
	}

	for(i=0; i<4; i++) {
		p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
		if(p == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		p->id     = node_id[i];
		p->weight = r_d_node[i] / r_d_sum;
		p->next   = weight_list->next;
		weight_list->next = p;
	}

	return 0;
}



static int
s2n_dist_node(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_inter_iftable *inter_tbl, double *coord,
		struct hecmw_couple_weight *weight_info)
{
	struct link_list *weight_list = NULL, *p;
	struct hecmw_couple_vertex coord_dst;
	int elem, n_item, id, rtc, n, i;

	n_item = inter_tbl->export_index[inter_tbl->n_neighbor_pe_export] + 1;
	weight_list = (struct link_list *)HECMW_malloc(sizeof(struct link_list)*n_item);
	if(weight_list == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<n_item; i++) {
		weight_list[i].id     = -1;
		weight_list[i].weight = 0.0;
		weight_list[i].next   = NULL;
	}

	for(i=0; i<inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]; i++) {
		coord_dst.x = coord[3*i];
		coord_dst.y = coord[3*i+1];
		coord_dst.z = coord[3*i+2];

		id = inter_tbl->export_item[i];
		elem = boundary_src->surf->item[2*id];
		if(mesh_src->elem_type[elem-1] == HECMW_ETYPE_TET1) {
			rtc = s2n_dist_node_tet1(mesh_src, boundary_src, id, &coord_dst, &weight_list[i]);
			if(rtc != HECMW_SUCCESS) goto error;
		} else if(mesh_src->elem_type[elem-1] == HECMW_ETYPE_HEX1) {
			rtc = s2n_dist_node_hex1(mesh_src, boundary_src, id, &coord_dst, &weight_list[i]);
			if(rtc != HECMW_SUCCESS) goto error;
		} else {
			HECMW_set_error(HECMWCPL_E_NONSUPPORT_ETYPE, "");
			goto error;
		}
	}


	weight_info->n = inter_tbl->export_index[inter_tbl->n_neighbor_pe_export];
	weight_info->type = HECMW_COUPLE_IP_NODE_TO_NODE;

	weight_info->index = (int *)HECMW_calloc(weight_info->n+1, sizeof(int));
	if(weight_info->index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]; i++) {
		for(p=weight_list[i].next; p; p=p->next) {
			n++;
		}
		weight_info->index[i+1] = n;
	}

	n_item = weight_info->index[weight_info->n];
	weight_info->id = (int *)HECMW_malloc(sizeof(int)*n_item);
	if(weight_info->id == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	weight_info->weight = (double *)HECMW_malloc(sizeof(double)*n_item);
	if(weight_info->weight == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]; i++) {
		for(p=weight_list[i].next; p; p=p->next) {
			weight_info->id[n]     = p->id;
			weight_info->weight[n] = p->weight;
			n++;
		}
	}

	for(i=0; i<inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]; i++) {
		free_link_list(weight_list[i].next);
	}
	HECMW_free(weight_list);

	return 0;

error:
	for(i=0; i<inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]; i++) {
		free_link_list(weight_list[i].next);
	}
	HECMW_free(weight_list);

	return -1;
}



extern struct hecmw_couple_weight_list *
HECMW_couple_s2n_dist_node(
		const struct hecmwST_local_mesh *mesh_src, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_boundary *boundary_dst,
		const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_inter_iftable *inter_tbl)
{
	struct hecmw_couple_weight_list *weight_info_list = NULL;
	struct hecmw_couple_weight *weight_info = NULL;
	double *coord = NULL;
	int rtc, i;

	if((weight_info_list = HECMW_couple_alloc_weight_list()) == NULL) return NULL;

	rtc = intercomm_d2s_coord(mapped_point, inter_tbl, comm_src, comm_dst, intercomm, &coord);
	if(rtc) goto error;

	if(comm_src->is_member) {
		if((weight_info = HECMW_couple_alloc_weight()) == NULL) goto error;
		weight_info_list->info = weight_info;

		if(s2n_dist_node(mesh_src, boundary_src, inter_tbl, coord, weight_info)) goto error;
	}

	return weight_info_list;

error:
	return NULL;
}

