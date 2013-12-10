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

#include "hecmw_struct.h"
#include "hecmw_msgno.h"
#include "hecmw_common_define.h"
#include "hecmw_error.h"
#include "hecmw_etype.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_mapped_point.h"

/*================================================================================================*/

extern struct hecmw_couple_mapped_point *
HECMW_couple_alloc_mapped_point(void)
{
	struct hecmw_couple_mapped_point *p = NULL;

	p = (struct hecmw_couple_mapped_point *)HECMW_malloc(sizeof(struct hecmw_couple_mapped_point));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	p->n     = 0;
	p->type  = HECMW_COUPLE_MAP_UNDEF;
	p->item  = NULL;
	p->id    = NULL;
	p->coord = NULL;

	return p;
}


extern void
HECMW_couple_free_mapped_point(struct hecmw_couple_mapped_point *p)
{
	if(p == NULL) return;

	HECMW_free(p->item);
	HECMW_free(p->id);
	HECMW_free(p->coord);
	HECMW_free(p);
	p = NULL;
}

/*================================================================================================*/

static struct hecmw_couple_mapped_point *
set_mapped_point_by_node(const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst)
{
	struct hecmw_couple_mapped_point *mapped_point = NULL;
	int node, i;

	mapped_point = HECMW_couple_alloc_mapped_point();
	if(mapped_point == NULL) return NULL;

	mapped_point->n    = boundary_dst->node->n;
	mapped_point->type = HECMW_COUPLE_NODE_GROUP;

	mapped_point->item = (int *)HECMW_malloc(sizeof(int)*mapped_point->n);
	if(mapped_point->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	mapped_point->id = (int *)HECMW_malloc(sizeof(int)*mapped_point->n);
	if(mapped_point->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	mapped_point->coord = (double *)HECMW_malloc(sizeof(double)*mapped_point->n*3);
	if(mapped_point->coord == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	for(i=0; i<boundary_dst->node->n; i++) {
		node = boundary_dst->node->item[i];
		mapped_point->item[i] = node;
		mapped_point->id[i]   = i;
		mapped_point->coord[3*i  ] = mesh_dst->node[3*(node-1)  ];
		mapped_point->coord[3*i+1] = mesh_dst->node[3*(node-1)+1];
		mapped_point->coord[3*i+2] = mesh_dst->node[3*(node-1)+2];
	}

	return mapped_point;

error:
	HECMW_couple_free_mapped_point(mapped_point);
	return NULL;
}


static struct hecmw_couple_mapped_point *
set_mapped_point_by_elem(const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst)
{
	struct hecmw_couple_mapped_point *mapped_point = NULL;
	double coord_x_sum, coord_y_sum, coord_z_sum;
	int elem, node, max_node, i, j;

	mapped_point = HECMW_couple_alloc_mapped_point();
	if(mapped_point == NULL) return NULL;

	mapped_point->n    = boundary_dst->elem->n;
	mapped_point->type = HECMW_COUPLE_ELEMENT_GROUP;

	mapped_point->item = (int *)HECMW_malloc(sizeof(int)*mapped_point->n);
	if(mapped_point->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	mapped_point->id = (int *)HECMW_malloc(sizeof(int)*mapped_point->n);
	if(mapped_point->id == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	mapped_point->coord = (double *)HECMW_malloc(sizeof(double)*mapped_point->n);
	if(mapped_point->coord == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	for(i=0; i<mapped_point->n; i++) {
		coord_x_sum = coord_y_sum = coord_z_sum = 0.0;
		elem     = boundary_dst->elem->item[i];
		max_node = HECMW_get_max_node(mesh_dst->elem_type[elem-1]);
		for(j=mesh_dst->elem_node_index[elem-1]; j<mesh_dst->elem_node_index[elem]; j++) {
			node         = mesh_dst->elem_node_item[j];
			coord_x_sum += mesh_dst->node[3*(node-1)  ];
			coord_y_sum += mesh_dst->node[3*(node-1)+1];
			coord_z_sum += mesh_dst->node[3*(node-1)+2];
		}
		mapped_point->item[i]      = elem;
		mapped_point->id[i]        = i;
		mapped_point->coord[3*i]   = coord_x_sum / max_node;
		mapped_point->coord[3*i+1] = coord_x_sum / max_node;
		mapped_point->coord[3*i+2] = coord_x_sum / max_node;
	}

	return mapped_point;

error:
	HECMW_couple_free_mapped_point(mapped_point);
	return NULL;
}


static struct hecmw_couple_mapped_point *
set_mapped_point_by_surf(const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst)
{
	struct hecmw_couple_mapped_point *mapped_point = NULL;

	/*@@@ NOT implementing @@@*/
	goto error;

	return mapped_point;

error:
	HECMW_couple_free_mapped_point(mapped_point);
	return NULL;
}


extern struct hecmw_couple_mapped_point *
HECMW_couple_set_mapped_point(const char *boundary_id, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst)
{
	struct hecmw_couple_mapped_point *mapped_point = NULL;
	int map_type;

	if(mesh_dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (mesh_dst)");
		return NULL;
	}
	if(boundary_dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_dst)");
		return NULL;
	}

	map_type = HECMW_COUPLE_MAP_NODE_TO_SURF;

	if(map_type == HECMW_COUPLE_MAP_NODE_TO_SURF) {
		if((mapped_point = set_mapped_point_by_node(mesh_dst, boundary_dst)) == NULL) goto error;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_MAPTYPE, "");
		goto error;
	}

	return mapped_point;

error:
	HECMW_couple_free_mapped_point(mapped_point);
	return NULL;
}

