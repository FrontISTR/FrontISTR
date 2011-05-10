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
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_s2n_average.h"


#define FRAC_1_3 (0.33333333333333333)


#define FRAC_1_4 (0.25)



struct link_list {
	int id;						
	double weight;				
	struct link_list *next;		
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
s2n_average_tet1(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary, int id, struct link_list *weight_list)
{
	struct link_list *p;
	int node_id, n, i;

	for(n=0, i=boundary->elem_node_index[id]; i<boundary->elem_node_index[id+1]; i++) {
		node_id = boundary->elem_node_item[i];

		p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
		if(p == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		p->id     = id;
		p->weight = FRAC_1_3;
		p->next   = weight_list[node_id].next;
		weight_list[node_id].next = p;
	}

	return 0;
}



static int
s2n_average_hex1(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary, int id, struct link_list *weight_list)
{
	struct link_list *p;
	int node_id, n, i;

	for(n=0, i=boundary->elem_node_index[id]; i<boundary->elem_node_index[id+1]; i++) {
		node_id = boundary->elem_node_item[i];

		p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
		if(p == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		p->id     = id;
		p->weight = FRAC_1_4;
		p->next   = weight_list[node_id].next;
		weight_list[node_id].next = p;
	}

	return 0;
}



static int
s2n_average(const struct hecmwST_local_mesh *mesh, const struct hecmw_couple_boundary *boundary,
		struct hecmw_couple_weight *weight_info)
{
	struct link_list *weight_list = NULL, *p;
	int elem, n_item, size, n, i;

	size = sizeof(struct link_list)*boundary->node->n;
	weight_list = (struct link_list *)HECMW_malloc(size);
	if(weight_list == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<boundary->node->n; i++) {
		weight_list[i].id     = -1;
		weight_list[i].weight = 0.0;
		weight_list[i].next   = NULL;
	}

	/*
	 * calculate weight
	 */
	for(i=0; i<boundary->surf->n; i++) {
		elem = boundary->surf->item[2*i];

		if(mesh->elem_type[elem-1] == HECMW_ETYPE_TET1) {
			if(s2n_average_tet1(mesh, boundary, i, weight_list)) goto error;
		} else if(mesh->elem_type[elem-1] == HECMW_ETYPE_HEX1) {
			if(s2n_average_hex1(mesh, boundary, i, weight_list)) goto error;
		} else {
			HECMW_set_error(HECMWCPL_E_NONSUPPORT_ETYPE, "");
			goto error;
		}
	}

	/*
	 * make interpolating information
	 */
	/* number of nodes */
	weight_info->n = boundary->node->n;

	/* interpolating type */
	weight_info->type = HECMW_COUPLE_IP_SURF_TO_NODE;

	/* index of list */
	weight_info->index = (int *)HECMW_calloc(weight_info->n+1, sizeof(int));
	if(weight_info->index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<boundary->node->n; i++) {
		for(p=weight_list[i].next; p; p=p->next) {
			n++;
		}
		weight_info->index[i+1] = n;
	}

	/* id which interpolates nodes and its weight */
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
	for(n=0, i=0; i<boundary->node->n; i++) {
		for(p=weight_list[i].next; p; p=p->next) {
			weight_info->id[n]     = p->id;
			weight_info->weight[n] = p->weight;
			n++;
		}
	}

	/*
	 * free linked list
	 */
	for(i=0; i<boundary->node->n; i++) {
		free_link_list(weight_list[i].next);
	}
	HECMW_free(weight_list);

	return 0;

error:
	for(i=0; i<boundary->node->n; i++) {
		free_link_list(weight_list[i].next);
	}
	HECMW_free(weight_list);
	return -1;
}



extern struct hecmw_couple_weight_list *
HECMW_couple_s2n_average(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary)
{
	struct hecmw_couple_weight_list *weight_info_list = NULL;
	struct hecmw_couple_weight *weight_info = NULL;
	struct link_list *weight_list = NULL, *p;
	int elem, n_item, size, n, i;

	if(mesh == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_s2n_by_area(): 'mesh' is NULL");
		return NULL;
	}
	if(boundary == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_s2n_by_area(): 'boundary' is NULL");
		return NULL;
	}

	if((weight_info_list = HECMW_couple_alloc_weight_list()) == NULL) return NULL;

	if((weight_info = HECMW_couple_alloc_weight()) == NULL) goto error;
	weight_info_list->info = weight_info;

	if(s2n_average(mesh, boundary, weight_info)) goto error;

	return weight_info_list;

error:
	HECMW_couple_free_weight(weight_info);
	weight_info_list->info = NULL;
	HECMW_couple_free_weight_list(weight_info_list);
	return NULL;
}
