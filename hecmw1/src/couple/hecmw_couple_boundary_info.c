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
#include "hecmw_common_define.h"
#include "hecmw_struct.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_table.h"
#include "hecmw_couple_control.h"
#include "hecmw_couple_info.h"
#include "hecmw_couple_boundary_info.h"



struct link_list {
	int item;					
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



extern void
HECMW_couple_free_boundary_info(struct hecmw_couple_boundary *boundary)
{
	if(boundary == NULL)  return;

	if(boundary->node) {
		HECMW_free(boundary->node->item);
		HECMW_free(boundary->node);
	}
	if(boundary->elem) {
		HECMW_free(boundary->elem->item);
		HECMW_free(boundary->elem);
	}
	if(boundary->surf) {
		HECMW_free(boundary->surf->item);
		HECMW_free(boundary->surf);
	}
	HECMW_free(boundary->elem_node_index);
	HECMW_free(boundary->elem_node_item);

	HECMW_free(boundary);
	boundary = NULL;
}



extern struct hecmw_couple_boundary *
HECMW_couple_alloc_boundary_info(void)
{
	struct hecmw_couple_boundary *boundary = NULL;
	int size;

	size = sizeof(struct hecmw_couple_boundary);
	boundary = (struct hecmw_couple_boundary *)HECMW_malloc(size);
	if(boundary == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	boundary->node = NULL;
	boundary->elem = NULL;
	boundary->surf = NULL;

	boundary->geom_type = HECMW_COUPLE_GROUP_UNDEF;
	boundary->data_type = HECMW_COUPLE_GROUP_UNDEF;

	size = sizeof(struct hecmw_couple_boundary_item);
	boundary->node = (struct hecmw_couple_boundary_item *)HECMW_malloc(size);
	if(boundary->node == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	boundary->node->n = 0;
	boundary->node->item = NULL;

	boundary->elem = (struct hecmw_couple_boundary_item *)HECMW_malloc(size);
	if(boundary->elem == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	boundary->elem->n = 0;
	boundary->elem->item = NULL;

	boundary->surf = (struct hecmw_couple_boundary_item *)HECMW_malloc(size);
	if(boundary->surf == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	boundary->surf->n = 0;
	boundary->surf->item = NULL;

	boundary->elem_node_item  = NULL;
	boundary->elem_node_index = NULL;

	return boundary;

error:
	HECMW_couple_free_boundary_info(boundary);
	return NULL;
}


/*================================================================================================*/

static int
check_group_name(int n_grp_mesh, char **grp_name_mesh, const char *grp_name_ctrl)
{
	int i;

	for(i=0; i<n_grp_mesh; i++) {
		if((strcmp(grp_name_ctrl, grp_name_mesh[i])) == 0 ) return i;
	}

	return -1;
}



static int
check_node_group_name(const struct hecmwST_local_mesh *mesh, const struct hecmw_couple_group *group)
{
	struct hecmwST_node_grp *node_grp = mesh->node_group;
	int i;

	for(i=0; i<group->n_grp; i++) {
		if(check_group_name(node_grp->n_grp, node_grp->grp_name, group->grp_name[i]) < 0) {
			HECMW_set_error(HECMWCPL_E_UNDEF_GRPNAME, "node group: %s", group->grp_name[i]);
			return -1;
		}
	}

	return 0;
}



static int
check_elem_group_name(const struct hecmwST_local_mesh *mesh, const struct hecmw_couple_group *group)
{
	struct hecmwST_elem_grp *elem_grp = mesh->elem_group;
	int i;

	for(i=0; i<group->n_grp; i++) {
		if(check_group_name(elem_grp->n_grp, elem_grp->grp_name, group->grp_name[i]) < 0) {
			HECMW_set_error(HECMWCPL_E_UNDEF_GRPNAME, "element group: %s", group->grp_name[i]);
			return -1;
		}
	}

	return 0;
}



static int
check_surf_group_name(const struct hecmwST_local_mesh *mesh, const struct hecmw_couple_group *group)
{
	struct hecmwST_surf_grp *surf_grp = mesh->surf_group;
	int i;

	for(i=0; i<group->n_grp; i++) {
		if(check_group_name(surf_grp->n_grp, surf_grp->grp_name, group->grp_name[i]) < 0) {
			HECMW_set_error(HECMWCPL_E_UNDEF_GRPNAME, "surface group: %s", group->grp_name[i]);
			return HECMW_ERROR;
		}
	}

	return 0;
}


/*------------------------------------------------------------------------------------------------*/

static int
set_boundary_node_by_node(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_group *group,
		struct hecmw_couple_boundary *boundary)
{
	int *mask = NULL;
	int node, index, n, i, j;

	HECMW_assert(group->geom_type == HECMW_COUPLE_NODE_GROUP);
	if(group->n_grp == 0) return 0;

	/* mask boundary nodes */
	mask = (int *)HECMW_malloc(sizeof(int)*mesh->n_node);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<mesh->n_node; i++) {
		mask[i] = HECMW_COUPLE_FALSE;
	}

	for(i=0; i<group->n_grp; i++) {
		index = check_group_name(mesh->node_group->n_grp, mesh->node_group->grp_name,
				group->grp_name[i]);
		HECMW_assert(index >= 0);

		for(j=mesh->node_group->grp_index[index]; j<mesh->node_group->grp_index[index+1]; j++) {
			node = mesh->node_group->grp_item[j];
			mask[node-1] = HECMW_COUPLE_TRUE;
		}
	}

	/* number of boundary nodes */
	for(n=0, i=0; i<mesh->n_node; i++) {
		if(mask[i] == HECMW_COUPLE_TRUE) n++;
	}
	boundary->node->n = n;
	if(boundary->node->n == 0) {
		HECMW_free(mask);
		return 0;
	}

	/* ids of boundary node */
	boundary->node->item = (int *)HECMW_calloc(boundary->node->n, sizeof(int));
	if(boundary->node->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<mesh->n_node; i++) {
		if(mask[i] == HECMW_COUPLE_TRUE) boundary->node->item[n++] = j+1;
	}
	HECMW_assert(n == boundary->node->n);

	HECMW_free(mask);
	return 0;

error:
	HECMW_free(mask);
	return -1;
}



static int
set_boundary_elem_by_elem(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_group *group, struct hecmw_couple_boundary *boundary)
{
	int *mask = NULL;
	int elem, index, n, i, j;

	HECMW_assert(group->geom_type == HECMW_COUPLE_ELEMENT_GROUP);
	if(group->n_grp == 0) return 0;

	/* mask boundary nodes */
	mask = (int *)HECMW_malloc(sizeof(int)*mesh->n_elem);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<mesh->n_elem; i++) {
		mask[i] = HECMW_COUPLE_FALSE;
	}

	for(i=0; i<group->n_grp; i++) {
		index = check_group_name(mesh->elem_group->n_grp, mesh->elem_group->grp_name,
				group->grp_name[i]);
		HECMW_assert(index >= 0);

		for(j=mesh->elem_group->grp_index[index]; j<mesh->elem_group->grp_index[index+1]; j++) {
			elem = mesh->elem_group->grp_item[j];
			mask[elem-1] = HECMW_COUPLE_TRUE;
		}
	}

	/* number of boundary nodes */
	for(n=0, i=0; i<mesh->n_elem; i++) {
		if(mask[i] == HECMW_COUPLE_TRUE) n++;
	}
	boundary->elem->n = n;
	if(boundary->elem->n == 0) {
		HECMW_free(mask);
		return 0;
	}

	/* ids of boundary node */
	boundary->elem->item = (int *)HECMW_calloc(boundary->elem->n, sizeof(int));
	if(boundary->elem->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<mesh->n_elem; i++) {
		if(mask[i] == HECMW_COUPLE_TRUE) boundary->elem->item[n++] = j+1;
	}
	HECMW_assert(n == boundary->elem->n);

	HECMW_free(mask);
	return 0;

error:
	HECMW_free(mask);
	return -1;
}



static int
set_boundary_node_by_elem(const struct hecmwST_local_mesh *mesh,
		struct hecmw_couple_boundary *boundary)
{
	int *mask = NULL;
	int elem, node, size, n, i, j;

	/* mask boundary nodes */
	mask = (int *)HECMW_malloc(sizeof(int)*mesh->n_node);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<mesh->n_node; i++) {
		mask[i] = HECMW_COUPLE_FALSE;
	}

	for(n=0, i=0; i<boundary->elem->n; i++) {
		elem = boundary->elem->item[i];
		for(j=mesh->elem_node_index[elem-1]; j<mesh->elem_node_index[elem]; j++) {
			node = mesh->elem_node_item[j];
			if(mask[node-1] < 0) {
				mask[node-1] = n++;
			}
		}
	}

	/* number of boundary nodes */
	for(n=0, i=0; i<mesh->n_node; i++) {
		if(mask[i] >= 0) n++;
	}
	boundary->node->n = n;
	if(boundary->node->n == 0) {
		HECMW_free(mask);
		return 0;
	}

	/* ids of boundary node */
	boundary->node->item = (int *)HECMW_malloc(sizeof(int)*boundary->node->n);
	if(boundary->node->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<mesh->n_node; i++) {
		if(mask[i] >= 0) {
			boundary->node->item[mask[i]] = i+1;
		}
	}

	/* connectivity of component nodes */
	boundary->elem_node_index = (int *)HECMW_calloc(boundary->elem->n+1, sizeof(int));
	if(boundary->elem_node_index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<boundary->elem->n; i++) {
		elem = boundary->elem->item[i];
		boundary->elem_node_index[i+1] =
			boundary->elem_node_index[i] + HECMW_get_max_node(mesh->elem_type[elem-1]);
	}

	size = sizeof(int) * boundary->elem_node_index[boundary->elem->n];
	boundary->elem_node_item = (int *)HECMW_malloc(size);
	if(boundary->elem_node_item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<boundary->elem->n; i++) {
		elem = boundary->elem->item[i];
		for(j=mesh->elem_node_index[elem-1]; j<mesh->elem_node_index[elem]; j++) {
			node = mesh->elem_node_item[j];
			boundary->elem_node_item[n] = mask[node-1];
		}
	}

	HECMW_free(mask);
	return 0;

error:
	HECMW_free(mask);
	return -1;
}



static int
set_boundary_surf_by_surf(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_group *group, struct hecmw_couple_boundary *boundary)
{
	struct link_list *mask, *p;
	int size, index, elem, surf, is_exist, n, i, j;

	HECMW_assert(group->geom_type == HECMW_COUPLE_SURFACE_GROUP);
	if(group->n_grp == 0)  return 0;

	/* mask boundary surfaces */
	mask = (struct link_list *)HECMW_malloc(sizeof(struct link_list)*mesh->n_elem);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<mesh->n_elem; i++) {
		mask[i].item = 0;
		mask[i].next = NULL;
	}

	for(n=0, i=0; i<group->n_grp; i++) {
		index = check_group_name(mesh->surf_group->n_grp, mesh->surf_group->grp_name,
				group->grp_name[i]);
		HECMW_assert(index >= 0);

		for(j=mesh->surf_group->grp_index[index]; j<mesh->surf_group->grp_index[index+1]; j++) {
			elem = mesh->surf_group->grp_item[2*j];
			surf = mesh->surf_group->grp_item[2*j+1];

			is_exist = HECMW_COUPLE_FALSE;
			p        = &mask[elem-1];
			while(p->next) {
				if(p->next->item == surf) {
					is_exist = HECMW_COUPLE_TRUE;
					break;
				}
			}

			if(is_exist == HECMW_COUPLE_FALSE) {
				p->next = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
				if(p->next == NULL) {
					HECMW_set_error(errno, "");
					goto error;
				}
				p->next->item = surf;
				p->next->next = NULL;
				n++;
			}
		}
	}

	/* number of boundary surfaces */
	boundary->surf->n = n;
	if(boundary->surf->n == 0) {
		HECMW_free(mask);
		return 0;
	}

	/* ids of boundary surface */
	boundary->surf->item = (int *)HECMW_calloc(boundary->surf->n*2, sizeof(int));
	if(boundary->surf->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<mesh->n_elem; i++) {
		for(p=mask[i].next; p; p=p->next) {
			boundary->surf->item[2*n  ] = i+1;
			boundary->surf->item[2*n+1] = p->item;
			n++;
		}
	}

	for(i=0; i<mesh->n_elem; i++) {
		free_link_list(mask[i].next);
	}
	HECMW_free(mask);

	return 0;

error:
	for(i=0; i<mesh->n_elem; i++) {
		free_link_list(mask[i].next);
	}
	HECMW_free(mask);

	return -1;
}



static int
set_boundary_node_by_surf(const struct hecmwST_local_mesh *mesh,
		struct hecmw_couple_boundary *boundary)
{
	int *mask = NULL;
	int elem, surf, node, node_index, offset, size, n, i, j;

	/* mask boundary nodes */
	mask = (int *)HECMW_malloc(sizeof(int)*mesh->n_node);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<mesh->n_node; i++) {
		mask[i] = -1;
	}
	for(n=0, i=0; i<boundary->surf->n; i++) {
		elem       = boundary->surf->item[2*i];
		surf       = boundary->surf->item[2*i+1];
		node_index = mesh->elem_node_index[elem-1];

		if(mesh->elem_type[elem-1] == HECMW_ETYPE_TET1) {			/* 1st-order Tetra	*/
			for(j=0; j<3; j++) {
				offset = hecmw_surf_node_table_tet1[surf-1][j];
				node   = mesh->elem_node_item[node_index+offset-1];
				if(mask[node-1] < 0) {
					mask[node-1] = n++;
				}
			}
		} else if(mesh->elem_type[elem-1] == HECMW_ETYPE_HEX1) {	/* 1st-order Hexa	*/
			for(j=0; j<4; j++) {
				offset = hecmw_surf_node_table_hex1[surf-1][j];
				node   = mesh->elem_node_item[node_index+offset-1];
				if(mask[node-1] < 0) {
					mask[node-1] = n++;
				}
			}
		} else {													/* error			*/
			HECMW_set_error(HECMWCPL_E_NONSUPPORT_ETYPE, "%d", mesh->elem_type[elem-1]);
			goto error;
		}
	}

	/* number of boundary nodes */
	for(n=0, i=0; i<mesh->n_node; i++) {
		if(mask[i] >= 0) n++;
	}
	boundary->node->n = n;
	if(boundary->node->n == 0) {
		HECMW_free(mask);
		return 0;
	}

	/* ids of boundary node */
	boundary->node->item = (int *)HECMW_malloc(sizeof(int)*boundary->node->n);
	if(boundary->node->item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<mesh->n_node; i++) {
		if(mask[i] >= 0) {
			boundary->node->item[mask[i]] = i+1;
		}
	}

	/* connectivity of component nodes */
	boundary->elem_node_index = (int *)HECMW_calloc(boundary->surf->n+1, sizeof(int));
	if(boundary->elem_node_index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<boundary->surf->n; i++) {
		elem = boundary->surf->item[2*i];
		if(mesh->elem_type[elem-1] == HECMW_ETYPE_TET1) {
			boundary->elem_node_index[i+1] = boundary->elem_node_index[i] + 3;
		} else if(mesh->elem_type[elem-1] == HECMW_ETYPE_HEX1) {
			boundary->elem_node_index[i+1] = boundary->elem_node_index[i] + 4;
		}
	}

	size = sizeof(int) * boundary->elem_node_index[boundary->surf->n];
	boundary->elem_node_item = (int *)HECMW_malloc(size);
	if(boundary->elem_node_item == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(n=0, i=0; i<boundary->surf->n; i++) {
		elem = boundary->surf->item[2*i];
		surf = boundary->surf->item[2*i+1];
		node_index = mesh->elem_node_index[elem-1];
		if(mesh->elem_type[elem-1] == HECMW_ETYPE_TET1) {			/* 1st-order Tetra	*/
			for(j=0; j<3; j++) {
				offset = hecmw_surf_node_table_tet1[surf-1][j];
				node   = mesh->elem_node_item[node_index+offset-1];
				boundary->elem_node_item[n++] = mask[node-1];
			}
		} else if(mesh->elem_type[elem-1] == HECMW_ETYPE_HEX1) {	/* 1st-order Hexa	*/
			for(j=0; j<4; j++) {
				offset = hecmw_surf_node_table_hex1[surf-1][j];
				node   = mesh->elem_node_item[node_index+offset-1];
				boundary->elem_node_item[n++] = mask[node-1];
			}
		} else {													/* error			*/
			HECMW_set_error(HECMWCPL_E_NONSUPPORT_ETYPE, "%d", mesh->elem_type[elem-1]);
			goto error;
		}
	}

	HECMW_free(mask);
	return 0;

error:
	HECMW_free(mask);
	return -1;
}


/*================================================================================================*/

extern struct hecmw_couple_boundary *
HECMW_couple_set_boundary_info(const char *boundary_id, int unit_specifier,
		const struct hecmwST_local_mesh *mesh)
{
	struct hecmw_couple_boundary *boundary = NULL;
	struct hecmw_couple_group *group = NULL;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_boundary_info(): 'boundary_id' is NULL");
		return NULL;
	}
	if(mesh == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_boundary_info(): 'mesh' is NULL");
		return NULL;
	}

	if((boundary = HECMW_couple_alloc_boundary_info()) == NULL) return NULL;

	if((group = HECMW_couple_ctrl_get_group(boundary_id, unit_specifier)) == NULL) goto error;

	if(group->geom_type == HECMW_COUPLE_NODE_GROUP) {				/* Node Group		 */
		HECMW_set_error(HECMWCPL_E_NONSUPPORT_GEOMTYPE,
				"In current version, node group is not supported");
		goto error;

		/*
		boundary->geom_type = group->geom_type;
		boundary->data_type = group->data_type;
		if(check_node_group_name(mesh, group)) goto error;
		if(set_boundary_node_by_node(mesh, group, boundary)) goto error;
		*/

	} else if(group->geom_type == HECMW_COUPLE_ELEMENT_GROUP) {		/* Element Group	*/
		HECMW_set_error(HECMWCPL_E_NONSUPPORT_GEOMTYPE,
				"In current version, element group is not supported");
		goto error;

		/*
		boundary->geom_type = group->geom_type;
		boundary->data_type = group->data_type;
		if(check_elem_group_name(mesh, group)) goto error;
		if(set_boundary_elem_by_elem(mesh, group, boundary)) goto error;
		if(set_boundary_node_by_elem(mesh, boundary)) goto error; 
		*/

	} else if(group->geom_type == HECMW_COUPLE_SURFACE_GROUP) {		/* Surface Group	*/
		boundary->geom_type = group->geom_type;
		boundary->data_type = group->data_type;
		if(check_surf_group_name(mesh, group)) goto error;
		if(set_boundary_surf_by_surf(mesh, group, boundary)) goto error;
		if(set_boundary_node_by_surf(mesh, boundary)) goto error;

	} else {														/* Error			*/
		HECMW_set_error(HECMWCPL_E_INVALID_GEOMTYPE, "");
		goto error;
	}

	HECMW_couple_ctrl_free_group(group);
	return boundary;

error:
	HECMW_couple_free_boundary_info(boundary);
	HECMW_couple_ctrl_free_group(group);
	return NULL;
}
