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
#include "hecmw_struct.h"
#include "hecmw_error.h"
#include "hecmw_comm.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_control.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_bounding_box.h"
#include "hecmw_couple_background_cell.h"


#define INFINITE (1.0E+37)


#define EPS (1.0E-06)

/*================================================================================================*/

extern void
HECMW_couple_free_background_cell(struct hecmw_couple_background_cell *bgcell)
{
	if(bgcell == NULL)  return;

	HECMW_free(bgcell);
	bgcell = NULL;
}



static struct hecmw_couple_background_cell *
alloc_struct_bgcell(void)
{
	struct hecmw_couple_background_cell *bgcell = NULL;
	int size;

	size = sizeof(struct hecmw_couple_background_cell);
	bgcell = (struct hecmw_couple_background_cell *)HECMW_malloc(size);
	if(bgcell == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	bgcell->n    = 0;
	bgcell->coef = 0.0;
	bgcell->nx   = 0;
	bgcell->ny   = 0;
	bgcell->nz   = 0;
	bgcell->dx   = 0.0;
	bgcell->dy   = 0.0;
	bgcell->dz   = 0.0;

	return bgcell;
}



static int
elem_size_by_elem(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary,
		double *min_dx, double *min_dy, double *min_dz,
		double *max_dx, double *max_dy, double *max_dz)
{
	double min_x, min_y, min_z, max_x, max_y, max_z, dx, dy, dz, coord_x, coord_y, coord_z;
	int elem, node, i, j;

	*min_dx = *min_dy = *min_dz = INFINITE;
	*max_dx = *max_dy = *max_dz = 0.0;

	for(i=0; i<boundary->elem->n; i++) {
		elem = boundary->elem->item[i];

		min_x = min_y = min_z = +INFINITE;
		max_x = max_y = max_z = -INFINITE;

		for(j=mesh->elem_node_index[elem-1]; j<mesh->elem_node_index[elem]; j++) {
			node = mesh->elem_node_item[j];
			coord_x = mesh->node[3*(node-1)  ];
			coord_y = mesh->node[3*(node-1)+1];
			coord_z = mesh->node[3*(node-1)+2];

			if(coord_x < min_x) min_x = coord_x;
			if(coord_y < min_y) min_y = coord_y;
			if(coord_z < min_z) min_z = coord_z;
			if(coord_x > max_x) max_x = coord_x;
			if(coord_y > max_y) max_y = coord_y;
			if(coord_z > max_z) max_z = coord_z;
		}

		dx = max_x - min_x;
		dy = max_y - min_y;
		dz = max_z - min_z;

		if(dx < *min_dx) *min_dx = dx;
		if(dy < *min_dy) *min_dy = dy;
		if(dz < *min_dz) *min_dz = dz;
		if(dx > *max_dx) *max_dx = dx;
		if(dy > *max_dy) *max_dy = dy;
		if(dz > *max_dz) *max_dz = dz;
	}

	return 0;
}



static int
elem_size_by_surf(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary,
		double *min_dx, double *min_dy, double *min_dz,
		double *max_dx, double *max_dy, double *max_dz)
{
	double min_x, min_y, min_z, max_x, max_y, max_z, dx, dy, dz, coord_x, coord_y, coord_z;
	int elem, node, i, j;

	*min_dx = *min_dy = *min_dz = INFINITE;
	*max_dx = *max_dy = *max_dz = 0.0;

	for(i=0; i<boundary->surf->n; i++) {
		elem = boundary->surf->item[2*i];

		min_x = min_y = min_z = +INFINITE;
		max_x = max_y = max_z = -INFINITE;

		for(j=mesh->elem_node_index[elem-1]; j<mesh->elem_node_index[elem]; j++) {
			node = mesh->elem_node_item[j];
			coord_x = mesh->node[3*(node-1)  ];
			coord_y = mesh->node[3*(node-1)+1];
			coord_z = mesh->node[3*(node-1)+2];

			if(coord_x < min_x) min_x = coord_x;
			if(coord_y < min_y) min_y = coord_y;
			if(coord_z < min_z) min_z = coord_z;
			if(coord_x > max_x) max_x = coord_x;
			if(coord_y > max_y) max_y = coord_y;
			if(coord_z > max_z) max_z = coord_z;
		}

		dx = max_x - min_x;
		dy = max_y - min_y;
		dz = max_z - min_z;

		if(dx < *min_dx) *min_dx = dx;
		if(dy < *min_dy) *min_dy = dy;
		if(dz < *min_dz) *min_dz = dz;
		if(dx > *max_dx) *max_dx = dx;
		if(dy > *max_dy) *max_dy = dy;
		if(dz > *max_dz) *max_dz = dz;
	}

	return 0;
}


/*================================================================================================*/

extern struct hecmw_couple_background_cell *
HECMW_couple_set_background_cell(const char *boundary_id,
		const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_bounding_box *bbox,
		const struct hecmw_couple_boundary *boundary)
{
	struct hecmw_couple_background_cell *bgcell;
	double min_dx, min_dy, min_dz, max_dx, max_dy, max_dz, dx, dy, dz;
	double bbox_size_x, bbox_size_y, bbox_size_z;
	int nx, ny, nz;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_background_cell(): 'boundary_id' is NULL");
		return NULL;
	}
	if(mesh == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_background_cell(): 'mesh' is NULL");
		return NULL;
	}
	if(bbox == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_background_cell(): 'bbox' is NULL");
		return NULL;
	}
	if(boundary == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_background_cell(): 'boundary' is NULL");
		return NULL;
	}

	if((bgcell = alloc_struct_bgcell()) == NULL) return NULL;

	if(HECMW_couple_ctrl_get_bgcoef(boundary_id, &bgcell->coef) != HECMW_SUCCESS) goto error;

	if(boundary->geom_type == HECMW_COUPLE_NODE_GROUP) {				/* node group		*/
		HECMW_set_error(HECMWCPL_E_NONSUPPORT_GEOMTYPE,
				"In current version, node group is not supported.");
		goto error;

	} else if(boundary->geom_type == HECMW_COUPLE_ELEMENT_GROUP) {		/* element group	*/
		HECMW_set_error(HECMWCPL_E_NONSUPPORT_GEOMTYPE,
				"In current version, element group is not supported.");
		goto error;
/*		elem_size_by_elem(mesh, boundary, &min_dx, &min_dy, &min_dz, &max_dx, &max_dy, &max_dz); */

	} else if(boundary->geom_type == HECMW_COUPLE_SURFACE_GROUP) {		/* surface group	*/
		elem_size_by_surf(mesh, boundary, &min_dx, &min_dy, &min_dz, &max_dx, &max_dy, &max_dz);

	} else {															/* error			*/
		HECMW_set_error(HECMWCPL_E_INVALID_GEOMTYPE, "");
		goto error;
	}

	bbox_size_x = bbox->enlarged->max_x - bbox->enlarged->min_x;
	bbox_size_y = bbox->enlarged->max_y - bbox->enlarged->min_y;
	bbox_size_z = bbox->enlarged->max_z - bbox->enlarged->min_z;

	if(max_dx*bgcell->coef > bbox->tolerance+EPS) {
		dx = max_dx * bgcell->coef;
	} else {
		dx = bbox->tolerance + EPS;
	}
	if(max_dy*bgcell->coef > bbox->tolerance+EPS) {
		dy = max_dy * bgcell->coef;
	} else {
		dy = bbox->tolerance + EPS;
	}
	if(max_dz*bgcell->coef > bbox->tolerance+EPS) {
		dz = max_dz * bgcell->coef;
	} else {
		dz = bbox->tolerance + EPS;
	}

	nx = bbox_size_x / dx;
	ny = bbox_size_y / dy;
	nz = bbox_size_z / dz;

	bgcell->nx = (nx > 0) ? nx : 1;
	bgcell->ny = (ny > 0) ? ny : 1;
	bgcell->nz = (nz > 0) ? nz : 1;

	bgcell->dx = bbox_size_x / bgcell->nx;
	bgcell->dy = bbox_size_y / bgcell->ny;
	bgcell->dz = bbox_size_z / bgcell->nz;

	bgcell->n = bgcell->nx * bgcell->ny * bgcell->nz;

	return bgcell;

error:
	HECMW_couple_free_background_cell(bgcell);
	return NULL;
}
