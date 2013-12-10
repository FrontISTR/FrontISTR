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

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_control.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_bounding_box.h"


#define INFINITE (1.0E+37)

/*================================================================================================*/

extern void
HECMW_couple_free_bounding_box(struct hecmw_couple_bounding_box *bbox)
{
	if(bbox == NULL) return;

	HECMW_free(bbox->just);
	HECMW_free(bbox->enlarged);
	HECMW_free(bbox);
	bbox = NULL;
}



static struct hecmw_couple_bounding_box *
alloc_struct_bbox(void)
{
	struct hecmw_couple_bounding_box *bbox = NULL;
	int size;

	size = sizeof(struct hecmw_couple_bounding_box);
	bbox = (struct hecmw_couple_bounding_box *)HECMW_malloc(size);
	if(bbox == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	bbox->just = NULL;
	bbox->enlarged = NULL;

	bbox->just = (struct hecmw_couple_box *)HECMW_malloc(sizeof(struct hecmw_couple_box));
	if(bbox->just == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	bbox->enlarged = (struct hecmw_couple_box *)HECMW_malloc(sizeof(struct hecmw_couple_box));
	if(bbox->enlarged == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	bbox->coef = 0.0;
	bbox->tolerance = 0.0;

	bbox->just->min_x = 0.0;
	bbox->just->min_y = 0.0;
	bbox->just->min_z = 0.0;
	bbox->just->max_x = 0.0;
	bbox->just->max_y = 0.0;
	bbox->just->max_z = 0.0;

	bbox->enlarged->min_x = 0.0;
	bbox->enlarged->min_y = 0.0;
	bbox->enlarged->min_z = 0.0;
	bbox->enlarged->max_x = 0.0;
	bbox->enlarged->max_y = 0.0;
	bbox->enlarged->max_z = 0.0;

	return bbox;

error:
	HECMW_couple_free_bounding_box(bbox);
	return NULL;
}


/*================================================================================================*/

extern struct hecmw_couple_bounding_box *
HECMW_couple_set_bounding_box(const char *boundary_id,
		const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary)
{
	struct hecmw_couple_bounding_box *bbox = NULL;
	double length_x, length_y, length_z, coord_x, coord_y, coord_z, half_coef;
	double min_x, min_y, min_z, max_x, max_y, max_z;
	int elem, node, i, j;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_bounding_box(): 'boundary_id' is NULL");
		return NULL;
	}
	if(mesh == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_bounding_box(): 'mesh' is NULL");
		return NULL;
	}
	if(boundary == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_set_bounding_box(): 'boundary' is NULL");
		return NULL;
	}

	/* allocation & initialization */
	if((bbox = alloc_struct_bbox()) == NULL) return NULL;

	HECMW_couple_ctrl_get_tolerance(boundary_id, &bbox->tolerance);
	if(bbox->tolerance < 0.0) goto error;
	HECMW_couple_ctrl_get_bbcoef(boundary_id, &bbox->coef);
	if(bbox->coef < 0.0) goto error;

	if(boundary->node->n == 0) return bbox;


	min_x = min_y = min_z = +INFINITE;
	max_x = max_y = max_z = -INFINITE;

	if(boundary->geom_type == HECMW_COUPLE_NODE_GROUP) {			/* node group 		*/
		HECMW_set_error(HECMWCPL_E_NONSUPPORT_GEOMTYPE,
				"In current version, node group is not supported");
		goto error;

		/*
		for(i=0; i<boundary->node->n; i++) {
			node = boundary->node->item[i];
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
		*/

	} else if(boundary->geom_type == HECMW_COUPLE_ELEMENT_GROUP) {	/* element group	*/
		HECMW_set_error(HECMWCPL_E_NONSUPPORT_GEOMTYPE,
				"In current version, element group is not supported");
		goto error;

		/*
		for(i=0; i<boundary->node->n; i++) {
			node = boundary->node->item[i];
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
		*/

	} else if(boundary->geom_type == HECMW_COUPLE_SURFACE_GROUP) {	/* surface group	*/
		for(i=0; i<boundary->surf->n; i++) {
			elem = boundary->surf->item[2*i];
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
		}

	} else {														/* error			*/
		HECMW_set_error(HECMWCPL_E_INVALID_GEOMTYPE, "");
		goto error;
	}

	/* just size bounding box */
	bbox->just->min_x = min_x;
	bbox->just->min_y = min_y;
	bbox->just->min_z = min_z;
	bbox->just->max_x = max_x;
	bbox->just->max_y = max_y;
	bbox->just->max_z = max_z;

	/* enlarged size bounding box */
	half_coef = (bbox->coef-1.0) * 0.5;
	length_x = bbox->just->max_x - bbox->just->min_x;
	length_y = bbox->just->max_y - bbox->just->min_y;
	length_z = bbox->just->max_z - bbox->just->min_z;

	if(length_x > half_coef) {
		bbox->enlarged->min_x = min_x - length_x * half_coef;
		bbox->enlarged->max_x = max_x + length_x * half_coef;
	} else {
		bbox->enlarged->min_x = min_x - half_coef;
		bbox->enlarged->max_x = max_x + half_coef;
	}
	if(length_y > half_coef) {
		bbox->enlarged->min_y = min_y - length_y * half_coef;
		bbox->enlarged->max_y = max_y + length_y * half_coef;
	} else {
		bbox->enlarged->min_y = min_y - half_coef;
		bbox->enlarged->max_y = max_y + half_coef;
	}
	if(length_z > half_coef) {
		bbox->enlarged->min_z = min_z - length_z * half_coef;
		bbox->enlarged->max_z = max_z + length_z * half_coef;
	} else {
		bbox->enlarged->min_z = min_z - half_coef;
		bbox->enlarged->max_z = max_z + half_coef;
	}

	return bbox;

error:
	HECMW_couple_free_bounding_box(bbox);
	return NULL;
}
