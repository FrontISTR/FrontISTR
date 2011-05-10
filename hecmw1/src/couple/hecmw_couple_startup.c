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
#include "hecmw_malloc.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_control.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_info.h"
#include "hecmw_couple_init.h"
#include "hecmw_couple_startup.h"

/*================================================================================================*/

extern void
HECMW_couple_free_couple_value(struct hecmw_couple_value *couple_value)
{
	if(couple_value == NULL) return;

	HECMW_free(couple_value->item);
	HECMW_free(couple_value->value);
	HECMW_free(couple_value);
	couple_value = NULL;
}



extern struct hecmw_couple_value *
HECMW_couple_alloc_couple_value(void)
{
	struct hecmw_couple_value *p = NULL;

	p = (struct hecmw_couple_value *)HECMW_malloc(sizeof(struct hecmw_couple_value));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	p->n         = 0;
	p->n_dof     = 0;
	p->item_type = HECMW_COUPLE_GROUP_UNDEF;
	p->item      = NULL;
	p->value     = NULL;

	return p;
}



extern void
HECMW_couple_print_couple_value(const struct hecmw_couple_value *couple_value, FILE *fp)
{
	int i, j;

	if(couple_value == NULL || fp == NULL) return;

	fprintf(fp, "*** Value of coupling area\n");

	fprintf(fp, "number of item: %d\n", couple_value->n);

	if(couple_value->item_type == HECMW_COUPLE_NODE_GROUP) {
		fprintf(fp, "item type: NODE GROUP\n");
	} else if(couple_value->item_type == HECMW_COUPLE_ELEMENT_GROUP) {
		fprintf(fp, "item type: ELEMENT GROUP\n");
	} else if(couple_value->item_type == HECMW_COUPLE_SURFACE_GROUP) {
		fprintf(fp, "item type: SURFACE GROUP\n");
	} else {
		fprintf(fp, "item type: UNKNOWN\n");
	}

	fprintf(fp, "number of DOF: %d\n", couple_value->n_dof);

	fprintf(fp, "ID & value:\n");
	for(i=0; i<couple_value->n; i++) {
		if(couple_value->item_type == HECMW_COUPLE_SURFACE_GROUP) {
			fprintf(fp, " %d %d", couple_value->item[2*i], couple_value->item[2*i+1]);
		} else {
			fprintf(fp, " %d", couple_value->item[i]);
		}

		for(j=0; j<couple_value->n_dof; j++) {
			fprintf(fp, " %lE", couple_value->value[i*couple_value->n_dof+j]);
		}
		fprintf(fp, "\n");
	}
}


/*------------------------------------------------------------------------------------------------*/

extern struct hecmw_couple_value *
HECMW_couple_startup(const char *boundary_id)
{
	struct hecmw_couple_value *couple_value = NULL;
	struct hecmw_couple_info *couple_info = NULL;
	struct hecmw_couple_boundary *boundary = NULL;
	int i;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple_startup(): 'boundary_id' is NULL");
		return NULL;
	}

	if((couple_value = HECMW_couple_alloc_couple_value()) == NULL) goto error;
	if((couple_info = HECMW_couple_get_info(boundary_id)) == NULL) goto error;

	if(couple_info->comm_src->is_member) {
		boundary = couple_info->boundary_src;

		couple_value->n_dof     = 0;
		couple_value->item_type = boundary->data_type;

		/* node group */
		if(boundary->data_type == HECMW_COUPLE_NODE_GROUP) {
			couple_value->n = boundary->node->n;

			if(couple_value->n > 0) {
				couple_value->item = (int *)HECMW_malloc(sizeof(int)*couple_value->n);
				if(couple_value->item == NULL) {
					HECMW_set_error(errno, "");
					goto error;
				}
				for(i=0; i<couple_value->n; i++) {
					couple_value->item[i] = boundary->node->item[i];
				}
			}

		/* element group */
		} else if(boundary->data_type == HECMW_COUPLE_ELEMENT_GROUP) {
			couple_value->n = boundary->elem->n;

			if(couple_value->n > 0) {
				couple_value->item = (int *)HECMW_malloc(sizeof(int)*couple_value->n);
				if(couple_value->item == NULL) {
					HECMW_set_error(errno, "");
					goto error;
				}
				for(i=0; i<couple_value->n; i++) {
					couple_value->item[i] = boundary->elem->item[i];
				}
			}

		/* surface group */
		} else if(boundary->data_type == HECMW_COUPLE_SURFACE_GROUP) {
			couple_value->n = boundary->surf->n;

			if(couple_value->n > 0) {
				couple_value->item = (int *)HECMW_malloc(sizeof(int)*couple_value->n*2);
				if(couple_value->item == NULL) {
					HECMW_set_error(errno, "");
					goto error;
				}
				for(i=0; i<couple_value->n; i++) {
					couple_value->item[2*i]   = boundary->surf->item[2*i];
					couple_value->item[2*i+1] = boundary->surf->item[2*i+1];
				}
			}

		/* invalid group type */
		} else {
			HECMW_set_error(HECMWCPL_E_INVALID_GRPTYPE, "");
			goto error;
		}
	}

	return couple_value;

error:
	HECMW_couple_free_couple_value(couple_value);
	return NULL;
}


/*------------------------------------------------------------------------------------------------*/

extern void
HECMW_couple_cleanup(struct hecmw_couple_value *couple_value)
{
	HECMW_couple_free_couple_value(couple_value);
}
