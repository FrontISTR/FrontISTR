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
#include "hecmw_malloc.h"
#include "hecmw_error.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_comm.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_intra_iftable.h"

/*================================================================================================*/

extern void
HECMW_couple_free_intra_iftable(struct hecmw_couple_intra_iftable *p)
{
	if(p == NULL) return;

	HECMW_free(p->neighbor_pe);
	HECMW_free(p->import_index);
	HECMW_free(p->import_item);
	HECMW_free(p->export_index);
	HECMW_free(p->export_item);
	HECMW_free(p);
	p = NULL;
}



extern struct hecmw_couple_intra_iftable *
HECMW_couple_alloc_intra_iftable(void)
{
	struct hecmw_couple_intra_iftable *p = NULL;
	int size;

	size = sizeof(struct hecmw_couple_intra_iftable);
	p = (struct hecmw_couple_intra_iftable *)HECMW_malloc(size);
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	p->n_neighbor_pe = 0;
	p->neighbor_pe   = NULL;
	p->import_index  = NULL;
	p->import_item   = NULL;
	p->export_index  = NULL;
	p->export_item   = NULL;

	return p;
}


/*================================================================================================*/

static int
set_n_neighbor_pe(const struct hecmwST_local_mesh *mesh,
		struct hecmw_couple_intra_iftable *intra_tbl)
{
	intra_tbl->n_neighbor_pe = mesh->n_neighbor_pe;

	return 0;
}



static int
set_neighbor_pe(const struct hecmwST_local_mesh *mesh,
		struct hecmw_couple_intra_iftable *intra_tbl)
{
	int i;

	if(intra_tbl->n_neighbor_pe == 0) return 0;

	intra_tbl->neighbor_pe = (int *)HECMW_malloc(sizeof(int)*intra_tbl->n_neighbor_pe);
	if(intra_tbl->neighbor_pe == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<mesh->n_neighbor_pe; i++) {
		intra_tbl->neighbor_pe[i] = mesh->neighbor_pe[i];
	}

	return 0;
}



static int *
mask_boundary_node(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary)
{
	int *is_boundary_node = NULL;
	int node, i;

	is_boundary_node = (int *)HECMW_malloc(sizeof(int)*mesh->n_node);
	if(is_boundary_node == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	for(i=0; i<mesh->n_node; i++) {
		is_boundary_node[i] = -1;
	}
	for(i=0; i<boundary->node->n; i++) {
		node = boundary->node->item[i];
		is_boundary_node[node-1] = i;
	}

	return is_boundary_node;
}



static int *
mask_import_node(const struct hecmwST_local_mesh *mesh, const int *is_boundary_node)
{
	int *is_import_node = NULL;
	int node, i;


	HECMW_assert(mesh->n_neighbor_pe > 0);

	is_import_node = (int *)HECMW_calloc(mesh->import_index[mesh->n_neighbor_pe]+1, sizeof(int));
	if(is_import_node == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	for(i=0; i<mesh->import_index[mesh->n_neighbor_pe]; i++) {
		is_import_node[i] = -1;
	}
	for(i=0; i<mesh->import_index[mesh->n_neighbor_pe]; i++) {
		node = mesh->import_item[i];
		if(is_boundary_node[node-1] >= 0) is_import_node[i] = is_boundary_node[node-1];
	}

	return is_import_node;
}



static int *
mask_export_node(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_comm *intracomm, int *is_import_node)
{
	int *is_export_node = NULL;
	int rtc;

	is_export_node = (int *)HECMW_calloc(mesh->export_index[mesh->n_neighbor_pe]+1, sizeof(int));
	if(is_export_node == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	rtc = HECMW_couple_intra_send_recv(mesh->n_neighbor_pe, mesh->neighbor_pe,
			mesh->import_index, is_import_node,
			mesh->export_index, is_export_node, HECMW_INT, intracomm->comm);
	if(rtc != 0) {
		HECMW_free(is_export_node);
		return NULL;
	}

	return is_export_node;
}



static int
set_intracomm_import_node(const struct hecmwST_local_mesh *mesh,
		const int *is_boundary_node, const int *is_import_node,
		struct hecmw_couple_intra_iftable *intra_tbl)
{
	int n_import, node, size, n, i, j;

	HECMW_assert(intra_tbl->n_neighbor_pe > 0);

	/* index */
	intra_tbl->import_index = (int *)HECMW_calloc(intra_tbl->n_neighbor_pe+1, sizeof(int));
	if(intra_tbl->import_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(n=0, i=0; i<mesh->n_neighbor_pe; i++) {
		for(j=mesh->import_index[i]; j<mesh->import_index[i+1]; j++) {
			if(is_import_node[j] >= 0) n++;
		}
		intra_tbl->import_index[i+1] = n;
	}

	/* item */
	n_import = intra_tbl->import_index[intra_tbl->n_neighbor_pe];
	if(n_import == 0) return 0;

	intra_tbl->import_item = (int *)HECMW_malloc(sizeof(int)*n_import);
	if(intra_tbl->import_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(n=0, i=0; i<mesh->n_neighbor_pe; i++) {
		for(j=mesh->import_index[i]; j<mesh->import_index[i+1]; j++) {
			node = mesh->import_item[j];
			if(is_import_node[j] >= 0) {
				intra_tbl->import_item[n++] = is_boundary_node[node-1];
			}
		}
		HECMW_assert(n == intra_tbl->import_index[i+1]);
	}

	return 0;
}



static int
set_intracomm_export_node(const struct hecmwST_local_mesh *mesh,
		const int *is_boundary_node, const int *is_export_node,
		struct hecmw_couple_intra_iftable *intra_tbl)
{
	int n_export, size, node, n, i, j;

	HECMW_assert(intra_tbl->n_neighbor_pe > 0);

	/* index */
	intra_tbl->export_index = (int *)HECMW_calloc(intra_tbl->n_neighbor_pe+1, sizeof(int));
	if(intra_tbl->export_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(n=0, i=0; i<mesh->n_neighbor_pe; i++) {
		for(j=mesh->export_index[i]; j<mesh->export_index[i+1]; j++) {
			node = mesh->export_item[j];
			if(is_export_node[j] >= 0) n++;
		}
		intra_tbl->export_index[i+1] = n;
	}

	/* item */
	n_export = intra_tbl->export_index[intra_tbl->n_neighbor_pe];
	if(n_export == 0) return 0;

	intra_tbl->export_item = (int *)HECMW_malloc(sizeof(int)*n_export);
	if(intra_tbl->export_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(n=0, i=0; i<mesh->n_neighbor_pe; i++) {
		for(j=mesh->export_index[i]; j<mesh->export_index[i+1]; j++) {
			node = mesh->export_item[j];
			if(is_export_node[j] >= 0) {
				intra_tbl->export_item[n++] = is_boundary_node[node-1];
			}
		}
		HECMW_assert(n == intra_tbl->export_index[i+1]);
	}

	return 0;
}


/*================================================================================================*/

extern struct hecmw_couple_intra_iftable *
HECMW_couple_make_intra_iftable(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary,
		const struct hecmw_couple_comm *intracomm)
{
	struct hecmw_couple_intra_iftable *intra_tbl;
	int *is_boundary_node = NULL, *is_import_node = NULL, *is_export_node = NULL;

	if(mesh == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_intra_iftable(): 'mesh' is NULL");
		return NULL;
	}
	if(boundary == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_intra_iftable(): 'boundary' is NULL");
		return NULL;
	}
	if(intracomm == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"HECMW_couple_make_intra_iftable(): 'intracomm' is NULL");
		return NULL;
	}

	/* allocation */
	if((intra_tbl = HECMW_couple_alloc_intra_iftable()) == NULL) return NULL;

	/* set neighboring process information */
	if(set_n_neighbor_pe(mesh, intra_tbl)) goto error;
	if(intra_tbl->n_neighbor_pe == 0) return intra_tbl;

	if(set_neighbor_pe(mesh, intra_tbl)) goto error;

	/* masking */
	if((is_boundary_node = mask_boundary_node(mesh, boundary)) == NULL) goto error;
	if((is_import_node = mask_import_node(mesh, is_boundary_node)) == NULL) goto error;
	if((is_export_node = mask_export_node(mesh, intracomm, is_import_node)) == NULL) goto error;

	/* set import/export node information */
	if(set_intracomm_import_node(mesh, is_boundary_node, is_import_node, intra_tbl)) goto error;
	if(set_intracomm_export_node(mesh, is_boundary_node, is_export_node, intra_tbl)) goto error;

	HECMW_free(is_import_node);
	HECMW_free(is_export_node);
	HECMW_free(is_boundary_node);
	return intra_tbl;

error:
	HECMW_free(is_import_node);
	HECMW_free(is_export_node);
	HECMW_free(is_boundary_node);
	HECMW_couple_free_intra_iftable(intra_tbl);
	return NULL;
}
