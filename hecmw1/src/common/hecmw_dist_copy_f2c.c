/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/05/02                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"

static struct hecmwST_local_mesh *mesh;


/*-----------------------------------------------------------------------------
 * SetFunc
 */

static int
set_hecmw_flag_adapt(void *src)
{
	mesh->hecmw_flag_adapt = *((int *)src);
	return 0;
}


static int
set_hecmw_flag_initcon(void *src)
{
	mesh->hecmw_flag_initcon = *((int *)src);
	return 0;
}


static int
set_hecmw_flag_parttype(void *src)
{
	mesh->hecmw_flag_parttype = *((int *)src);
	return 0;
}


static int
set_hecmw_flag_partdepth(void *src)
{
	mesh->hecmw_flag_partdepth = *((int *)src);
	return 0;
}


static int
set_hecmw_flag_version(void *src)
{
	mesh->hecmw_flag_version = *((int *)src);
	return 0;
}


static int
set_gridfile(void *src)
{
	void *dst = mesh->gridfile;
	HECMW_strcpy_f2c_r(src, HECMW_FILENAME_LEN, dst, HECMW_FILENAME_LEN+1);
	return 0;
}


static int
set_hecmw_n_file(void *src)
{
	mesh->hecmw_n_file = *((int *)src);
	return 0;
}


static int
set_files(void *src)
{
	int i;

	if(mesh->hecmw_n_file <= 0) return 0;

	mesh->files = HECMW_calloc(mesh->hecmw_n_file, sizeof(*mesh->files));
	if(mesh->files == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < mesh->hecmw_n_file; i++) {
		char *src_point = (char *)src + HECMW_FILENAME_LEN * i;
		mesh->files[i] = HECMW_strcpy_f2c(src_point, HECMW_FILENAME_LEN);
		if(mesh->files[i] == NULL) goto error;
	}
	return 0;
error:
	if(mesh->files) {
		for(i=0; i < mesh->hecmw_n_file; i++) {
			HECMW_free(mesh->files[i]);
		}
	}
	HECMW_free(mesh->files);
	mesh->files = NULL;
	return -1;
}


static int
set_header(void *src)
{
	void *dst = mesh->header;
	HECMW_strcpy_f2c_r(src, HECMW_HEADER_LEN, dst, HECMW_HEADER_LEN+1);
	return 0;
}


static int
set_zero_temp(void *src)
{
	mesh->zero_temp = *((double *)src);
	return 0;
}


static int
set_n_node(void *src)
{
	mesh->n_node = *((int *)src);
	return 0;
}


static int
set_n_node_gross(void *src)
{
	mesh->n_node_gross = *((int *)src);
	return 0;
}


static int
set_nn_internal(void *src)
{
	mesh->nn_internal = *((int *)src);
	return 0;
}


static int
set_node_internal_list(void *src)
{
	int size;

	if(mesh->nn_internal <= 0) return 0;
	size = sizeof(*mesh->node_internal_list)*mesh->nn_internal;
	mesh->node_internal_list = HECMW_malloc(size);
	if(mesh->node_internal_list == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_internal_list, src, size);
	return 0;
}


static int
set_node_ID(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->node_ID)*2*mesh->n_node_gross;
	mesh->node_ID = HECMW_malloc(size);
	if(mesh->node_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_ID, src, size);
	return 0;
}


static int
set_global_node_ID(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->global_node_ID)*mesh->n_node_gross;
	mesh->global_node_ID = HECMW_malloc(size);
	if(mesh->global_node_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->global_node_ID, src, size);
	return 0;
}


static int
set_node(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->node)*3*mesh->n_node_gross;
	mesh->node = HECMW_malloc(size);
	if(mesh->node == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node, src, size);
	return 0;
}


static int
set_n_dof(void *src)
{
	mesh->n_dof = *((int *)src);
	return 0;
}


static int
set_n_dof_grp(void *src)
{
	mesh->n_dof_grp = *((int *)src);
	return 0;
}


static int
set_n_dof_tot(void *src)
{
	mesh->n_dof_tot = *((int *)src);
	return 0;
}


static int
set_node_dof_index(void *src)
{
	int size;

	if(mesh->n_dof_grp <= 0) return 0;
	size = sizeof(*mesh->node_dof_index)*(mesh->n_dof_grp+1);
	mesh->node_dof_index = HECMW_malloc(size);
	if(mesh->node_dof_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_dof_index, src, size);
	return 0;
}


static int
set_node_dof_item(void *src)
{
	int size;

	if(mesh->n_dof_grp <= 0) return 0;
	size = sizeof(*mesh->node_dof_item)*mesh->n_dof_grp;
	mesh->node_dof_item = HECMW_malloc(size);
	if(mesh->node_dof_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_dof_item, src, size);
	return 0;
}


static int
set_node_val_index(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->node_val_index)*(mesh->n_node_gross+1);
	mesh->node_val_index = HECMW_malloc(size);
	if(mesh->node_val_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_val_index, src, size);
	return 0;
}


static int
set_node_val_item(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	if(mesh->node_val_index == NULL) return 0;
	if(mesh->node_val_index[mesh->n_node_gross] <= 0) return 0;

	size = sizeof(*mesh->node_val_item)*mesh->node_val_index[mesh->n_node_gross];
	mesh->node_val_item = HECMW_malloc(size);
	if(mesh->node_val_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_val_item, src, size);
	return 0;
}


static int
set_node_init_val_index(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->node_init_val_index)*(mesh->n_node_gross+1);
	mesh->node_init_val_index = HECMW_malloc(size);
	if(mesh->node_init_val_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_init_val_index, src, size);
	return 0;
}


static int
set_node_init_val_item(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	if(mesh->node_init_val_index == NULL) return 0;
	if(mesh->node_init_val_index[mesh->n_node_gross] <= 0) return 0;

	size = sizeof(*mesh->node_init_val_item)*mesh->node_init_val_index[mesh->n_node_gross];
	mesh->node_init_val_item = HECMW_malloc(size);
	if(mesh->node_init_val_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_init_val_item, src, size);
	return 0;
}


static int
set_n_elem(void *src)
{
	mesh->n_elem = *((int *)src);
	return 0;
}


static int
set_n_elem_gross(void *src)
{
	mesh->n_elem_gross = *((int *)src);
	return 0;
}


static int
set_ne_internal(void *src)
{
	mesh->ne_internal = *((int *)src);
	return 0;
}


static int
set_elem_internal_list(void *src)
{
	int size;

	if(mesh->ne_internal <= 0) return 0;
	size = sizeof(*mesh->elem_internal_list)*mesh->ne_internal;
	mesh->elem_internal_list = HECMW_malloc(size);
	if(mesh->elem_internal_list == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_internal_list, src, size);
	return 0;
}


static int
set_elem_ID(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_ID)*2*mesh->n_elem_gross;
	mesh->elem_ID = HECMW_malloc(size);
	if(mesh->elem_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_ID, src, size);
	return 0;
}


static int
set_global_elem_ID(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->global_elem_ID)*mesh->n_elem_gross;
	mesh->global_elem_ID = HECMW_malloc(size);
	if(mesh->global_elem_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->global_elem_ID, src, size);
	return 0;
}


static int
set_elem_type(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_type)*mesh->n_elem_gross;
	mesh->elem_type = HECMW_malloc(size);
	if(mesh->elem_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_type, src, size);
	return 0;
}


static int
set_n_elem_type(void *src)
{
	mesh->n_elem_type = *((int *)src);
	return 0;
}


static int
set_elem_type_index(void *src)
{
	int size;

	if(mesh->n_elem_type <= 0) return 0;
	size = sizeof(*mesh->elem_type_index)*(mesh->n_elem_type+1);
	mesh->elem_type_index = HECMW_malloc(size);
	if(mesh->elem_type_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_type_index, src, size);
	return 0;
}


static int
set_elem_type_item(void *src)
{
	int size;

	if(mesh->n_elem_type <= 0) return 0;
	size = sizeof(*mesh->elem_type_item)*mesh->n_elem_type;
	mesh->elem_type_item = HECMW_malloc(size);
	if(mesh->elem_type_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_type_item, src, size);
	return 0;
}


static int
set_elem_node_index(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_node_index)*(mesh->n_elem_gross+1);
	mesh->elem_node_index = HECMW_malloc(size);
	if(mesh->elem_node_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_node_index, src, size);
	return 0;
}


static int
set_elem_node_item(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	if(mesh->elem_node_index == NULL) return 0;
	if(mesh->elem_node_index[mesh->n_elem_gross] <= 0) return 0;

	size = sizeof(*mesh->elem_node_item)*mesh->elem_node_index[mesh->n_elem_gross];
	mesh->elem_node_item = HECMW_malloc(size);
	if(mesh->elem_node_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_node_item, src, size);
	return 0;
}


static int
set_section_ID(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->section_ID)*mesh->n_elem_gross;
	mesh->section_ID = HECMW_malloc(size);
	if(mesh->section_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section_ID, src, size);
	return 0;
}


static int
set_n_elem_mat_ID(void *src)
{
	mesh->n_elem_mat_ID = *((int *)src);
	return 0;
}


static int
set_elem_mat_ID_index(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_mat_ID_index)*(mesh->n_elem_gross+1);
	mesh->elem_mat_ID_index = HECMW_malloc(size);
	if(mesh->elem_mat_ID_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_mat_ID_index, src, size);
	return 0;
}


static int
set_elem_mat_ID_item(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	if(mesh->elem_mat_ID_index == NULL) return 0;
	if(mesh->elem_mat_ID_index[mesh->n_elem_gross] <= 0) return 0;

	size = sizeof(*mesh->elem_mat_ID_item)*mesh->elem_mat_ID_index[mesh->n_elem_gross];
	mesh->elem_mat_ID_item = HECMW_malloc(size);
	if(mesh->elem_mat_ID_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_mat_ID_item, src, size);
	return 0;
}


static int
set_elem_mat_int_index(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_mat_int_index)*(mesh->n_elem_gross+1);
	mesh->elem_mat_int_index = HECMW_malloc(size);
	if(mesh->elem_mat_int_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_mat_int_index, src, size);
	return 0;
}


static int
set_elem_mat_int_val(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	if(mesh->elem_mat_int_index == NULL) return 0;
	if(mesh->elem_mat_int_index[mesh->n_elem_gross] <= 0) return 0;

	size = sizeof(*mesh->elem_mat_int_val)*mesh->elem_mat_int_index[mesh->n_elem_gross];
	mesh->elem_mat_int_val = HECMW_malloc(size);
	if(mesh->elem_mat_int_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_mat_int_val, src, size);
	return 0;
}


static int
set_elem_val_index(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_val_index)*(mesh->n_elem_gross+1);
	mesh->elem_val_index = HECMW_malloc(size);
	if(mesh->elem_val_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_val_index, src, size);
	return 0;
}


static int
set_elem_val_item(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	if(mesh->elem_val_index == NULL) return 0;
	if(mesh->elem_val_index[mesh->n_elem_gross] <= 0) return 0;

	size = sizeof(*mesh->elem_val_item)*mesh->elem_val_index[mesh->n_elem_gross];
	mesh->elem_val_item = HECMW_malloc(size);
	if(mesh->elem_val_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_val_item, src, size);
	return 0;
}


static int
set_zero(void *src)
{
	mesh->zero = *((int *)src);
	return 0;
}


static int
set_HECMW_COMM(void *src)
{
	mesh->HECMW_COMM = HECMW_Comm_f2c(*((HECMW_Fint *)src));
	return 0;
}


static int
set_PETOT(void *src)
{
	mesh->PETOT = *((int *)src);
	return 0;
}


static int
set_PEsmpTOT(void *src)
{
	mesh->PEsmpTOT = *((int *)src);
	return 0;
}


static int
set_my_rank(void *src)
{
	mesh->my_rank = *((int *)src);
	return 0;
}


static int
set_errnof(void *src)
{
	mesh->errnof = *((int *)src);
	return 0;
}


static int
set_n_subdomain(void *src)
{
	mesh->n_subdomain = *((int *)src);
	return 0;
}


static int
set_n_neighbor_pe(void *src)
{
	mesh->n_neighbor_pe = *((int *)src);
	return 0;
}


static int
set_neighbor_pe(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	size = sizeof(*mesh->neighbor_pe)*mesh->n_neighbor_pe;
	mesh->neighbor_pe = HECMW_malloc(size);
	if(mesh->neighbor_pe == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->neighbor_pe, src, size);
	return 0;
}

static int
set_import_index(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	size = sizeof(*mesh->import_index)*(mesh->n_neighbor_pe+1);
	mesh->import_index = HECMW_malloc(size);
	if(mesh->import_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->import_index, src, size);
	return 0;
}


static int
set_import_item(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	if(mesh->import_index == NULL) return 0;
	if(mesh->import_index[mesh->n_neighbor_pe] <= 0) return 0;

	size = sizeof(*mesh->import_item)*mesh->import_index[mesh->n_neighbor_pe];
	mesh->import_item = HECMW_malloc(size);
	if(mesh->import_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->import_item, src, size);
	return 0;
}


static int
set_export_index(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	size = sizeof(*mesh->export_index)*(mesh->n_neighbor_pe+1);
	mesh->export_index = HECMW_malloc(size);
	if(mesh->export_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->export_index, src, size);
	return 0;
}


static int
set_export_item(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	if(mesh->export_index == NULL) return 0;
	if(mesh->export_index[mesh->n_neighbor_pe] <= 0) return 0;

	size = sizeof(*mesh->export_item)*mesh->export_index[mesh->n_neighbor_pe];
	mesh->export_item = HECMW_malloc(size);
	if(mesh->export_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->export_item, src, size);
	return 0;
}


static int
set_shared_index(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	size = sizeof(*mesh->shared_index)*(mesh->n_neighbor_pe+1);
	mesh->shared_index = HECMW_malloc(size);
	if(mesh->shared_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->shared_index, src, size);
	return 0;
}


static int
set_shared_item(void *src)
{
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;
	if(mesh->shared_index == NULL) return 0;
	if(mesh->shared_index[mesh->n_neighbor_pe] <= 0) return 0;

	size = sizeof(*mesh->shared_item)*mesh->shared_index[mesh->n_neighbor_pe];
	mesh->shared_item = HECMW_malloc(size);
	memcpy(mesh->shared_item, src, size);
	return 0;
}


static int
set_coarse_grid_level(void *src)
{
	mesh->coarse_grid_level = *((int *)src);
	return 0;
}


static int
set_n_adapt(void *src)
{
	mesh->n_adapt = *((int *)src);
	return 0;
}


static int
set_when_i_was_refined_node(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->when_i_was_refined_node)*mesh->n_node_gross;
	mesh->when_i_was_refined_node = HECMW_malloc(size);
	if(mesh->when_i_was_refined_node == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->when_i_was_refined_node, src, size);
	return 0;
}


static int
set_when_i_was_refined_elem(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->when_i_was_refined_elem)*mesh->n_elem_gross;
	mesh->when_i_was_refined_elem = HECMW_malloc(size);
	if(mesh->when_i_was_refined_elem == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->when_i_was_refined_elem, src, size);
	return 0;
}


static int
set_adapt_parent_type(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->adapt_parent_type)*mesh->n_elem_gross;
	mesh->adapt_parent_type = HECMW_malloc(size);
	if(mesh->adapt_parent_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->adapt_parent_type, src, size);
	return 0;
}


static int
set_adapt_type(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->adapt_type)*mesh->n_elem_gross;
	mesh->adapt_type = HECMW_malloc(size);
	if(mesh->adapt_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->adapt_type, src, size);
	return 0;
}


static int
set_adapt_level(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->adapt_level)*mesh->n_elem_gross;
	mesh->adapt_level = HECMW_malloc(size);
	if(mesh->adapt_level == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->adapt_level, src, size);
	return 0;
}


static int
set_adapt_parent(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->adapt_parent)*2*mesh->n_elem_gross;
	mesh->adapt_parent = HECMW_malloc(size);
	if(mesh->adapt_parent == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->adapt_parent, src, size);
	return 0;
}


static int
set_adapt_children_index(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->adapt_children_index)*(mesh->n_elem_gross+1);
	mesh->adapt_children_index = HECMW_malloc(size);
	if(mesh->adapt_children_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->adapt_children_index, src, size);
	return 0;
}


static int
set_adapt_children_item(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	if(mesh->adapt_children_index == NULL) return 0;
	if(mesh->adapt_children_index[mesh->n_elem_gross] <= 0) return 0;

	size = sizeof(*mesh->adapt_children_item)*2*mesh->adapt_children_index[mesh->n_elem_gross];
	mesh->adapt_children_item = HECMW_malloc(size);
	if(mesh->adapt_children_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->adapt_children_item, src, size);
	return 0;
}


static int
set_n_refine(void *src)
{
	mesh->n_refine = *((int *)src);
	return 0;
}


static int
set_node_old2new(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->node_old2new)*mesh->n_node_gross;
	mesh->node_old2new = HECMW_malloc(size);
	if(mesh->node_old2new == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_old2new, src, size);
	return 0;
}


static int
set_node_new2old(void *src)
{
	int size;

	if(mesh->n_node_gross <= 0) return 0;
	size = sizeof(*mesh->node_new2old)*mesh->n_node_gross;
	mesh->node_new2old = HECMW_malloc(size);
	if(mesh->node_new2old == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_new2old, src, size);
	return 0;
}


static int
set_elem_old2new(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_old2new)*mesh->n_elem_gross;
	mesh->elem_old2new = HECMW_malloc(size);
	if(mesh->elem_old2new == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_old2new, src, size);
	return 0;
}


static int
set_elem_new2old(void *src)
{
	int size;

	if(mesh->n_elem_gross <= 0) return 0;
	size = sizeof(*mesh->elem_new2old)*mesh->n_elem_gross;
	mesh->elem_new2old = HECMW_malloc(size);
	if(mesh->elem_new2old == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_new2old, src, size);
	return 0;
}


static int
set_n_sect(void *src)
{
	mesh->section->n_sect = *((int *)src);
	return 0;
}


static int
set_sect_type(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	size = sizeof(*mesh->section->sect_type)*mesh->section->n_sect;
	mesh->section->sect_type = HECMW_malloc(size);
	if(mesh->section->sect_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_type, src, size);
	return 0;
}


static int
set_sect_opt(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	size = sizeof(*mesh->section->sect_opt)*mesh->section->n_sect;
	mesh->section->sect_opt = HECMW_malloc(size);
	if(mesh->section->sect_opt == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_opt, src, size);
	return 0;
}


static int
set_sect_mat_ID_index(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	size = sizeof(*mesh->section->sect_mat_ID_index)*(mesh->section->n_sect+1);
	mesh->section->sect_mat_ID_index = HECMW_malloc(size);
	if(mesh->section->sect_mat_ID_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_mat_ID_index, src, size);
	return 0;
}


static int
set_sect_mat_ID_item(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	if(mesh->section->sect_mat_ID_index == NULL) return 0;
	if(mesh->section->sect_mat_ID_index[mesh->section->n_sect] <= 0) return 0;

	size = sizeof(*mesh->section->sect_mat_ID_item)*mesh->section->sect_mat_ID_index[mesh->section->n_sect];
	mesh->section->sect_mat_ID_item = HECMW_malloc(size);
	if(mesh->section->sect_mat_ID_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_mat_ID_item, src, size);
	return 0;
}


static int
set_sect_I_index(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	size = sizeof(*mesh->section->sect_I_index)*(mesh->section->n_sect+1);
	mesh->section->sect_I_index = HECMW_malloc(size);
	if(mesh->section->sect_I_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_I_index, src, size);
	return 0;
}


static int
set_sect_I_item(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	if(mesh->section->sect_I_index == NULL) return 0;
	if(mesh->section->sect_I_index[mesh->section->n_sect] <= 0) return 0;

	size = sizeof(*mesh->section->sect_I_item)*mesh->section->sect_I_index[mesh->section->n_sect];
	mesh->section->sect_I_item = HECMW_malloc(size);
	if(mesh->section->sect_I_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_I_item, src, size);
	return 0;
}


static int
set_sect_R_index(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	size = sizeof(*mesh->section->sect_R_index)*(mesh->section->n_sect+1);
	mesh->section->sect_R_index = HECMW_malloc(size);
	if(mesh->section->sect_R_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_R_index, src, size);
	return 0;
}


static int
set_sect_R_item(void *src)
{
	int size;

	if(mesh->section->n_sect <= 0) return 0;
	if(mesh->section->sect_R_index == NULL) return 0;
	if(mesh->section->sect_R_index[mesh->section->n_sect] <= 0) return 0;

	size = sizeof(*mesh->section->sect_R_item)*mesh->section->sect_R_index[mesh->section->n_sect];
	mesh->section->sect_R_item = HECMW_malloc(size);
	if(mesh->section->sect_R_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->section->sect_R_item, src, size);
	return 0;
}


static int
set_n_mat(void *src)
{
	mesh->material->n_mat = *((int *)src);
	return 0;
}


static int
set_n_mat_item(void *src)
{
	mesh->material->n_mat_item = *((int *)src);
	return 0;
}


static int
set_n_mat_subitem(void *src)
{
	mesh->material->n_mat_subitem = *((int *)src);
	return 0;
}


static int
set_n_mat_table(void *src)
{
	mesh->material->n_mat_table = *((int *)src);
	return 0;
}


static int
set_mat_name(void *src)
{
	int i;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat <= 0) return 0;

	mat->mat_name = HECMW_calloc(mat->n_mat, sizeof(*mat->mat_name));
	if(mat->mat_name == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < mat->n_mat; i++) {
		char *src_point = (char *)src + HECMW_NAME_LEN * i;
		mat->mat_name[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
		if(mat->mat_name[i] == NULL) goto error;
	}
	return 0;
error:
	if(mat->mat_name) {
		for(i=0; i < mat->n_mat; i++) {
			HECMW_free(mat->mat_name[i]);
		}
	}
	HECMW_free(mat->mat_name);
	mat->mat_name = NULL;
	return -1;
}


static int
set_mat_item_index(void *src)
{
	int size;

	if(mesh->material->n_mat <= 0) return 0;
	size = sizeof(*mesh->material->mat_item_index)*(mesh->material->n_mat+1);
	mesh->material->mat_item_index = HECMW_malloc(size);
	if(mesh->material->mat_item_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->material->mat_item_index, src, size);
	return 0;
}


static int
set_mat_subitem_index(void *src)
{
	int size;

	if(mesh->material->n_mat_item <= 0) return 0;
	size = sizeof(*mesh->material->mat_subitem_index)*(mesh->material->n_mat_item+1);
	mesh->material->mat_subitem_index = HECMW_malloc(size);
	if(mesh->material->mat_subitem_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->material->mat_subitem_index, src, size);
	return 0;
}


static int
set_mat_table_index(void *src)
{
	int size;

	if(mesh->material->n_mat_subitem <= 0) return 0;
	size = sizeof(*mesh->material->mat_table_index)*(mesh->material->n_mat_subitem+1);
	mesh->material->mat_table_index = HECMW_malloc(size);
	if(mesh->material->mat_table_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->material->mat_table_index, src, size);
	return 0;
}


static int
set_mat_val(void *src)
{
	int size;

	if(mesh->material->n_mat_table <= 0) return 0;
	size = sizeof(*mesh->material->mat_val)*mesh->material->n_mat_table;
	mesh->material->mat_val = HECMW_malloc(size);
	if(mesh->material->mat_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->material->mat_val, src, size);
	return 0;
}


static int
set_mat_temp(void *src)
{
	int size;

	if(mesh->material->n_mat_table <= 0) return 0;
	size = sizeof(*mesh->material->mat_temp)*mesh->material->n_mat_table;
	mesh->material->mat_temp = HECMW_malloc(size);
	if(mesh->material->mat_temp == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->material->mat_temp, src, size);
	return 0;
}


static int
set_n_mpc(void *src)
{
	mesh->mpc->n_mpc = *((int *)src);
	return 0;
}


static int
set_mpc_index(void *src)
{
	int size;

	if(mesh->mpc->n_mpc <= 0) return 0;
	size = sizeof(*mesh->mpc->mpc_index)*(mesh->mpc->n_mpc+1);
	mesh->mpc->mpc_index = HECMW_malloc(size);
	if(mesh->mpc->mpc_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->mpc->mpc_index, src, size);
	return 0;
}


static int
set_mpc_item(void *src)
{
	int size;

	if(mesh->mpc->n_mpc <= 0) return 0;
	if(mesh->mpc->mpc_index == NULL) return 0;
	if(mesh->mpc->mpc_index[mesh->mpc->n_mpc] <= 0) return 0;

	size = sizeof(*mesh->mpc->mpc_item)*mesh->mpc->mpc_index[mesh->mpc->n_mpc];
	mesh->mpc->mpc_item = HECMW_malloc(size);
	if(mesh->mpc->mpc_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->mpc->mpc_item, src, size);
	return 0;
}


static int
set_mpc_dof(void *src)
{
	int size;

	if(mesh->mpc->n_mpc <= 0) return 0;
	if(mesh->mpc->mpc_index == NULL) return 0;
	if(mesh->mpc->mpc_index[mesh->mpc->n_mpc] <= 0) return 0;

	size = sizeof(*mesh->mpc->mpc_dof)*mesh->mpc->mpc_index[mesh->mpc->n_mpc];
	mesh->mpc->mpc_dof = HECMW_malloc(size);
	if(mesh->mpc->mpc_dof == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->mpc->mpc_dof, src, size);
	return 0;
}


static int
set_mpc_val(void *src)
{
	int size;

	if(mesh->mpc->n_mpc <= 0) return 0;
	if(mesh->mpc->mpc_index == NULL) return 0;
	if(mesh->mpc->mpc_index[mesh->mpc->n_mpc] <= 0) return 0;

	size = sizeof(*mesh->mpc->mpc_val)*mesh->mpc->mpc_index[mesh->mpc->n_mpc];
	mesh->mpc->mpc_val = HECMW_malloc(size);
	if(mesh->mpc->mpc_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->mpc->mpc_val, src, size);
	return 0;
}


static int
set_mpc_const(void *src)
{
	int size;

	if(mesh->mpc->n_mpc <= 0) return 0;

	size = sizeof(*mesh->mpc->mpc_const)*mesh->mpc->n_mpc;
	mesh->mpc->mpc_const = HECMW_malloc(size);
	if(mesh->mpc->mpc_const == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->mpc->mpc_const, src, size);
	return 0;
}


static int
set_n_amp(void *src)
{
	mesh->amp->n_amp = *((int *)src);
	return 0;
}


static int
set_amp_name(void *src)
{
	int i;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	amp->amp_name = HECMW_calloc(amp->n_amp, sizeof(*amp->amp_name));
	if(amp->amp_name == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < amp->n_amp; i++) {
		char *src_point = (char *)src + HECMW_NAME_LEN * i;
		amp->amp_name[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
		if(amp->amp_name[i] == NULL) goto error;
	}
	return 0;
error:
	if(amp->amp_name) {
		for(i=0; i < amp->n_amp; i++) {
			HECMW_free(amp->amp_name[i]);
		}
	}
	HECMW_free(amp->amp_name);
	amp->amp_name = NULL;
	return -1;
}


static int
set_amp_type_definition(void *src)
{
	int size;

	if(mesh->amp->n_amp <= 0) return 0;
	size = sizeof(*mesh->amp->amp_type_definition)*mesh->amp->n_amp;
	mesh->amp->amp_type_definition = HECMW_malloc(size);
	if(mesh->amp->amp_type_definition == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->amp->amp_type_definition, src, size);
	return 0;
}


static int
set_amp_type_time(void *src)
{
	int size;

	if(mesh->amp->n_amp <= 0) return 0;
	size = sizeof(*mesh->amp->amp_type_time)*mesh->amp->n_amp;
	mesh->amp->amp_type_time = HECMW_malloc(size);
	if(mesh->amp->amp_type_time == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->amp->amp_type_time, src, size);
	return 0;
}


static int
set_amp_type_value(void *src)
{
	int size;

	if(mesh->amp->n_amp <= 0) return 0;
	size = sizeof(*mesh->amp->amp_type_value)*mesh->amp->n_amp;
	mesh->amp->amp_type_value = HECMW_malloc(size);
	if(mesh->amp->amp_type_value == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->amp->amp_type_value, src, size);
	return 0;
}


static int
set_amp_index(void *src)
{
	int size;

	if(mesh->amp->n_amp <= 0) return 0;
	size = sizeof(*mesh->amp->amp_index)*(mesh->amp->n_amp+1);
	mesh->amp->amp_index = HECMW_malloc(size);
	if(mesh->amp->amp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->amp->amp_index, src, size);
	return 0;
}


static int
set_amp_val(void *src)
{
	int size;

	if(mesh->amp->n_amp <= 0) return 0;
	if(mesh->amp->amp_index == NULL) return 0;
	if(mesh->amp->amp_index[mesh->amp->n_amp] <= 0) return 0;

	size = sizeof(*mesh->amp->amp_val)*mesh->amp->amp_index[mesh->amp->n_amp];
	mesh->amp->amp_val = HECMW_malloc(size);
	if(mesh->amp->amp_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->amp->amp_val, src, size);
	return 0;
}


static int
set_amp_table(void *src)
{
	int size;

	if(mesh->amp->n_amp <= 0) return 0;
	if(mesh->amp->amp_index == NULL) return 0;
	if(mesh->amp->amp_index[mesh->amp->n_amp] <= 0) return 0;

	size = sizeof(*mesh->amp->amp_table)*mesh->amp->amp_index[mesh->amp->n_amp];
	mesh->amp->amp_table = HECMW_malloc(size);
	if(mesh->amp->amp_table == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->amp->amp_table, src, size);
	return 0;
}


static int
set_ngrp_n_grp(void *src)
{
	mesh->node_group->n_grp = *((int *)src);
	return 0;
}


static int
set_ngrp_n_bc(void *src)
{
	mesh->node_group->n_bc = *((int *)src);
	return 0;
}


static int
set_ngrp_grp_name(void *src)
{
	int i;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_grp <= 0) return 0;

	grp->grp_name = HECMW_calloc(grp->n_grp, sizeof(*grp->grp_name));
	if(grp->grp_name == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < grp->n_grp; i++) {
		char *src_point = (char *)src + HECMW_NAME_LEN * i;
		grp->grp_name[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
		if(grp->grp_name[i] == NULL) goto error;
	}
	return 0;
error:
	if(grp->grp_name) {
		for(i=0; i < grp->n_grp; i++) {
			HECMW_free(grp->grp_name[i]);
		}
	}
	HECMW_free(grp->grp_name);
	grp->grp_name = NULL;
	return -1;
}


static int
set_ngrp_grp_index(void *src)
{
	int size;

	if(mesh->node_group->n_grp <= 0) return 0;
	size = sizeof(*mesh->node_group->grp_index)*(mesh->node_group->n_grp+1);
	mesh->node_group->grp_index = HECMW_malloc(size);
	if(mesh->node_group->grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->grp_index, src, size);
	return 0;
}


static int
set_ngrp_grp_item(void *src)
{
	int size;

	if(mesh->node_group->n_grp <= 0) return 0;
	if(mesh->node_group->grp_index == NULL) return 0;
	if(mesh->node_group->grp_index[mesh->node_group->n_grp] <= 0) return 0;

	size = sizeof(*mesh->node_group->grp_item)*mesh->node_group->grp_index[mesh->node_group->n_grp];
	mesh->node_group->grp_item = HECMW_malloc(size);
	if(mesh->node_group->grp_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->grp_item, src, size);
	return 0;
}


static int
set_ngrp_bc_grp_ID(void *src)
{
	int size;

	if(mesh->node_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->node_group->bc_grp_ID)*mesh->node_group->n_bc;
	mesh->node_group->bc_grp_ID = HECMW_malloc(size);
	if(mesh->node_group->bc_grp_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->bc_grp_ID, src, size);
	return 0;
}


static int
set_ngrp_bc_grp_type(void *src)
{
	int size;

	if(mesh->node_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->node_group->bc_grp_type)*mesh->node_group->n_bc;
	mesh->node_group->bc_grp_type = HECMW_malloc(size);
	if(mesh->node_group->bc_grp_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->bc_grp_type, src, size);
	return 0;
}


static int
set_ngrp_bc_grp_index(void *src)
{
	int size;

	if(mesh->node_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->node_group->bc_grp_index)*mesh->node_group->n_bc;
	mesh->node_group->bc_grp_index = HECMW_malloc(size);
	if(mesh->node_group->bc_grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->bc_grp_index, src, size);
	return 0;
}


static int
set_ngrp_bc_grp_dof(void *src)
{
	int size;

	if(mesh->node_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->node_group->bc_grp_dof)*mesh->node_group->n_bc;
	mesh->node_group->bc_grp_dof = HECMW_malloc(size);
	if(mesh->node_group->bc_grp_dof == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->bc_grp_dof, src, size);
	return 0;
}


static int
set_ngrp_bc_grp_val(void *src)
{
	int size;

	if(mesh->node_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->node_group->bc_grp_val)*mesh->node_group->n_bc;
	mesh->node_group->bc_grp_val = HECMW_malloc(size);
	if(mesh->node_group->bc_grp_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->node_group->bc_grp_val, src, size);
	return 0;
}


static int
set_egrp_n_grp(void *src)
{
	mesh->elem_group->n_grp = *((int *)src);
	return 0;
}


static int
set_egrp_n_bc(void *src)
{
	mesh->elem_group->n_bc = *((int *)src);
	return 0;
}


static int
set_egrp_grp_name(void *src)
{
	int i;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_grp <= 0) return 0;

	grp->grp_name = HECMW_calloc(grp->n_grp, sizeof(*grp->grp_name));
	if(grp->grp_name == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < grp->n_grp; i++) {
		char *src_point = (char *)src + HECMW_NAME_LEN * i;
		grp->grp_name[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
		if(grp->grp_name[i] == NULL) goto error;
	}
	return 0;
error:
	if(grp->grp_name) {
		for(i=0; i < grp->n_grp; i++) {
			HECMW_free(grp->grp_name[i]);
		}
	}
	HECMW_free(grp->grp_name);
	grp->grp_name = NULL;
	return -1;
}


static int
set_egrp_grp_index(void *src)
{
	int size;

	if(mesh->elem_group->n_grp <= 0) return 0;
	size = sizeof(*mesh->elem_group->grp_index)*(mesh->elem_group->n_grp+1);
	mesh->elem_group->grp_index = HECMW_malloc(size);
	if(mesh->elem_group->grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_group->grp_index, src, size);
	return 0;
}


static int
set_egrp_grp_item(void *src)
{
	int size;

	if(mesh->elem_group->n_grp <= 0) return 0;
	if(mesh->elem_group->grp_index == NULL) return 0;
	if(mesh->elem_group->grp_index[mesh->elem_group->n_grp] <= 0) return 0;

	size = sizeof(*mesh->elem_group->grp_item)*mesh->elem_group->grp_index[mesh->elem_group->n_grp];
	mesh->elem_group->grp_item = HECMW_malloc(size);
	if(mesh->elem_group->grp_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_group->grp_item, src, size);
	return 0;
}


static int
set_egrp_bc_grp_ID(void *src)
{
	int size;

	if(mesh->elem_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->elem_group->bc_grp_ID)*mesh->elem_group->n_bc;
	mesh->elem_group->bc_grp_ID = HECMW_malloc(size);
	if(mesh->elem_group->bc_grp_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_group->bc_grp_ID, src, size);
	return 0;
}


static int
set_egrp_bc_grp_type(void *src)
{
	int size;

	if(mesh->elem_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->elem_group->bc_grp_type)*mesh->elem_group->n_bc;
	mesh->elem_group->bc_grp_type = HECMW_malloc(size);
	if(mesh->elem_group->bc_grp_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_group->bc_grp_type, src, size);
	return 0;
}


static int
set_egrp_bc_grp_index(void *src)
{
	int size;

	if(mesh->elem_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->elem_group->bc_grp_index)*mesh->elem_group->n_bc;
	mesh->elem_group->bc_grp_index = HECMW_malloc(size);
	if(mesh->elem_group->bc_grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_group->bc_grp_index, src, size);
	return 0;
}


static int
set_egrp_bc_grp_val(void *src)
{
	int size;

	if(mesh->elem_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->elem_group->bc_grp_val)*mesh->elem_group->n_bc;
	mesh->elem_group->bc_grp_val = HECMW_malloc(size);
	if(mesh->elem_group->bc_grp_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->elem_group->bc_grp_val, src, size);
	return 0;
}


static int
set_sgrp_n_grp(void *src)
{
	mesh->surf_group->n_grp = *((int *)src);
	return 0;
}


static int
set_sgrp_n_bc(void *src)
{
	mesh->surf_group->n_bc = *((int *)src);
	return 0;
}


static int
set_sgrp_grp_name(void *src)
{
	int i;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_grp <= 0) return 0;

	grp->grp_name = HECMW_calloc(grp->n_grp, sizeof(*grp->grp_name));
	if(grp->grp_name == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < grp->n_grp; i++) {
		char *src_point = (char *)src + HECMW_NAME_LEN * i;
		grp->grp_name[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
		if(grp->grp_name[i] == NULL) goto error;
	}
	return 0;
error:
	if(grp->grp_name) {
		for(i=0; i < grp->n_grp; i++) {
			HECMW_free(grp->grp_name[i]);
		}
	}
	HECMW_free(grp->grp_name);
	grp->grp_name = NULL;
	return -1;
}


static int
set_sgrp_grp_index(void *src)
{
	int size;

	if(mesh->surf_group->n_grp <= 0) return 0;
	size = sizeof(*mesh->surf_group->grp_index)*(mesh->surf_group->n_grp+1);
	mesh->surf_group->grp_index = HECMW_malloc(size);
	if(mesh->surf_group->grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->surf_group->grp_index, src, size);
	return 0;
}


static int
set_sgrp_grp_item(void *src)
{
	int size;

	if(mesh->surf_group->n_grp <= 0) return 0;
	if(mesh->surf_group->grp_index == NULL) return 0;
	if(mesh->surf_group->grp_index[mesh->surf_group->n_grp] <= 0) return 0;

	size = sizeof(*mesh->surf_group->grp_item)*2*mesh->surf_group->grp_index[mesh->surf_group->n_grp];
	mesh->surf_group->grp_item = HECMW_malloc(size);
	if(mesh->surf_group->grp_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->surf_group->grp_item, src, size);
	return 0;
}


static int
set_sgrp_bc_grp_ID(void *src)
{
	int size;

	if(mesh->surf_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->surf_group->bc_grp_ID)*mesh->surf_group->n_bc;
	mesh->surf_group->bc_grp_ID = HECMW_malloc(size);
	if(mesh->surf_group->bc_grp_ID == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->surf_group->bc_grp_ID, src, size);
	return 0;
}


static int
set_sgrp_bc_grp_type(void *src)
{
	int size;

	if(mesh->surf_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->surf_group->bc_grp_type)*mesh->surf_group->n_bc;
	mesh->surf_group->bc_grp_type = HECMW_malloc(size);
	if(mesh->surf_group->bc_grp_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->surf_group->bc_grp_type, src, size);
	return 0;
}


static int
set_sgrp_bc_grp_index(void *src)
{
	int size;

	if(mesh->surf_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->surf_group->bc_grp_index)*2*mesh->surf_group->n_bc;
	mesh->surf_group->bc_grp_index = HECMW_malloc(size);
	if(mesh->surf_group->bc_grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->surf_group->bc_grp_index, src, size);
	return 0;
}


static int
set_sgrp_bc_grp_val(void *src)
{
	int size;

	if(mesh->surf_group->n_bc <= 0) return 0;
	size = sizeof(*mesh->surf_group->bc_grp_val)*mesh->surf_group->n_bc;
	mesh->surf_group->bc_grp_val = HECMW_malloc(size);
	if(mesh->surf_group->bc_grp_val == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->surf_group->bc_grp_val, src, size);
	return 0;
}


static int
set_contact_pair_n_pair(void *src)
{
	mesh->contact_pair->n_pair = *((int *)src);
	return 0;
}


static int
set_contact_pair_name(void *src)
{
	int i;
	struct hecmwST_contact_pair *cpair = mesh->contact_pair;

	if(cpair->n_pair <= 0) return 0;

	cpair->name = HECMW_calloc(cpair->n_pair, sizeof(*cpair->name));
	if(cpair->name == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i < cpair->n_pair; i++) {
		char *src_point = (char *)src + HECMW_NAME_LEN * i;
		cpair->name[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
		if(cpair->name[i] == NULL) goto error;
	}
	return 0;
error:
	if(cpair->name) {
		for(i=0; i < cpair->n_pair; i++) {
			HECMW_free(cpair->name[i]);
		}
	}
	HECMW_free(cpair->name);
	cpair->name = NULL;
	return -1;
}


static int
set_contact_pair_type(void *src)
{
	int size;

	if(mesh->contact_pair->n_pair <= 0) return 0;
	size = sizeof(*mesh->contact_pair->type)*(mesh->contact_pair->n_pair);
	mesh->contact_pair->type = HECMW_malloc(size);
	if(mesh->contact_pair->type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->contact_pair->type, src, size);
	return 0;
}


static int
set_contact_pair_slave_grp_id(void *src)
{
	int size;

	if(mesh->contact_pair->n_pair <= 0) return 0;
	size = sizeof(*mesh->contact_pair->slave_grp_id)*(mesh->contact_pair->n_pair);
	mesh->contact_pair->slave_grp_id = HECMW_malloc(size);
	if(mesh->contact_pair->slave_grp_id == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->contact_pair->slave_grp_id, src, size);
	return 0;
}


static int
set_contact_pair_master_grp_id(void *src)
{
	int size;

	if(mesh->contact_pair->n_pair <= 0) return 0;
	size = sizeof(*mesh->contact_pair->master_grp_id)*(mesh->contact_pair->n_pair);
	mesh->contact_pair->master_grp_id = HECMW_malloc(size);
	if(mesh->contact_pair->master_grp_id == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(mesh->contact_pair->master_grp_id, src, size);
	return 0;
}


/*-----------------------------------------------------------------------------
 * SetFunc table
 */

typedef int (*SetFunc)(void *);

static struct func_table {
	char *struct_name;
	char *var_name;
	SetFunc set_func;
} functions[] = {
/*  { Struct name, Variable name, memcpy function } */
	{"hecmwST_local_mesh","hecmw_flag_adapt",    set_hecmw_flag_adapt    },
	{"hecmwST_local_mesh","hecmw_flag_initcon",  set_hecmw_flag_initcon  },
	{"hecmwST_local_mesh","hecmw_flag_parttype", set_hecmw_flag_parttype },
	{"hecmwST_local_mesh","hecmw_flag_partdepth",set_hecmw_flag_partdepth},
	{"hecmwST_local_mesh","hecmw_flag_version",  set_hecmw_flag_version  },

	{"hecmwST_local_mesh","gridfile",            set_gridfile            },
	{"hecmwST_local_mesh","hecmw_n_file",        set_hecmw_n_file        },
	{"hecmwST_local_mesh","files",               set_files               },
	{"hecmwST_local_mesh","header",              set_header              },
	{"hecmwST_local_mesh","zero_temp",           set_zero_temp           },

	{"hecmwST_local_mesh","n_node",              set_n_node              },
	{"hecmwST_local_mesh","n_node_gross",        set_n_node_gross        },
	{"hecmwST_local_mesh","nn_internal",         set_nn_internal         },
	{"hecmwST_local_mesh","node_internal_list",  set_node_internal_list  },
	{"hecmwST_local_mesh","node_ID",             set_node_ID             },
	{"hecmwST_local_mesh","global_node_ID",      set_global_node_ID      },
	{"hecmwST_local_mesh","node",                set_node                },
	{"hecmwST_local_mesh","n_dof",               set_n_dof               },
	{"hecmwST_local_mesh","n_dof_grp",           set_n_dof_grp           },
	{"hecmwST_local_mesh","n_dof_tot",           set_n_dof_tot           },
	{"hecmwST_local_mesh","node_dof_index",      set_node_dof_index      },
	{"hecmwST_local_mesh","node_dof_item",       set_node_dof_item       },
	{"hecmwST_local_mesh","node_val_index",      set_node_val_index      },
	{"hecmwST_local_mesh","node_val_item",       set_node_val_item       },
	{"hecmwST_local_mesh","node_init_val_index", set_node_init_val_index },
	{"hecmwST_local_mesh","node_init_val_item",  set_node_init_val_item  },
	{"hecmwST_local_mesh","n_elem",              set_n_elem              },
	{"hecmwST_local_mesh","n_elem_gross",        set_n_elem_gross        },
	{"hecmwST_local_mesh","ne_internal",         set_ne_internal         },
	{"hecmwST_local_mesh","elem_internal_list",  set_elem_internal_list  },
	{"hecmwST_local_mesh","elem_ID",             set_elem_ID             },
	{"hecmwST_local_mesh","global_elem_ID",      set_global_elem_ID      },
	{"hecmwST_local_mesh","elem_type",           set_elem_type           },
	{"hecmwST_local_mesh","n_elem_type",         set_n_elem_type         },
	{"hecmwST_local_mesh","elem_type_index",     set_elem_type_index     },
	{"hecmwST_local_mesh","elem_type_item",      set_elem_type_item      },
	{"hecmwST_local_mesh","elem_node_index",     set_elem_node_index     },
	{"hecmwST_local_mesh","elem_node_item",      set_elem_node_item      },
	{"hecmwST_local_mesh","section_ID",          set_section_ID          },
	{"hecmwST_local_mesh","n_elem_mat_ID",       set_n_elem_mat_ID       },
	{"hecmwST_local_mesh","elem_mat_ID_index",   set_elem_mat_ID_index   },
	{"hecmwST_local_mesh","elem_mat_ID_item",    set_elem_mat_ID_item    },
	{"hecmwST_local_mesh","elem_mat_int_index",  set_elem_mat_int_index  },
	{"hecmwST_local_mesh","elem_mat_int_val",    set_elem_mat_int_val    },
	{"hecmwST_local_mesh","elem_val_index",      set_elem_val_index      },
	{"hecmwST_local_mesh","elem_val_item",       set_elem_val_item       },

	{"hecmwST_local_mesh","zero",                set_zero                },
	{"hecmwST_local_mesh","HECMW_COMM",          set_HECMW_COMM          },
	{"hecmwST_local_mesh","PETOT",               set_PETOT               },
	{"hecmwST_local_mesh","PEsmpTOT",            set_PEsmpTOT            },
	{"hecmwST_local_mesh","my_rank",             set_my_rank             },
	{"hecmwST_local_mesh","errnof",              set_errnof              },
	{"hecmwST_local_mesh","n_subdomain",         set_n_subdomain         },
	{"hecmwST_local_mesh","n_neighbor_pe",       set_n_neighbor_pe       },
	{"hecmwST_local_mesh","neighbor_pe",         set_neighbor_pe         },
	{"hecmwST_local_mesh","import_index",        set_import_index        },
	{"hecmwST_local_mesh","import_item",         set_import_item         },
	{"hecmwST_local_mesh","export_index",        set_export_index        },
	{"hecmwST_local_mesh","export_item",         set_export_item         },
	{"hecmwST_local_mesh","shared_index",        set_shared_index        },
	{"hecmwST_local_mesh","shared_item",         set_shared_item         },

	{"hecmwST_local_mesh","coarse_grid_level",   set_coarse_grid_level   },
	{"hecmwST_local_mesh","n_adapt",             set_n_adapt             },
	{"hecmwST_local_mesh","when_i_was_refined_node",set_when_i_was_refined_node},
	{"hecmwST_local_mesh","when_i_was_refined_elem",set_when_i_was_refined_elem},
	{"hecmwST_local_mesh","adapt_parent_type",   set_adapt_parent_type   },
	{"hecmwST_local_mesh","adapt_type",          set_adapt_type          },
	{"hecmwST_local_mesh","adapt_level",         set_adapt_level         },
	{"hecmwST_local_mesh","adapt_parent",        set_adapt_parent        },
	{"hecmwST_local_mesh","adapt_children_index",set_adapt_children_index},
	{"hecmwST_local_mesh","adapt_children_item", set_adapt_children_item },

	{"hecmwST_local_mesh","n_refine",            set_n_refine            },
	{"hecmwST_local_mesh","node_old2new",        set_node_old2new        },
	{"hecmwST_local_mesh","node_new2old",        set_node_new2old        },
	{"hecmwST_local_mesh","elem_old2new",        set_elem_old2new        },
	{"hecmwST_local_mesh","elem_new2old",        set_elem_new2old        },

	{"hecmwST_section","n_sect",                 set_n_sect              },
	{"hecmwST_section","sect_type",              set_sect_type           },
	{"hecmwST_section","sect_opt",               set_sect_opt            },
	{"hecmwST_section","sect_mat_ID_index",      set_sect_mat_ID_index   },
	{"hecmwST_section","sect_mat_ID_item",       set_sect_mat_ID_item    },
	{"hecmwST_section","sect_I_index",           set_sect_I_index        },
	{"hecmwST_section","sect_I_item",            set_sect_I_item         },
	{"hecmwST_section","sect_R_index",           set_sect_R_index        },
	{"hecmwST_section","sect_R_item",            set_sect_R_item         },

	{"hecmwST_material","n_mat",                 set_n_mat               },
	{"hecmwST_material","n_mat_item",            set_n_mat_item          },
	{"hecmwST_material","n_mat_subitem",         set_n_mat_subitem       },
	{"hecmwST_material","n_mat_table",           set_n_mat_table         },
	{"hecmwST_material","mat_name",              set_mat_name            },
	{"hecmwST_material","mat_item_index",        set_mat_item_index      },
	{"hecmwST_material","mat_subitem_index",     set_mat_subitem_index   },
	{"hecmwST_material","mat_table_index",       set_mat_table_index     },
	{"hecmwST_material","mat_val",               set_mat_val             },
	{"hecmwST_material","mat_temp",              set_mat_temp            },

	{"hecmwST_mpc","n_mpc",                      set_n_mpc               },
	{"hecmwST_mpc","mpc_index",                  set_mpc_index           },
	{"hecmwST_mpc","mpc_item",                   set_mpc_item            },
	{"hecmwST_mpc","mpc_dof",                    set_mpc_dof             },
	{"hecmwST_mpc","mpc_val",                    set_mpc_val             },
	{"hecmwST_mpc","mpc_const",                  set_mpc_const           },

	{"hecmwST_amplitude","n_amp",                set_n_amp               },
	{"hecmwST_amplitude","amp_name",             set_amp_name            },
	{"hecmwST_amplitude","amp_type_definition",  set_amp_type_definition },
	{"hecmwST_amplitude","amp_type_time",        set_amp_type_time       },
	{"hecmwST_amplitude","amp_type_value",       set_amp_type_value      },
	{"hecmwST_amplitude","amp_index",            set_amp_index           },
	{"hecmwST_amplitude","amp_val",              set_amp_val             },
	{"hecmwST_amplitude","amp_table",            set_amp_table           },

	{"hecmwST_node_grp","n_grp",                 set_ngrp_n_grp          },
	{"hecmwST_node_grp","grp_name",              set_ngrp_grp_name       },
	{"hecmwST_node_grp","grp_index",             set_ngrp_grp_index      },
	{"hecmwST_node_grp","grp_item",              set_ngrp_grp_item       },
	{"hecmwST_node_grp","n_bc",                  set_ngrp_n_bc           },
	{"hecmwST_node_grp","bc_grp_ID",             set_ngrp_bc_grp_ID      },
	{"hecmwST_node_grp","bc_grp_type",           set_ngrp_bc_grp_type    },
	{"hecmwST_node_grp","bc_grp_index",          set_ngrp_bc_grp_index   },
	{"hecmwST_node_grp","bc_grp_dof",            set_ngrp_bc_grp_dof     },
	{"hecmwST_node_grp","bc_grp_val",            set_ngrp_bc_grp_val     },

	{"hecmwST_elem_grp","n_grp",                 set_egrp_n_grp          },
	{"hecmwST_elem_grp","grp_name",              set_egrp_grp_name       },
	{"hecmwST_elem_grp","grp_index",             set_egrp_grp_index      },
	{"hecmwST_elem_grp","grp_item",              set_egrp_grp_item       },
	{"hecmwST_elem_grp","n_bc",                  set_egrp_n_bc           },
	{"hecmwST_elem_grp","bc_grp_ID",             set_egrp_bc_grp_ID      },
	{"hecmwST_elem_grp","bc_grp_type",           set_egrp_bc_grp_type    },
	{"hecmwST_elem_grp","bc_grp_index",          set_egrp_bc_grp_index   },
	{"hecmwST_elem_grp","bc_grp_val",            set_egrp_bc_grp_val     },

	{"hecmwST_surf_grp","n_grp",                 set_sgrp_n_grp          },
	{"hecmwST_surf_grp","grp_name",              set_sgrp_grp_name       },
	{"hecmwST_surf_grp","grp_index",             set_sgrp_grp_index      },
	{"hecmwST_surf_grp","grp_item",              set_sgrp_grp_item       },
	{"hecmwST_surf_grp","n_bc",                  set_sgrp_n_bc           },
	{"hecmwST_surf_grp","bc_grp_ID",             set_sgrp_bc_grp_ID      },
	{"hecmwST_surf_grp","bc_grp_type",           set_sgrp_bc_grp_type    },
	{"hecmwST_surf_grp","bc_grp_index",          set_sgrp_bc_grp_index   },
	{"hecmwST_surf_grp","bc_grp_val",            set_sgrp_bc_grp_val     },

	{"hecmwST_contact_pair","n_pair",            set_contact_pair_n_pair },
	{"hecmwST_contact_pair","name",              set_contact_pair_name},
	{"hecmwST_contact_pair","type",              set_contact_pair_type},
	{"hecmwST_contact_pair","slave_grp_id",      set_contact_pair_slave_grp_id},
	{"hecmwST_contact_pair","master_grp_id",     set_contact_pair_master_grp_id},
};

static const int NFUNC = sizeof(functions) / sizeof(functions[0]);



static SetFunc
get_set_func(char *struct_name, char *var_name)
{
	int i;

	for(i=0; i < NFUNC; i++) {
		if(strcmp(functions[i].struct_name, struct_name) == 0
		&& strcmp(functions[i].var_name, var_name) == 0) {
			return functions[i].set_func;
		}
	}
	return NULL;
}




int
HECMW_dist_copy_f2c_init(struct hecmwST_local_mesh *local_mesh)
{
	struct hecmwST_local_mesh tmp;

	tmp = *local_mesh;
	memset(local_mesh, 0, sizeof(*local_mesh));
	local_mesh->section = tmp.section;
	memset(local_mesh->section, 0, sizeof(*local_mesh->section));
	local_mesh->material = tmp.material;
	memset(local_mesh->material, 0, sizeof(*local_mesh->material));
	local_mesh->mpc = tmp.mpc;
	memset(local_mesh->mpc, 0, sizeof(*local_mesh->mpc));
	local_mesh->amp = tmp.amp;
	memset(local_mesh->amp, 0, sizeof(*local_mesh->amp));
	local_mesh->node_group = tmp.node_group;
	memset(local_mesh->node_group, 0, sizeof(*local_mesh->node_group));
	local_mesh->elem_group = tmp.elem_group;
	memset(local_mesh->elem_group, 0, sizeof(*local_mesh->elem_group));
	local_mesh->surf_group = tmp.surf_group;
	memset(local_mesh->surf_group, 0, sizeof(*local_mesh->surf_group));
	local_mesh->contact_pair = tmp.contact_pair;
	memset(local_mesh->contact_pair, 0, sizeof(*local_mesh->contact_pair));

	mesh = local_mesh;
	return 0;
}


int
HECMW_dist_copy_f2c_finalize(void)
{
	mesh = NULL;
	return 0;
}


/*----------------------------------------------------------------------------*/


void
hecmw_dist_copy_f2c_set_if(char *struct_name, char *var_name, void *src, int *err, int slen, int vlen)
{
	SetFunc func;
	char sname[HECMW_NAME_LEN+1];
	char vname[HECMW_NAME_LEN+1];

	*err = 1;

	if(mesh == NULL) {
		HECMW_set_error(HECMW_ALL_E0102, "HECMW_dist_copy_f2c_set_if(): 'mesh' has not initialized yet");
		return;
	}
	if(struct_name == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "HECMW_dist_copy_f2c_set_if(): 'sname' is NULL");
		return;
	}
	if(var_name == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "HECMW_dist_copy_f2c_set_if(): 'vname' is NULL");
		return;
	}

	if(HECMW_strcpy_f2c_r(struct_name, slen, sname, sizeof(sname)) == NULL) {
		return;
	}
	if(HECMW_strcpy_f2c_r(var_name, vlen, vname, sizeof(vname)) == NULL) {
		return;
	}

	func = get_set_func(sname, vname);
	if(func == NULL) {
		HECMW_set_error(HECMW_ALL_E0102, "HECMW_dist_copy_f2c_set_if(): SetFunc not found");
		return;
	}

	if((*func)(src)) return;

	*err = 0;
}



void
hecmw_dist_copy_f2c_set_if_(char *struct_name, char *var_name,
									void *src, int *err, int slen, int vlen)
{
	hecmw_dist_copy_f2c_set_if(struct_name, var_name, src, err, slen, vlen);
}



void
hecmw_dist_copy_f2c_set_if__(char *struct_name, char *var_name,
									void *src, int *err, int slen, int vlen)
{
	hecmw_dist_copy_f2c_set_if(struct_name, var_name, src, err, slen, vlen);
}



void
HECMW_DIST_COPY_F2C_SET_IF(char *struct_name, char *var_name,
									void *src, int *err, int slen, int vlen)
{
	hecmw_dist_copy_f2c_set_if(struct_name, var_name, src, err, slen, vlen);
}
