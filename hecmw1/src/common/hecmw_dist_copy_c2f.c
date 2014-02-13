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
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_dist_free.h"
#include "hecmw_io_get_mesh.h"
#include "hecmw_dist_copy_c2f.h"

static struct hecmwST_local_mesh *mesh;


/*-----------------------------------------------------------------------------
 * SetFunc
 */

static int
set_hecmw_flag_adapt(void *dst)
{
	void *src;
	int size;

	src = &mesh->hecmw_flag_adapt;
	size = sizeof(mesh->hecmw_flag_adapt);
	memcpy(dst, src, size);

	return 0;
}


static int
set_hecmw_flag_initcon(void *dst)
{
	void *src;
	int size;

	src = &mesh->hecmw_flag_initcon;
	size = sizeof(mesh->hecmw_flag_initcon);
	memcpy(dst, src, size);

	return 0;
}


static int
set_hecmw_flag_parttype(void *dst)
{
	void *src;
	int size;

	src = &mesh->hecmw_flag_parttype;
	size = sizeof(mesh->hecmw_flag_parttype);
	memcpy(dst, src, size);

	return 0;
}


static int
set_hecmw_flag_partdepth(void *dst)
{
	void *src;
	int size;

	src = &mesh->hecmw_flag_partdepth;
	size = sizeof(mesh->hecmw_flag_partdepth);
	memcpy(dst, src, size);
	return 0;
}


static int
set_hecmw_flag_version(void *dst)
{
	void *src;
	int size;

	src = &mesh->hecmw_flag_version;
	size = sizeof(mesh->hecmw_flag_version);
	memcpy(dst, src, size);

	return 0;
}


static int
set_gridfile(void *dst)
{
	void *src;

	src = mesh->gridfile;
	HECMW_strcpy_c2f(src, dst, HECMW_FILENAME_LEN);

	return 0;
}


static int
set_hecmw_n_file(void *dst)
{
	void *src;
	int size;

	src = &mesh->hecmw_n_file;
	size = sizeof(mesh->hecmw_n_file);
	memcpy(dst, src, size);

	return 0;
}


static int
set_files(void *dst)
{
	int i;

	if(mesh->hecmw_n_file <= 0) return 0;

	for(i=0; i < mesh->hecmw_n_file; i++) {
		char *dst_point = (char *)dst + HECMW_FILENAME_LEN * i;
		char *src = mesh->files[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_FILENAME_LEN);
	}

	return 0;
}


static int
set_header(void *dst)
{
	void *src;

	src = mesh->header;
	HECMW_strcpy_c2f(src, dst, HECMW_HEADER_LEN);

	return 0;
}


static int
set_zero_temp(void *dst)
{
	void *src;
	int size;

	src = &mesh->zero_temp;
	size = sizeof(mesh->zero_temp);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_node(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_node;
	size = sizeof(mesh->n_node);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_node_gross(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_node_gross;
	size = sizeof(mesh->n_node_gross);
	memcpy(dst, src, size);

	return 0;
}


static int
set_nn_internal(void *dst)
{
	void *src;
	int size;

	src = &mesh->nn_internal;
	size = sizeof(mesh->nn_internal);
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_internal_list(void *dst)
{
	void *src;
	int size;

	if(mesh->nn_internal <= 0) return 0;

	src = mesh->node_internal_list;
	size = sizeof(*mesh->node_internal_list) * mesh->nn_internal;
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_ID(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_ID;
	size = sizeof(*mesh->node_ID) * mesh->n_node_gross * 2;
	memcpy(dst, src, size);

	return 0;
}


static int
set_global_node_ID(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->global_node_ID;
	size = sizeof(*mesh->global_node_ID) * mesh->n_node_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_node(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node;
	size = sizeof(*mesh->node) * mesh->n_node_gross * 3;
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_dof(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_dof;
	size = sizeof(mesh->n_dof);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_dof_grp(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_dof_grp;
	size = sizeof(mesh->n_dof_grp);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_dof_tot(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_dof_tot;
	size = sizeof(mesh->n_dof_tot);
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_dof_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_dof_grp <= 0) return 0;

	src = mesh->node_dof_index;
	size = sizeof(*mesh->node_dof_index) * (mesh->n_dof_grp+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_dof_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_dof_grp <= 0) return 0;

	src = mesh->node_dof_item;
	size = sizeof(*mesh->node_dof_item) * mesh->n_dof_grp;
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_val_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_val_index;
	size = sizeof(*mesh->node_val_index) * (mesh->n_node_gross+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_val_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_val_item;
	size = sizeof(*mesh->node_val_item) * mesh->node_val_index[mesh->n_node_gross];
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_init_val_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_init_val_index;
	size = sizeof(*mesh->node_init_val_index) * (mesh->n_node_gross+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_init_val_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_init_val_item;
	size = sizeof(*mesh->node_init_val_item) *
					mesh->node_init_val_index[mesh->n_node_gross];
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_elem(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_elem;
	size = sizeof(mesh->n_elem);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_elem_gross(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_elem_gross;
	size = sizeof(mesh->n_elem_gross);
	memcpy(dst, src, size);

	return 0;
}


static int
set_ne_internal(void *dst)
{
	void *src;
	int size;

	src = &mesh->ne_internal;
	size = sizeof(mesh->ne_internal);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_internal_list(void *dst)
{
	void *src;
	int size;

	if(mesh->ne_internal <= 0) return 0;

	src = mesh->elem_internal_list;
	size = sizeof(*mesh->elem_internal_list) * mesh->ne_internal;
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_ID(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_ID;
	size = sizeof(*mesh->elem_ID) * mesh->n_elem_gross * 2;
	memcpy(dst, src, size);

	return 0;
}


static int
set_global_elem_ID(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->global_elem_ID;
	size = sizeof(*mesh->global_elem_ID) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_type(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_type;
	size = sizeof(*mesh->elem_type) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_elem_type(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_elem_type;
	size = sizeof(mesh->n_elem_type);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_type_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_type <= 0) return 0;

	src = mesh->elem_type_index;
	size = sizeof(*mesh->elem_type_index) * (mesh->n_elem_type+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_type_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_type <= 0) return 0;

	src = mesh->elem_type_item;
	size = sizeof(*mesh->elem_type_item) * mesh->n_elem_type;
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_node_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_node_index;
	size = sizeof(*mesh->elem_node_index) * (mesh->n_elem_gross+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_node_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_node_item;
	size = sizeof(*mesh->elem_node_item) *
					mesh->elem_node_index[mesh->n_elem_gross];
	memcpy(dst, src, size);

	return 0;
}


static int
set_section_ID(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->section_ID;
	size = sizeof(*mesh->section_ID) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_elem_mat_ID(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_elem_mat_ID;
	size = sizeof(mesh->n_elem_mat_ID);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_mat_ID_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_mat_ID_index;
	size = sizeof(*mesh->elem_mat_ID_index) * (mesh->n_elem_gross+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_mat_ID_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_mat_ID_item;
	size = sizeof(*mesh->elem_mat_ID_item) *
					mesh->elem_mat_ID_index[mesh->n_elem_gross];
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_mat_int_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_mat_ID <= 0) return 0;

	src = mesh->elem_mat_int_index;
	size = sizeof(*mesh->elem_mat_int_index) * (mesh->n_elem_mat_ID+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_mat_int_val(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_mat_ID <= 0) return 0;

	src = mesh->elem_mat_ID_item;
	size = sizeof(*mesh->elem_mat_int_val) *
					mesh->elem_mat_int_index[mesh->n_elem_mat_ID];
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_val_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_val_index;
	size = sizeof(*mesh->elem_val_index) * (mesh->n_elem_gross+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_val_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_val_item;
	size = sizeof(*mesh->elem_val_item) *
					mesh->elem_val_index[mesh->n_elem_gross];
	memcpy(dst, src, size);

	return 0;
}


static int
set_zero(void *dst)
{
	void *src;
	int size;

	src = &mesh->zero;
	size = sizeof(mesh->zero);
	memcpy(dst, src, size);

	return 0;
}


static int
set_HECMW_COMM(void *dst)
{
	HECMW_Fint comm;
	void *src;
	int size;

	comm = HECMW_Comm_c2f(mesh->HECMW_COMM);
	src = &comm;
	size = sizeof(comm);
	memcpy(dst, src, size);

	return 0;
}


static int
set_PETOT(void *dst)
{
	void *src;
	int size;

	src = &mesh->PETOT;
	size = sizeof(mesh->PETOT);
	memcpy(dst, src, size);

	return 0;
}


static int
set_PEsmpTOT(void *dst)
{
	void *src;
	int size;

	src = &mesh->PEsmpTOT;
	size = sizeof(mesh->PEsmpTOT);
	memcpy(dst, src, size);

	return 0;
}


static int
set_my_rank(void *dst)
{
	void *src;
	int size;

	src = &mesh->my_rank;
	size = sizeof(mesh->my_rank);
	memcpy(dst, src, size);

	return 0;
}


static int
set_errnof(void *dst)
{
	void *src;
	int size;

	src = &mesh->errnof;
	size = sizeof(mesh->errnof);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_subdomain(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_subdomain;
	size = sizeof(mesh->n_subdomain);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_neighbor_pe(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_neighbor_pe;
	size = sizeof(mesh->n_neighbor_pe);
	memcpy(dst, src, size);

	return 0;
}


static int
set_neighbor_pe(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->neighbor_pe;
	size = sizeof(*mesh->neighbor_pe) * mesh->n_neighbor_pe;
	memcpy(dst, src, size);

	return 0;
}

static int
set_import_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->import_index;
	size = sizeof(*mesh->import_index) * (mesh->n_neighbor_pe+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_import_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->import_item;
	size = sizeof(*mesh->import_item) *
					mesh->import_index[mesh->n_neighbor_pe];
	memcpy(dst, src, size);

	return 0;
}


static int
set_export_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->export_index;
	size = sizeof(*mesh->export_index) * (mesh->n_neighbor_pe+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_export_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->export_item;
	size = sizeof(*mesh->export_item) *
					mesh->export_index[mesh->n_neighbor_pe];
	memcpy(dst, src, size);

	return 0;
}


static int
set_shared_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->shared_index;
	size = sizeof(*mesh->shared_index) * (mesh->n_neighbor_pe+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_shared_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_neighbor_pe <= 0) return 0;

	src = mesh->shared_item;
	size = sizeof(*mesh->shared_item) *
					mesh->shared_index[mesh->n_neighbor_pe];
	memcpy(dst, src, size);

	return 0;
}


static int
set_coarse_grid_level(void *dst)
{
	void *src;
	int size;

	src = &mesh->coarse_grid_level;
	size = sizeof(mesh->coarse_grid_level);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_adapt(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_adapt;
	size = sizeof(mesh->n_adapt);
	memcpy(dst, src, size);

	return 0;
}


static int
set_when_i_was_refined_node(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->when_i_was_refined_node;
	size = sizeof(*mesh->when_i_was_refined_node) * mesh->n_node_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_when_i_was_refined_elem(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->when_i_was_refined_elem;
	size = sizeof(*mesh->when_i_was_refined_elem) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_adapt_parent_type(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->adapt_parent_type;
	size = sizeof(*mesh->adapt_parent_type) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_adapt_type(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->adapt_type;
	size = sizeof(*mesh->adapt_type) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_adapt_level(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->adapt_level;
	size = sizeof(*mesh->adapt_level) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_adapt_parent(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->adapt_parent;
	size = sizeof(*mesh->adapt_parent) * mesh->n_elem_gross * 2;
	memcpy(dst, src, size);

	return 0;
}


static int
set_adapt_children_index(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->adapt_children_index;
	size = sizeof(*mesh->adapt_children_index) * (mesh->n_elem_gross+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_adapt_children_item(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->adapt_children_item;
	size = sizeof(*mesh->adapt_children_item) *
					mesh->adapt_children_index[mesh->n_elem_gross] * 2;
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_refine(void *dst)
{
	void *src;
	int size;

	src = &mesh->n_refine;
	size = sizeof(mesh->n_refine);
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_old2new(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_old2new;
	size = sizeof(*mesh->node_old2new) * mesh->n_node_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_node_new2old(void *dst)
{
	void *src;
	int size;

	if(mesh->n_node_gross <= 0) return 0;

	src = mesh->node_new2old;
	size = sizeof(*mesh->node_new2old) * mesh->n_node_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_old2new(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_old2new;
	size = sizeof(*mesh->elem_old2new) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_elem_new2old(void *dst)
{
	void *src;
	int size;

	if(mesh->n_elem_gross <= 0) return 0;

	src = mesh->elem_new2old;
	size = sizeof(*mesh->elem_new2old) * mesh->n_elem_gross;
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_sect(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	src = &sect->n_sect;
	size = sizeof(sect->n_sect);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_type(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_type;
	size = sizeof(*sect->sect_type) * sect->n_sect;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_opt(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_opt;
	size = sizeof(*sect->sect_opt) * sect->n_sect;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_mat_ID_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_mat_ID_index;
	size = sizeof(*sect->sect_mat_ID_index) * (sect->n_sect+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_mat_ID_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_mat_ID_item;
	size = sizeof(*sect->sect_mat_ID_item) *
					sect->sect_mat_ID_index[sect->n_sect];
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_I_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_I_index;
	size = sizeof(*sect->sect_I_index) * (sect->n_sect+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_I_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;
    if( sect->sect_I_item==NULL ) return 0;

	src = sect->sect_I_item;
	size = sizeof(*sect->sect_I_item) * sect->sect_I_index[sect->n_sect];
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_R_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_R_index;
	size = sizeof(*sect->sect_R_index) * (sect->n_sect+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sect_R_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_section *sect = mesh->section;

	if(sect->n_sect <= 0) return 0;

	src = sect->sect_R_item;
	size = sizeof(*sect->sect_R_item) * sect->sect_R_index[sect->n_sect];
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_mat(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	src = &mat->n_mat;
	size = sizeof(mat->n_mat);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_mat_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	src = &mat->n_mat_item;
	size = sizeof(mat->n_mat_item);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_mat_subitem(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	src = &mat->n_mat_subitem;
	size = sizeof(mat->n_mat_subitem);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_mat_table(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	src = &mat->n_mat_table;
	size = sizeof(mat->n_mat_table);
	memcpy(dst, src, size);

	return 0;
}


static int
set_mat_name(void *dst)
{
	int i;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat <= 0) return 0;

	for(i=0; i < mat->n_mat; i++) {
		char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
		char *src = mat->mat_name[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
	}

	return 0;
}


static int
set_mat_item_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat <= 0) return 0;

	src = mat->mat_item_index;
	size = sizeof(*mat->mat_item_index) * (mat->n_mat+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_mat_subitem_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat_item <= 0) return 0;

	src = mat->mat_subitem_index;
	size = sizeof(*mat->mat_subitem_index) * (mat->n_mat_item+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_mat_table_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat_subitem <= 0) return 0;

	src = mat->mat_table_index;
	size = sizeof(*mat->mat_table_index) * (mat->n_mat_subitem+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_mat_val(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat <= 0) return 0;

	src = mat->mat_val;
	size = sizeof(*mat->mat_val) * mat->mat_table_index[mat->n_mat_subitem];
	memcpy(dst, src, size);

	return 0;
}


static int
set_mat_temp(void *dst)
{
	void *src;
	int size;
	struct hecmwST_material *mat = mesh->material;

	if(mat->n_mat <= 0) return 0;

	src = mat->mat_temp;
	size = sizeof(*mat->mat_temp) * mat->mat_table_index[mat->n_mat_subitem];
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_mpc(void *dst)
{
	void *src;
	int size;
	struct hecmwST_mpc *mpc = mesh->mpc;

	src = &mpc->n_mpc;
	size = sizeof(mpc->n_mpc);
	memcpy(dst, src, size);

	return 0;
}


static int
set_mpc_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_mpc *mpc = mesh->mpc;

	if(mpc->n_mpc <= 0) return 0;

	src = mpc->mpc_index;
	size = sizeof(*mpc->mpc_index) * (mpc->n_mpc+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_mpc_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_mpc *mpc = mesh->mpc;

	if(mpc->n_mpc <= 0) return 0;

	src = mpc->mpc_item;
	size = sizeof(*mpc->mpc_item) * mpc->mpc_index[mpc->n_mpc];
	memcpy(dst, src, size);

	return 0;
}


static int
set_mpc_dof(void *dst)
{
	void *src;
	int size;
	struct hecmwST_mpc *mpc = mesh->mpc;

	if(mpc->n_mpc <= 0) return 0;

	src = mpc->mpc_dof;
	size = sizeof(*mpc->mpc_dof) * mpc->mpc_index[mpc->n_mpc];
	memcpy(dst, src, size);

	return 0;
}


static int
set_mpc_val(void *dst)
{
	void *src;
	int size;
	struct hecmwST_mpc *mpc = mesh->mpc;

	if(mpc->n_mpc <= 0) return 0;

	src = mpc->mpc_val;
	size = sizeof(*mpc->mpc_val) * mpc->mpc_index[mpc->n_mpc];
	memcpy(dst, src, size);

	return 0;
}


static int
set_mpc_const(void *dst)
{
	void *src;
	int size;
	struct hecmwST_mpc *mpc = mesh->mpc;

	if(mpc->n_mpc <= 0) return 0;

	src = mpc->mpc_const;
	size = sizeof(*mpc->mpc_const) * mpc->n_mpc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_amp(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	src = &amp->n_amp;
	size = sizeof(amp->n_amp);
	memcpy(dst, src, size);

	return 0;
}


static int
set_amp_name(void *dst)
{
	int i;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	for(i=0; i < amp->n_amp; i++) {
		char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
		char *src = amp->amp_name[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
	}

	return 0;
}


static int
set_amp_type_definition(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	src = amp->amp_type_definition;
	size = sizeof(*amp->amp_type_definition) * amp->n_amp;
	memcpy(dst, src, size);

	return 0;
}


static int
set_amp_type_time(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	src = amp->amp_type_time;
	size = sizeof(*amp->amp_type_time) * amp->n_amp;
	memcpy(dst, src, size);

	return 0;
}


static int
set_amp_type_value(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	src = amp->amp_type_value;
	size = sizeof(*amp->amp_type_value) * amp->n_amp;
	memcpy(dst, src, size);

	return 0;
}


static int
set_amp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	src = amp->amp_index;
	size = sizeof(*amp->amp_index) * (amp->n_amp+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_amp_val(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	src = amp->amp_val;
	size = sizeof(*amp->amp_val) * amp->amp_index[amp->n_amp];
	memcpy(dst, src, size);

	return 0;
}


static int
set_amp_table(void *dst)
{
	void *src;
	int size;
	struct hecmwST_amplitude *amp = mesh->amp;

	if(amp->n_amp <= 0) return 0;

	src = amp->amp_table;
	size = sizeof(*amp->amp_table) * amp->amp_index[amp->n_amp];
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_n_grp(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	src = &grp->n_grp;
	size = sizeof(grp->n_grp);
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_n_bc(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	src = &grp->n_bc;
	size = sizeof(grp->n_bc);
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_grp_name(void *dst)
{
	int i;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_grp <= 0) return 0;

	for(i=0; i < grp->n_grp; i++) {
		char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
		char *src = grp->grp_name[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
	}

	return 0;
}


static int
set_ngrp_grp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_grp <= 0) return 0;

	src = grp->grp_index;
	size = sizeof(*grp->grp_index) * (grp->n_grp+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_grp_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_grp <= 0) return 0;

	src = grp->grp_item;
	size = sizeof(*grp->grp_item) * grp->grp_index[grp->n_grp];
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_bc_grp_ID(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_ID;
	size = sizeof(*grp->bc_grp_ID) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_bc_grp_type(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_type;
	size = sizeof(*grp->bc_grp_type) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_bc_grp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_index;
	size = sizeof(*grp->bc_grp_index) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_bc_grp_dof(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_dof;
	size = sizeof(*grp->bc_grp_dof) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_ngrp_bc_grp_val(void *dst)
{
	void *src;
	int size;
	struct hecmwST_node_grp *grp = mesh->node_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_val;
	size = sizeof(*grp->bc_grp_val) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_n_grp(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	src = &grp->n_grp;
	size = sizeof(grp->n_grp);
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_n_bc(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	src = &grp->n_bc;
	size = sizeof(grp->n_bc);
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_grp_name(void *dst)
{
	int i;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_grp <= 0) return 0;

	for(i=0; i < grp->n_grp; i++) {
		char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
		char *src = grp->grp_name[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
	}

	return 0;
}


static int
set_egrp_grp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_grp <= 0) return 0;

	src = grp->grp_index;
	size = sizeof(*grp->grp_index) * (grp->n_grp+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_grp_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_grp <= 0) return 0;

	src = grp->grp_item;
	size = sizeof(*grp->grp_item) * grp->grp_index[grp->n_grp];
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_bc_grp_ID(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_ID;
	size = sizeof(*grp->bc_grp_ID) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_bc_grp_type(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_type;
	size = sizeof(*grp->bc_grp_type) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_bc_grp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_index;
	size = sizeof(*grp->bc_grp_index) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_egrp_bc_grp_val(void *dst)
{
	void *src;
	int size;
	struct hecmwST_elem_grp *grp = mesh->elem_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_val;
	size = sizeof(*grp->bc_grp_val) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_n_grp(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	src = &grp->n_grp;
	size = sizeof(grp->n_grp);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_n_bc(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	src = &grp->n_bc;
	size = sizeof(grp->n_bc);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_grp_name(void *dst)
{
	int i;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_grp <= 0) return 0;

	for(i=0; i < grp->n_grp; i++) {
		char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
		char *src = grp->grp_name[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
	}

	return 0;
}


static int
set_sgrp_grp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_grp <= 0) return 0;

	src = grp->grp_index;
	size = sizeof(*grp->grp_index) * (grp->n_grp+1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_grp_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_grp <= 0) return 0;

	src = grp->grp_item;
	size = sizeof(*grp->grp_item) * grp->grp_index[grp->n_grp] * 2;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_bc_grp_ID(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_ID;
	size = sizeof(*grp->bc_grp_ID) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_bc_grp_type(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_type;
	size = sizeof(*grp->bc_grp_type) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_bc_grp_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_index;
	size = sizeof(*grp->bc_grp_index) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_sgrp_bc_grp_val(void *dst)
{
	void *src;
	int size;
	struct hecmwST_surf_grp *grp = mesh->surf_group;

	if(grp->n_bc <= 0) return 0;

	src = grp->bc_grp_val;
	size = sizeof(*grp->bc_grp_val) * grp->n_bc;
	memcpy(dst, src, size);

	return 0;
}


static int
set_contact_pair_n_pair(void *dst)
{
	void *src;
	int size;
	struct hecmwST_contact_pair *cpair = mesh->contact_pair;

	src = &cpair->n_pair;
	size = sizeof(cpair->n_pair);
	memcpy(dst, src, size);

	return 0;
}


static int
set_contact_pair_name(void *dst)
{
	int i;
	struct hecmwST_contact_pair *cpair = mesh->contact_pair;

	if(cpair->n_pair <= 0) return 0;

	for(i=0; i < cpair->n_pair; i++) {
		char *dst_point = (char *)dst + HECMW_NAME_LEN * i;
		char *src = cpair->name[i];
		HECMW_strcpy_c2f(src, dst_point, HECMW_NAME_LEN);
	}

	return 0;
}


static int
set_contact_pair_type(void *dst)
{
	void *src;
	int size;
	struct hecmwST_contact_pair *cpair = mesh->contact_pair;

	if(cpair->n_pair <= 0) return 0;

	src = cpair->type;
	size = sizeof(*cpair->type) * (cpair->n_pair);
	memcpy(dst, src, size);

	return 0;
}


static int
set_contact_pair_slave_grp_id(void *dst)
{
	void *src;
	int size;
	struct hecmwST_contact_pair *cpair = mesh->contact_pair;

	if(cpair->n_pair <= 0) return 0;

	src = cpair->slave_grp_id;
	size = sizeof(*cpair->slave_grp_id) * (cpair->n_pair);
	memcpy(dst, src, size);

	return 0;
}


static int
set_contact_pair_master_grp_id(void *dst)
{
	void *src;
	int size;
	struct hecmwST_contact_pair *cpair = mesh->contact_pair;

	if(cpair->n_pair <= 0) return 0;

	src = cpair->master_grp_id;
	size = sizeof(*cpair->master_grp_id) * (cpair->n_pair);
	memcpy(dst, src, size);

	return 0;
}


static int
set_refine_origin_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_refine_origin *reforg = mesh->refine_origin;

	src = reforg->index;
	size = sizeof(*reforg->index) * (mesh->n_refine + 1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_refine_origin_item_index(void *dst)
{
	void *src;
	int size;
	struct hecmwST_refine_origin *reforg = mesh->refine_origin;

	src = reforg->item_index;
	size = sizeof(*reforg->item_index) * (reforg->index[mesh->n_refine] + 1);
	memcpy(dst, src, size);

	return 0;
}


static int
set_refine_origin_item_item(void *dst)
{
	void *src;
	int size;
	struct hecmwST_refine_origin *reforg = mesh->refine_origin;

	src = reforg->item_item;
	size = sizeof(*reforg->item_item) * (reforg->item_index[reforg->index[mesh->n_refine]]);
	memcpy(dst, src, size);

	return 0;
}


/*---------------------------------------------------------------------------*/

static int
is_alloc_node_internal_list(void)
{
	return mesh->node_internal_list ? 1 : 0;
}


static int
is_alloc_node_val_index(void)
{
	return mesh->node_val_index ? 1 : 0;
}


static int
is_alloc_node_val_item(void)
{
	return mesh->node_val_item ? 1 : 0;
}


static int
is_alloc_node_init_val_index(void)
{
	return mesh->node_init_val_index ? 1 : 0;
}


static int
is_alloc_node_init_val_item(void)
{
	return mesh->node_init_val_item ? 1 : 0;
}


static int
is_alloc_elem_internal_list(void)
{
	return mesh->elem_internal_list ? 1 : 0;
}


static int
is_alloc_elem_mat_int_index(void)
{
	return mesh->elem_mat_int_index ? 1 : 0;
}


static int
is_alloc_elem_mat_int_val(void)
{
	return mesh->elem_mat_int_val ? 1 : 0;
}


static int
is_alloc_elem_val_index(void)
{
	return mesh->elem_val_index ? 1 : 0;
}


static int
is_alloc_elem_val_item(void)
{
	return mesh->elem_val_item ? 1 : 0;
}


static int
is_alloc_node_old2new(void)
{
	return mesh->node_old2new ? 1 : 0;
}


static int
is_alloc_node_new2old(void)
{
	return mesh->node_new2old ? 1 : 0;
}


static int
is_alloc_elem_old2new(void)
{
	return mesh->elem_old2new ? 1 : 0;
}


static int
is_alloc_elem_new2old(void)
{
	return mesh->elem_new2old ? 1 : 0;
}


static int
is_alloc_when_i_was_refined_node(void)
{
	return mesh->when_i_was_refined_node ? 1 : 0;
}


static int
is_alloc_when_i_was_refined_elem(void)
{
	return mesh->when_i_was_refined_elem ? 1 : 0;
}


static int
is_alloc_adapt_parent_type(void)
{
	return mesh->adapt_parent_type ? 1 : 0;
}


static int
is_alloc_adapt_type(void)
{
	return mesh->adapt_type ? 1 : 0;
}


static int
is_alloc_adapt_level(void)
{
	return mesh->adapt_level ? 1 : 0;
}


static int
is_alloc_adapt_parent(void)
{
	return mesh->adapt_parent ? 1 : 0;
}


static int
is_alloc_adapt_children_index(void)
{
	return mesh->adapt_children_index ? 1 : 0;
}


static int
is_alloc_adapt_children_item(void)
{
	return mesh->adapt_children_item ? 1 : 0;
}


static int
is_alloc_ngrp_bc_grp_ID(void)
{
	return mesh->node_group->bc_grp_ID ? 1 : 0;
}


static int
is_alloc_ngrp_bc_grp_type(void)
{
	return mesh->node_group->bc_grp_type ? 1 : 0;
}


static int
is_alloc_ngrp_bc_grp_index(void)
{
	return mesh->node_group->bc_grp_index ? 1 : 0;
}


static int
is_alloc_ngrp_bc_grp_dof(void)
{
	return mesh->node_group->bc_grp_dof ? 1 : 0;
}


static int
is_alloc_ngrp_bc_grp_val(void)
{
	return mesh->node_group->bc_grp_val ? 1 : 0;
}


static int
is_alloc_egrp_bc_grp_ID(void)
{
	return mesh->elem_group->bc_grp_ID ? 1 : 0;
}


static int
is_alloc_egrp_bc_grp_type(void)
{
	return mesh->elem_group->bc_grp_type ? 1 : 0;
}


static int
is_alloc_egrp_bc_grp_index(void)
{
	return mesh->elem_group->bc_grp_index ? 1 : 0;
}


static int
is_alloc_egrp_bc_grp_val(void)
{
	return mesh->elem_group->bc_grp_val ? 1 : 0;
}


static int
is_alloc_sgrp_bc_grp_ID(void)
{
	return mesh->surf_group->bc_grp_ID ? 1 : 0;
}


static int
is_alloc_sgrp_bc_grp_type(void)
{
	return mesh->surf_group->bc_grp_type ? 1 : 0;
}


static int
is_alloc_sgrp_bc_grp_index(void)
{
	return mesh->surf_group->bc_grp_index ? 1 : 0;
}


static int
is_alloc_sgrp_bc_grp_val(void)
{
	return mesh->surf_group->bc_grp_val ? 1 : 0;
}


/*-----------------------------------------------------------------------------
 * SetFunc table
 */

typedef int (*SetFunc)(void *);
typedef int (*IsAllocatedFunc)(void);

static struct func_table {
	char *struct_name;
	char *var_name;
	SetFunc set_func;
	IsAllocatedFunc is_allocated_func;
} functions[] = {
/*  { Struct name, Variable name, memcpy function, check allocation function }*/
	{"hecmwST_local_mesh","hecmw_flag_adapt",    set_hecmw_flag_adapt,    NULL},
	{"hecmwST_local_mesh","hecmw_flag_initcon",  set_hecmw_flag_initcon,  NULL},
	{"hecmwST_local_mesh","hecmw_flag_parttype", set_hecmw_flag_parttype, NULL},
	{"hecmwST_local_mesh","hecmw_flag_partdepth",set_hecmw_flag_partdepth,NULL},
	{"hecmwST_local_mesh","hecmw_flag_version",  set_hecmw_flag_version,  NULL},

	{"hecmwST_local_mesh","gridfile",            set_gridfile,            NULL},
	{"hecmwST_local_mesh","hecmw_n_file",        set_hecmw_n_file,        NULL},
	{"hecmwST_local_mesh","files",               set_files,               NULL},
	{"hecmwST_local_mesh","header",              set_header,              NULL},
	{"hecmwST_local_mesh","zero_temp",           set_zero_temp,           NULL},

	{"hecmwST_local_mesh","n_node",              set_n_node,              NULL},
	{"hecmwST_local_mesh","n_node_gross",        set_n_node_gross,        NULL},
	{"hecmwST_local_mesh","nn_internal",         set_nn_internal,         NULL},
	{"hecmwST_local_mesh","node_internal_list",  set_node_internal_list,
		                                         is_alloc_node_internal_list  },
	{"hecmwST_local_mesh","node_ID",             set_node_ID,             NULL},
	{"hecmwST_local_mesh","global_node_ID",      set_global_node_ID,      NULL},
	{"hecmwST_local_mesh","node",                set_node,                NULL},
	{"hecmwST_local_mesh","n_dof",               set_n_dof,               NULL},
	{"hecmwST_local_mesh","n_dof_grp",           set_n_dof_grp,           NULL},
	{"hecmwST_local_mesh","n_dof_tot",           set_n_dof_tot,           NULL},
	{"hecmwST_local_mesh","node_dof_index",      set_node_dof_index,      NULL},
	{"hecmwST_local_mesh","node_dof_item",       set_node_dof_item,       NULL},
	{"hecmwST_local_mesh","node_val_index",      set_node_val_index,
	                                             is_alloc_node_val_index      },
	{"hecmwST_local_mesh","node_val_item",       set_node_val_item,
	                                             is_alloc_node_val_item       },
	{"hecmwST_local_mesh","node_init_val_index", set_node_init_val_index,
	                                             is_alloc_node_init_val_index },
	{"hecmwST_local_mesh","node_init_val_item",  set_node_init_val_item,
	                                             is_alloc_node_init_val_item  },
	{"hecmwST_local_mesh","n_elem",              set_n_elem,              NULL},
	{"hecmwST_local_mesh","n_elem_gross",        set_n_elem_gross,        NULL},
	{"hecmwST_local_mesh","ne_internal",         set_ne_internal,         NULL},
	{"hecmwST_local_mesh","elem_internal_list",  set_elem_internal_list,
	                                             is_alloc_elem_internal_list  },
	{"hecmwST_local_mesh","elem_ID",             set_elem_ID,             NULL},
	{"hecmwST_local_mesh","global_elem_ID",      set_global_elem_ID,      NULL},
	{"hecmwST_local_mesh","elem_type",           set_elem_type,           NULL},
	{"hecmwST_local_mesh","n_elem_type",         set_n_elem_type,         NULL},
	{"hecmwST_local_mesh","elem_type_index",     set_elem_type_index,     NULL},
	{"hecmwST_local_mesh","elem_type_item",      set_elem_type_item,      NULL},
	{"hecmwST_local_mesh","elem_node_index",     set_elem_node_index,     NULL},
	{"hecmwST_local_mesh","elem_node_item",      set_elem_node_item,      NULL},
	{"hecmwST_local_mesh","section_ID",          set_section_ID,          NULL},
	{"hecmwST_local_mesh","n_elem_mat_ID",       set_n_elem_mat_ID,       NULL},
	{"hecmwST_local_mesh","elem_mat_ID_index",   set_elem_mat_ID_index,   NULL},
	{"hecmwST_local_mesh","elem_mat_ID_item",    set_elem_mat_ID_item,    NULL},
	{"hecmwST_local_mesh","elem_mat_int_index",  set_elem_mat_int_index,
	                                             is_alloc_elem_mat_int_index  },
	{"hecmwST_local_mesh","elem_mat_int_val",    set_elem_mat_int_val,
	                                             is_alloc_elem_mat_int_val    },
	{"hecmwST_local_mesh","elem_val_index",      set_elem_val_index,
	                                             is_alloc_elem_val_index      },
	{"hecmwST_local_mesh","elem_val_item",       set_elem_val_item,
	                                             is_alloc_elem_val_item       },

	{"hecmwST_local_mesh","zero",                set_zero,                NULL},
	{"hecmwST_local_mesh","HECMW_COMM",          set_HECMW_COMM,          NULL},
	{"hecmwST_local_mesh","PETOT",               set_PETOT,               NULL},
	{"hecmwST_local_mesh","PEsmpTOT",            set_PEsmpTOT,            NULL},
	{"hecmwST_local_mesh","my_rank",             set_my_rank,             NULL},
	{"hecmwST_local_mesh","errnof",              set_errnof,              NULL},
	{"hecmwST_local_mesh","n_subdomain",         set_n_subdomain,         NULL},
	{"hecmwST_local_mesh","n_neighbor_pe",       set_n_neighbor_pe,       NULL},
	{"hecmwST_local_mesh","neighbor_pe",         set_neighbor_pe,         NULL},
	{"hecmwST_local_mesh","import_index",        set_import_index,        NULL},
	{"hecmwST_local_mesh","import_item",         set_import_item,         NULL},
	{"hecmwST_local_mesh","export_index",        set_export_index,        NULL},
	{"hecmwST_local_mesh","export_item",         set_export_item,         NULL},
	{"hecmwST_local_mesh","shared_index",        set_shared_index,        NULL},
	{"hecmwST_local_mesh","shared_item",         set_shared_item,         NULL},

	{"hecmwST_local_mesh","coarse_grid_level",   set_coarse_grid_level,   NULL},
	{"hecmwST_local_mesh","n_adapt",             set_n_adapt,             NULL},
	{"hecmwST_local_mesh","when_i_was_refined_node",set_when_i_was_refined_node,
	                                                is_alloc_when_i_was_refined_node},
	{"hecmwST_local_mesh","when_i_was_refined_elem",set_when_i_was_refined_elem,
	                                                is_alloc_when_i_was_refined_elem},
	{"hecmwST_local_mesh","adapt_parent_type",   set_adapt_parent_type,
	                                             is_alloc_adapt_parent_type   },
	{"hecmwST_local_mesh","adapt_type",          set_adapt_type,
	                                             is_alloc_adapt_type          },
	{"hecmwST_local_mesh","adapt_level",         set_adapt_level,
	                                             is_alloc_adapt_level         },
	{"hecmwST_local_mesh","adapt_parent",        set_adapt_parent,
	                                             is_alloc_adapt_parent        },
	{"hecmwST_local_mesh","adapt_children_index",set_adapt_children_index,
	                                             is_alloc_adapt_children_index},
	{"hecmwST_local_mesh","adapt_children_item", set_adapt_children_item,
	                                             is_alloc_adapt_children_item},

	{"hecmwST_local_mesh","n_refine",            set_n_refine,            NULL},
	{"hecmwST_local_mesh","node_old2new",        set_node_old2new,
	                                             is_alloc_node_old2new        },
	{"hecmwST_local_mesh","node_new2old",        set_node_new2old,
	                                             is_alloc_node_new2old        },
	{"hecmwST_local_mesh","elem_old2new",        set_elem_old2new,
	                                             is_alloc_elem_old2new        },
	{"hecmwST_local_mesh","elem_new2old",        set_elem_new2old,
	                                             is_alloc_elem_new2old        },

	{"hecmwST_section","n_sect",                 set_n_sect,              NULL},
	{"hecmwST_section","sect_type",              set_sect_type,           NULL},
	{"hecmwST_section","sect_opt",               set_sect_opt,            NULL},
	{"hecmwST_section","sect_mat_ID_index",      set_sect_mat_ID_index,   NULL},
	{"hecmwST_section","sect_mat_ID_item",       set_sect_mat_ID_item,    NULL},
	{"hecmwST_section","sect_I_index",           set_sect_I_index,        NULL},
	{"hecmwST_section","sect_I_item",            set_sect_I_item,         NULL},
	{"hecmwST_section","sect_R_index",           set_sect_R_index,        NULL},
	{"hecmwST_section","sect_R_item",            set_sect_R_item,         NULL},

	{"hecmwST_material","n_mat",                 set_n_mat,               NULL},
	{"hecmwST_material","n_mat_item",            set_n_mat_item,          NULL},
	{"hecmwST_material","n_mat_subitem",         set_n_mat_subitem,       NULL},
	{"hecmwST_material","n_mat_table",           set_n_mat_table,         NULL},
	{"hecmwST_material","mat_name",              set_mat_name,            NULL},
	{"hecmwST_material","mat_item_index",        set_mat_item_index,      NULL},
	{"hecmwST_material","mat_subitem_index",     set_mat_subitem_index,   NULL},
	{"hecmwST_material","mat_table_index",       set_mat_table_index,     NULL},
	{"hecmwST_material","mat_val",               set_mat_val,             NULL},
	{"hecmwST_material","mat_temp",              set_mat_temp,            NULL},

	{"hecmwST_mpc","n_mpc",                      set_n_mpc,               NULL},
	{"hecmwST_mpc","mpc_index",                  set_mpc_index,           NULL},
	{"hecmwST_mpc","mpc_item",                   set_mpc_item,            NULL},
	{"hecmwST_mpc","mpc_dof",                    set_mpc_dof,             NULL},
	{"hecmwST_mpc","mpc_val",                    set_mpc_val,             NULL},
	{"hecmwST_mpc","mpc_const",                  set_mpc_const,           NULL},

	{"hecmwST_amplitude","n_amp",                set_n_amp,               NULL},
	{"hecmwST_amplitude","amp_name",             set_amp_name,            NULL},
	{"hecmwST_amplitude","amp_type_definition",  set_amp_type_definition, NULL},
	{"hecmwST_amplitude","amp_type_time",        set_amp_type_time,       NULL},
	{"hecmwST_amplitude","amp_type_value",       set_amp_type_value,      NULL},
	{"hecmwST_amplitude","amp_index",            set_amp_index,           NULL},
	{"hecmwST_amplitude","amp_val",              set_amp_val,             NULL},
	{"hecmwST_amplitude","amp_table",            set_amp_table,           NULL},

	{"hecmwST_node_grp","n_grp",                 set_ngrp_n_grp,          NULL},
	{"hecmwST_node_grp","grp_name",              set_ngrp_grp_name,       NULL},
	{"hecmwST_node_grp","grp_index",             set_ngrp_grp_index,      NULL},
	{"hecmwST_node_grp","grp_item",              set_ngrp_grp_item,       NULL},
	{"hecmwST_node_grp","n_bc",                  set_ngrp_n_bc,           NULL},
	{"hecmwST_node_grp","bc_grp_ID",             set_ngrp_bc_grp_ID,
                                                 is_alloc_ngrp_bc_grp_ID      },
	{"hecmwST_node_grp","bc_grp_type",           set_ngrp_bc_grp_type,
                                                 is_alloc_ngrp_bc_grp_type    },
	{"hecmwST_node_grp","bc_grp_index",          set_ngrp_bc_grp_index,
                                                 is_alloc_ngrp_bc_grp_index   },
	{"hecmwST_node_grp","bc_grp_dof",            set_ngrp_bc_grp_dof,
                                                 is_alloc_ngrp_bc_grp_dof     },
	{"hecmwST_node_grp","bc_grp_val",            set_ngrp_bc_grp_val,
                                                 is_alloc_ngrp_bc_grp_val     },

	{"hecmwST_elem_grp","n_grp",                 set_egrp_n_grp,          NULL},
	{"hecmwST_elem_grp","grp_name",              set_egrp_grp_name,       NULL},
	{"hecmwST_elem_grp","grp_index",             set_egrp_grp_index,      NULL},
	{"hecmwST_elem_grp","grp_item",              set_egrp_grp_item,       NULL},
	{"hecmwST_elem_grp","n_bc",                  set_egrp_n_bc,           NULL},
	{"hecmwST_elem_grp","bc_grp_ID",             set_egrp_bc_grp_ID,
                                                 is_alloc_egrp_bc_grp_ID      },
	{"hecmwST_elem_grp","bc_grp_type",           set_egrp_bc_grp_type,
                                                 is_alloc_egrp_bc_grp_type    },
	{"hecmwST_elem_grp","bc_grp_index",          set_egrp_bc_grp_index,
                                                 is_alloc_egrp_bc_grp_index   },
	{"hecmwST_elem_grp","bc_grp_val",            set_egrp_bc_grp_val,
                                                 is_alloc_egrp_bc_grp_val     },

	{"hecmwST_surf_grp","n_grp",                 set_sgrp_n_grp,          NULL},
	{"hecmwST_surf_grp","grp_name",              set_sgrp_grp_name,       NULL},
	{"hecmwST_surf_grp","grp_index",             set_sgrp_grp_index,      NULL},
	{"hecmwST_surf_grp","grp_item",              set_sgrp_grp_item,       NULL},
	{"hecmwST_surf_grp","n_bc",                  set_sgrp_n_bc,           NULL},
	{"hecmwST_surf_grp","bc_grp_ID",             set_sgrp_bc_grp_ID,
                                                 is_alloc_sgrp_bc_grp_ID      },
	{"hecmwST_surf_grp","bc_grp_type",           set_sgrp_bc_grp_type,
                                                 is_alloc_sgrp_bc_grp_type    },
	{"hecmwST_surf_grp","bc_grp_index",          set_sgrp_bc_grp_index,
                                                 is_alloc_sgrp_bc_grp_index   },
	{"hecmwST_surf_grp","bc_grp_val",            set_sgrp_bc_grp_val,
                                                 is_alloc_sgrp_bc_grp_val     },

	{"hecmwST_contact_pair","n_pair",            set_contact_pair_n_pair, NULL},
	{"hecmwST_contact_pair","name",              set_contact_pair_name,   NULL},
	{"hecmwST_contact_pair","type",              set_contact_pair_type,   NULL},
	{"hecmwST_contact_pair","slave_grp_id",      set_contact_pair_slave_grp_id,NULL},
	{"hecmwST_contact_pair","master_grp_id",     set_contact_pair_master_grp_id,NULL},
	{"hecmwST_refine_origin","index",            set_refine_origin_index, NULL},
	{"hecmwST_refine_origin","item_index",       set_refine_origin_item_index,NULL},
	{"hecmwST_refine_origin","item_item",        set_refine_origin_item_item,NULL},
};

static const int NFUNC = sizeof(functions) / sizeof(functions[0]);



static IsAllocatedFunc
get_is_allocated_func(char *struct_name, char *var_name)
{
	int i;

	for(i=0; i < NFUNC; i++) {
		if(strcmp(functions[i].struct_name, struct_name) == 0
		&& strcmp(functions[i].var_name, var_name) == 0) {
			return functions[i].is_allocated_func;
		}
	}
	return NULL;
}



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


/*----------------------------------------------------------------------------*/

int
HECMW_dist_copy_c2f_init(struct hecmwST_local_mesh *local_mesh)
{
	mesh = local_mesh;
	return 0;
}


int
HECMW_dist_copy_c2f_finalize(void)
{
	mesh = NULL;
	return 0;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/


void
hecmw_dist_copy_c2f_isalloc_if(char *struct_name, char *var_name,
						int *is_allocated, int *err, int len_struct, int len_var)
{
	IsAllocatedFunc func;
	char sname[HECMW_NAME_LEN+1];
	char vname[HECMW_NAME_LEN+1];

	*err = 1;

	if(mesh == NULL) {
		HECMW_set_error(HECMW_ALL_E0102, "mesh copy(C->Fortran): confirm allocation : 'mesh' has not initialized yet");
		return;
	}
	if(struct_name == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "mesh copy(C->Fortran): confirm allocation: 'struct_name' is NULL");
		return;
	}
	if(var_name == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "mesh copy(C->Fortran): confirm allocation: 'var_name' is NULL");
		return;
	}
	if(is_allocated == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "mesh copy(C->Fortran): confirm allocation: 'is_allocated' is NULL");
		return;
	}

	if(HECMW_strcpy_f2c_r(struct_name, len_struct, sname, sizeof(sname)) == NULL) {
		return;
	}
	if(HECMW_strcpy_f2c_r(var_name, len_var, vname, sizeof(vname)) == NULL) {
		return;
	}

	func = get_is_allocated_func(sname, vname);
	if(func == NULL) {
		HECMW_set_error(HECMW_ALL_E0102, "HECMW_mesh_is_allocated(): IsAllocatedFunc not found");
		return;
	}

	if((*func)()) {
		*is_allocated = 1;
	} else {
		*is_allocated = 0;
	}

	*err = 0;	/* no error */
}


void
hecmw_dist_copy_c2f_isalloc_if_(char *struct_name, char *var_name,
					int *is_allocated, int *err, int len_struct, int len_var)
{
	hecmw_dist_copy_c2f_isalloc_if(struct_name, var_name,
						is_allocated, err, len_struct, len_var);
}


void
hecmw_dist_copy_c2f_isalloc_if__(char *struct_name, char *var_name,
					int *is_allocated, int *err, int len_struct, int len_var)
{
	hecmw_dist_copy_c2f_isalloc_if(struct_name, var_name,
						is_allocated, err, len_struct, len_var);
}


void
HECMW_DIST_COPY_C2F_ISALLOC_IF(char *struct_name, char *var_name,
					int *is_allocated, int *err, int len_struct, int len_var)
{
	hecmw_dist_copy_c2f_isalloc_if(struct_name, var_name,
						is_allocated, err, len_struct, len_var);
}


/*----------------------------------------------------------------------------*/


void
hecmw_dist_copy_c2f_set_if(char *struct_name, char *var_name,
						void *dst, int *err, int len_struct, int len_var)
{
	SetFunc func;
	char sname[HECMW_NAME_LEN+1];
	char vname[HECMW_NAME_LEN+1];

	*err = 1;

	if(mesh == NULL) {
		HECMW_set_error(HECMW_ALL_E0102, "mesh copy(C->Fortran): 'mesh' has not initialized yet");
		return;
	}
	if(struct_name == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "mesh copy(C->Fortran): 'sname' is NULL");
		return;
	}
	if(var_name == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "mesh copy(C->Fotran): 'vname' is NULL");
		return;
	}
	if(dst == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "mesh copy(C->Fortran): 'dst' is NULL");
		return;
	}

	if(HECMW_strcpy_f2c_r(struct_name, len_struct, sname, sizeof(sname)) == NULL) {
		return;
	}
	if(HECMW_strcpy_f2c_r(var_name, len_var, vname, sizeof(vname)) == NULL) {
		return;
	}

	func = get_set_func(sname, vname);
	if(func == NULL) {
		HECMW_set_error(HECMW_ALL_E0102, "mesh copy(C->Fortran): SetFunc not found");
		return;
	}

	if((*func)(dst)) {
		return;
	}

	*err = 0;	/* no error */
}



void
hecmw_dist_copy_c2f_set_if_(char *struct_name, char *var_name,
						void *dst, int *err, int len_struct, int len_var)
{
	hecmw_dist_copy_c2f_set_if(struct_name,var_name,dst,err,len_struct,len_var);
}



void
hecmw_dist_copy_c2f_set_if__(char *struct_name, char *var_name,
						void *dst, int *err, int len_struct, int len_var)
{
	hecmw_dist_copy_c2f_set_if(struct_name,var_name,dst,err,len_struct,len_var);
}



void
HECMW_DIST_COPY_C2F_SET_IF(char *struct_name, char *var_name,
						void *dst, int *err, int len_struct, int len_var)
{
	hecmw_dist_copy_c2f_set_if(struct_name,var_name,dst,err,len_struct,len_var);
}


