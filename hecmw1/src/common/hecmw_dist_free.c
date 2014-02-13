/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
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
#include "hecmw_dist_free.h"
#include "hecmw_dist_alloc.h"
#include "hecmw_util.h"
#include "hecmw_struct.h"

static void
free_section(struct hecmwST_section *section)
{
	if(section == NULL) return;
	HECMW_free(section->sect_type);
	HECMW_free(section->sect_opt);
	HECMW_free(section->sect_mat_ID_index);
	HECMW_free(section->sect_mat_ID_item);
	HECMW_free(section->sect_I_index);
	HECMW_free(section->sect_I_item);
	HECMW_free(section->sect_R_index);
	HECMW_free(section->sect_R_item);
}


static void
free_material(struct hecmwST_material *material)
{
	int i;

	if(material == NULL) return;
	for(i=0; i < material->n_mat; i++) {
		HECMW_free(material->mat_name[i]);
	}
	HECMW_free(material->mat_name);
	HECMW_free(material->mat_item_index);
	HECMW_free(material->mat_subitem_index);
	HECMW_free(material->mat_table_index);
	HECMW_free(material->mat_val);
	HECMW_free(material->mat_temp);
}


static void
free_mpc(struct hecmwST_mpc *mpc)
{
	if(mpc == NULL) return;
	HECMW_free(mpc->mpc_index);
	HECMW_free(mpc->mpc_item);
	HECMW_free(mpc->mpc_dof);
	HECMW_free(mpc->mpc_val);
	HECMW_free(mpc->mpc_const);
}

static void
free_amplitude(struct hecmwST_amplitude *amp)
{
	int i;

	if(amp == NULL) return;

	for(i=0; i < amp->n_amp; i++) {
		HECMW_free(amp->amp_name[i]);
	}
	HECMW_free(amp->amp_name);
	HECMW_free(amp->amp_type_definition);
	HECMW_free(amp->amp_type_time);
	HECMW_free(amp->amp_type_value);
	HECMW_free(amp->amp_index);
	HECMW_free(amp->amp_val);
	HECMW_free(amp->amp_table);
}


static void
free_ngrp(struct hecmwST_node_grp *grp)
{
	int i;

	if(grp == NULL) return;

	for(i=0; i < grp->n_grp; i++) {
		HECMW_free(grp->grp_name[i]);
	}
	HECMW_free(grp->grp_name);
	HECMW_free(grp->grp_index);
	HECMW_free(grp->grp_item);
	HECMW_free(grp->bc_grp_ID);
	HECMW_free(grp->bc_grp_type);
	HECMW_free(grp->bc_grp_index);
	HECMW_free(grp->bc_grp_dof);
	HECMW_free(grp->bc_grp_val);
}


static void
free_egrp(struct hecmwST_elem_grp *grp)
{
	int i;

	if(grp == NULL) return;

	for(i=0; i < grp->n_grp; i++) {
		HECMW_free(grp->grp_name[i]);
	}
	HECMW_free(grp->grp_name);
	HECMW_free(grp->grp_index);
	HECMW_free(grp->grp_item);
	HECMW_free(grp->bc_grp_ID);
	HECMW_free(grp->bc_grp_type);
	HECMW_free(grp->bc_grp_index);
	HECMW_free(grp->bc_grp_val);
}


static void
free_sgrp(struct hecmwST_surf_grp *grp)
{
	int i;

	if(grp == NULL) return;

	for(i=0; i < grp->n_grp; i++) {
		HECMW_free(grp->grp_name[i]);
	}
	HECMW_free(grp->grp_name);
	HECMW_free(grp->grp_index);
	HECMW_free(grp->grp_item);
	HECMW_free(grp->bc_grp_ID);
	HECMW_free(grp->bc_grp_type);
	HECMW_free(grp->bc_grp_index);
	HECMW_free(grp->bc_grp_val);
}


static void
free_contact_pair(struct hecmwST_contact_pair *cpair)
{
	int i;

	if(cpair == NULL) return;

	for(i=0; i < cpair->n_pair; i++) {
		HECMW_free(cpair->name[i]);
	}

	HECMW_free(cpair->type);
	HECMW_free(cpair->slave_grp_id);
	HECMW_free(cpair->master_grp_id);
}


static void
free_refine_origin(struct hecmwST_refine_origin *reforg)
{
	int i;

	if(reforg == NULL) return;

	HECMW_free(reforg->index);
	HECMW_free(reforg->item_index);
	HECMW_free(reforg->item_item);
}


void
HECMW_dist_clean(struct hecmwST_local_mesh *mesh)
{
	int i;

	if(mesh == NULL) return;

	for(i=0; i < mesh->hecmw_n_file; i++) {
		HECMW_free(mesh->files[i]);
	}
	HECMW_free(mesh->files);

	HECMW_free(mesh->node_internal_list);
	HECMW_free(mesh->node_ID);
	HECMW_free(mesh->global_node_ID);
	HECMW_free(mesh->node);
	HECMW_free(mesh->node_dof_index);
	HECMW_free(mesh->node_dof_item);
	HECMW_free(mesh->node_val_index);
	HECMW_free(mesh->node_val_item);
	HECMW_free(mesh->node_init_val_index);
	HECMW_free(mesh->node_init_val_item);

	HECMW_free(mesh->elem_internal_list);
	HECMW_free(mesh->elem_ID);
	HECMW_free(mesh->global_elem_ID);
	HECMW_free(mesh->elem_type);
	HECMW_free(mesh->elem_type_index);
	HECMW_free(mesh->elem_type_item);
	HECMW_free(mesh->elem_node_index);
	HECMW_free(mesh->elem_node_item);
	HECMW_free(mesh->section_ID);
	HECMW_free(mesh->elem_mat_ID_index);
	HECMW_free(mesh->elem_mat_ID_item);
	HECMW_free(mesh->elem_mat_int_index);
	HECMW_free(mesh->elem_mat_int_val);
	HECMW_free(mesh->elem_val_index);
	HECMW_free(mesh->elem_val_item);

	HECMW_free(mesh->neighbor_pe);
	HECMW_free(mesh->import_index);
	HECMW_free(mesh->import_item);
	HECMW_free(mesh->export_index);
	HECMW_free(mesh->export_item);
	HECMW_free(mesh->shared_index);
	HECMW_free(mesh->shared_item);

	HECMW_free(mesh->when_i_was_refined_node);
	HECMW_free(mesh->when_i_was_refined_elem);
	HECMW_free(mesh->adapt_parent_type);
	HECMW_free(mesh->adapt_type);
	HECMW_free(mesh->adapt_level);
	HECMW_free(mesh->adapt_parent);
	HECMW_free(mesh->adapt_children_index);
	HECMW_free(mesh->adapt_children_item);

	HECMW_free(mesh->node_old2new);
	HECMW_free(mesh->node_new2old);
	HECMW_free(mesh->elem_old2new);
	HECMW_free(mesh->elem_new2old);

	free_section(mesh->section);
	free_material(mesh->material);
	free_mpc(mesh->mpc);
	free_amplitude(mesh->amp);
	free_ngrp(mesh->node_group);
	free_egrp(mesh->elem_group);
	free_sgrp(mesh->surf_group);
	free_contact_pair(mesh->contact_pair);
	free_refine_origin(mesh->refine_origin);

	HECMW_dist_init(mesh);
}


void
HECMW_dist_free(struct hecmwST_local_mesh *mesh)
{
	if(mesh == NULL) return;

	HECMW_dist_clean(mesh);

	HECMW_free(mesh->section);
	HECMW_free(mesh->material);
	HECMW_free(mesh->mpc);
	HECMW_free(mesh->amp);
	HECMW_free(mesh->node_group);
	HECMW_free(mesh->elem_group);
	HECMW_free(mesh->surf_group);
	HECMW_free(mesh->contact_pair);
	HECMW_free(mesh->refine_origin);

	HECMW_free(mesh);
}

