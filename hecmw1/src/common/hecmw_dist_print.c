/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                       Shin'ichi Ezure (RIST)                        *
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
#include "hecmw_dist_print.h"


void
HECMW_dist_print_flags(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "FLAGS:\n");
	fprintf(fp, "hecmw_flag_adapt: %d\n", mesh->hecmw_flag_adapt);
	fprintf(fp, "hecmw_flag_initcon: %d\n", mesh->hecmw_flag_initcon);
	fprintf(fp, "hecmw_flag_parttype: %d\n", mesh->hecmw_flag_parttype);
	fprintf(fp, "hecmw_flag_version: %d\n", mesh->hecmw_flag_version);
	fprintf(fp, "END of FLAGS\n");
}


void
HECMW_dist_print_header(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "HEADER:\n");
	fprintf(fp, "%s\n", mesh->header);
	fprintf(fp, "END of HEADER\n");
}


void
HECMW_dist_print_gridfile(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "GRIDFILE:\n");
	fprintf(fp, "%s\n", mesh->gridfile);
	fprintf(fp, "END of GRIDFILE\n");
}


void
HECMW_dist_print_files(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	int i;

	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "FILES:\n");
	fprintf(fp, "hecmw_n_file: %d\n", mesh->hecmw_n_file);
	for(i=0; mesh->files && i < mesh->hecmw_n_file; i++) {
		fprintf(fp, "%s\n", mesh->files[i]);
	}
	fprintf(fp, "END of FILES\n");
}


void
HECMW_dist_print_zero_temp(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "ZERO:\n");
	fprintf(fp, "%E\n", mesh->zero_temp);
	fprintf(fp, "END of ZERO\n");
}


void
HECMW_dist_print_node(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "NODE:\n");
	fprintf(fp, "n_node: %d\n", mesh->n_node);
	fprintf(fp, "nn_internal: %d\n", mesh->nn_internal);

	fprintf(fp, "node_internal_list:%s\n",
			mesh->node_internal_list ? "" : " NULL");
	for(i=0; mesh->node_internal_list && i < mesh->nn_internal; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->node_internal_list[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node_ID:%s\n", mesh->node_ID ? "" : " NULL");
	for(i=0; mesh->node_ID && i < mesh->n_node; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d %d  ", mesh->node_ID[i*2], mesh->node_ID[i*2+1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "global_node_ID:%s\n", mesh->global_node_ID ? "" : " NULL");
	for(i=0; mesh->global_node_ID && i < mesh->n_node; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->global_node_ID[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node:%s\n", mesh->node ? "" : " NULL");
	for(i=0; mesh->node && i < mesh->n_node; i++) {
		fprintf(fp, "%E %E %E\n",
				mesh->node[i*3], mesh->node[i*3+1], mesh->node[i*3+2]);
	}

	fprintf(fp, "n_dof: %d\n", mesh->n_dof);
	fprintf(fp, "n_dof_grp: %d\n", mesh->n_dof_grp);
	fprintf(fp, "n_dof_tot: %d\n", mesh->n_dof_tot);

	fprintf(fp, "node_dof_index:%s\n", mesh->node_dof_index ? "" : " NULL");
	for(i=0; mesh->node_dof_index && i <= mesh->n_dof_grp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->node_dof_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node_dof_item:%s\n", mesh->node_dof_item ? "" : " NULL");
	for(i=0; mesh->node_dof_index && mesh->node_dof_item &&
			i < mesh->n_dof_grp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->node_dof_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node_val_index:%s\n", mesh->node_val_index ? "" : " NULL");
	for(i=0; mesh->node_val_index && i <= mesh->n_node; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->node_val_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node_val_item:%s\n", mesh->node_val_item ? "" : " NULL");
	for(i=0; mesh->node_val_index && mesh->node_val_item &&
			i < mesh->node_val_index[mesh->n_node]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", mesh->node_val_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node_init_val_index:%s\n",
			mesh->node_init_val_index ? "" : " NULL");
	for(i=0; mesh->node_init_val_index && i <= mesh->n_node; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->node_init_val_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "node_init_val_item:%s\n",
			mesh->node_init_val_item ? "" : " NULL");
	for(i=0; mesh->node_init_val_index && mesh->node_init_val_item &&
			i < mesh->node_init_val_index[mesh->n_node]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", mesh->node_init_val_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of NODE\n");
}


void
HECMW_dist_print_elem(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "ELEMENT:\n");

	fprintf(fp, "n_elem: %d\n", mesh->n_elem);
	fprintf(fp, "ne_internal: %d\n", mesh->ne_internal);

	fprintf(fp, "elem_internal_list:%s\n",
			mesh->elem_internal_list ? "" : " NULL");
	for(i=0; mesh->elem_internal_list && i < mesh->ne_internal; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_internal_list[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_ID:%s\n", mesh->elem_ID ? "" : " NULL");
	for(i=0; mesh->elem_ID && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d %d  ", mesh->elem_ID[2*i], mesh->elem_ID[2*i+1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "global_elem_ID:%s\n", mesh->global_elem_ID ? "" : " NULL");
	for(i=0; mesh->global_elem_ID && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->global_elem_ID[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_type:%s\n", mesh->elem_type ? "" : " NULL");
	for(i=0; mesh->elem_type && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_elem_type: %d\n", mesh->n_elem_type);

	fprintf(fp, "elem_type_index:%s\n", mesh->elem_type_index ? "" : " NULL");
	for(i=0; mesh->elem_type_index && i <= mesh->n_elem_type; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_type_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_type_item:%s\n", mesh->elem_type_item ? "" : " NULL");
	for(i=0; mesh->elem_type_index && mesh->elem_type_item &&
			i < mesh->n_elem_type; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_type_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_node_index:%s\n", mesh->elem_node_index ? "" : " NULL");
	for(i=0; mesh->elem_node_index && i <= mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_node_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_node_item:%s\n", mesh->elem_node_item ? "" : " NULL");
	for(i=0; mesh->elem_node_index && mesh->elem_node_item &&
			i < mesh->elem_node_index[mesh->n_elem]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_node_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "section_ID:%s\n", mesh->section_ID ? "" : " NULL");
	for(i=0; mesh->section_ID && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->section_ID[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_elem_mat_ID:%d\n", mesh->n_elem_mat_ID);

	fprintf(fp, "elem_mat_ID_index:%s\n",
			mesh->elem_mat_ID_index ? "" : " NULL");
	for(i=0; mesh->elem_mat_ID_index && i <= mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_mat_ID_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_mat_ID_item:%s\n", mesh->elem_mat_ID_item ? "" : " NULL");
	for(i=0; mesh->elem_mat_ID_index && mesh->elem_mat_ID_item &&
			i < mesh->elem_mat_ID_index[mesh->n_elem]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_mat_ID_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_mat_int_index:%s\n",
			mesh->elem_mat_int_index ? "" : " NULL");
	for(i=0; mesh->elem_mat_int_index && i <= mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_mat_int_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_mat_int_val:%s\n", mesh->elem_mat_int_val ? "" : " NULL");
	for(i=0; mesh->elem_mat_int_index && mesh->elem_mat_int_val &&
			i < mesh->elem_mat_int_index[mesh->n_elem]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", mesh->elem_mat_int_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_val_index:%s\n", mesh->elem_val_index ? "" : " NULL");
	for(i=0; mesh->elem_val_index && i <= mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->elem_val_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "elem_val_item:%s\n", mesh->elem_val_item ? "" : " NULL");
	for(i=0; mesh->elem_val_index && mesh->elem_val_item &&
			i < mesh->elem_val_index[mesh->n_elem]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", mesh->elem_val_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of ELEMENT\n");
}


void
HECMW_dist_print_pe(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "PE:\n");

	fprintf(fp, "zero: %d\n", mesh->zero);
	fprintf(fp, "HECMW_COMM: %ld\n", (long)mesh->HECMW_COMM);
	fprintf(fp, "PETOT: %d\n", mesh->PETOT);
	fprintf(fp, "PEsmpTOT: %d\n", mesh->PEsmpTOT);
	fprintf(fp, "my_rank: %d\n", mesh->my_rank);
	fprintf(fp, "errnof: %d\n", mesh->errnof);
	fprintf(fp, "n_subdomain: %d\n", mesh->n_subdomain);
	fprintf(fp, "n_neighbor_pe: %d\n", mesh->n_neighbor_pe);

	fprintf(fp, "neighbor_pe:%s\n", mesh->neighbor_pe ? "" : " NULL");
	for(i=0; mesh->neighbor_pe && i < mesh->n_neighbor_pe; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->neighbor_pe[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "import_index:%s\n", mesh->import_index ? "" : " NULL");
	for(i=0; mesh->import_index && i <= mesh->n_neighbor_pe; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->import_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "import_item:%s\n", mesh->import_item ? "" : " NULL");
	for(i=0; mesh->import_index && mesh->import_item &&
			i < mesh->import_index[mesh->n_neighbor_pe]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->import_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "export_index:%s\n", mesh->export_index ? "" : " NULL");
	for(i=0; mesh->export_index && i <= mesh->n_neighbor_pe; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->export_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "export_item:%s\n", mesh->export_item ? "" : " NULL");
	for(i=0; mesh->export_index && mesh->export_item &&
			i < mesh->export_index[mesh->n_neighbor_pe]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->export_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "shared_index:%s\n", mesh->shared_index ? "" : " NULL");
	for(i=0; mesh->shared_index && i <= mesh->n_neighbor_pe; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->shared_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "shared_item:%s\n", mesh->shared_item ? "" : " NULL");
	for(i=0; mesh->shared_index && mesh->shared_item &&
			i < mesh->shared_index[mesh->n_neighbor_pe]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->shared_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of PE\n");
}


void
HECMW_dist_print_adapt(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(mesh == NULL) return;
	if(fp == NULL) return;

	fprintf(fp, "ADAPTATION:\n");

	fprintf(fp, "coarse_grid_level: %d\n", mesh->coarse_grid_level);
	fprintf(fp, "n_adapt: %d\n", mesh->n_adapt);

	fprintf(fp, "when_i_was_refined_node:%s\n",
			mesh->when_i_was_refined_node ? "" : " NULL");
	for(i=0; mesh->when_i_was_refined_node && i < mesh->n_node; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->when_i_was_refined_node[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "when_i_was_refined_elem:%s\n",
			mesh->when_i_was_refined_elem ? "" : " NULL");
	for(i=0; mesh->when_i_was_refined_elem && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->when_i_was_refined_elem[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "adapt_parent_type:%s\n",
			mesh->adapt_parent_type ? "" : " NULL");
	for(i=0; mesh->adapt_parent_type && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->adapt_parent_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "adapt_type:%s\n",
			mesh->adapt_type ? "" : " NULL");
	for(i=0; mesh->adapt_type && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->adapt_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "adapt_level:%s\n",
			mesh->adapt_level ? "" : " NULL");
	for(i=0; mesh->adapt_level && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->adapt_level[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "adapt_parent:%s\n",
			mesh->adapt_parent ? "" : " NULL");
	for(i=0; mesh->adapt_parent && i < mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d %d  ",
				mesh->adapt_parent[2*i], mesh->adapt_parent[2*i+1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "adapt_children_index:%s\n",
			mesh->adapt_children_index ? "" : " NULL");
	for(i=0; mesh->adapt_children_index && i <= mesh->n_elem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mesh->adapt_children_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "adapt_children_item:%s\n",
			mesh->adapt_children_index ? "" : " NULL");
	for(i=0; mesh->adapt_children_index && mesh->adapt_children_item &&
			i < mesh->adapt_children_index[mesh->n_elem]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d %d  ",
				mesh->adapt_children_item[2*i],
				mesh->adapt_children_item[2*i+1]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of ADAPTATION:\n");
}


void
HECMW_dist_print_section(const struct hecmwST_section *sect, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "SECTION:%s\n", sect ? "" : " NULL");
	if(sect == NULL) return;

	fprintf(fp, "n_sect: %d\n", sect->n_sect);

	fprintf(fp, "sect_type:%s\n", sect->sect_type ? "" : " NULL");
	for(i=0; sect->sect_type && i < sect->n_sect; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_opt:%s\n", sect->sect_opt ? "" : " NULL");
	for(i=0;  sect->sect_opt && i < sect->n_sect; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_opt[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_mat_ID_index:%s\n",
			sect->sect_mat_ID_index ? "" : " NULL");
	for(i=0;  sect->sect_mat_ID_index && i <= sect->n_sect; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_mat_ID_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_mat_ID_item:%s\n", sect->sect_mat_ID_item ? "" : " NULL");
	for(i=0;  sect->sect_mat_ID_index && sect->sect_mat_ID_item &&
			i < sect->sect_mat_ID_index[sect->n_sect]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_mat_ID_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_I_index:%s\n", sect->sect_I_index ? "" : " NULL");
	for(i=0;  sect->sect_I_index && i <= sect->n_sect; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_I_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_I_item:%s\n", sect->sect_I_item ? "" : " NULL");
	for(i=0;  sect->sect_I_index && sect->sect_I_item &&
			i < sect->sect_I_index[sect->n_sect]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_I_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_R_index:%s\n", sect->sect_R_index ? "" : " NULL");
	for(i=0;  sect->sect_R_index && i <= sect->n_sect; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sect->sect_R_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "sect_R_item:%s\n", sect->sect_R_item ? "" : " NULL");
	for(i=0;  sect->sect_R_index && sect->sect_R_item &&
			i < sect->sect_R_index[sect->n_sect]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", sect->sect_R_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of SECTION\n");
}


void
HECMW_dist_print_material(const struct hecmwST_material *material, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "MATERIAL:%s\n", material ? "" : " NULL");
	if(material == NULL) return;

	fprintf(fp, "n_mat: %d\n", material->n_mat);

	fprintf(fp, "mat_name:%s\n", material->mat_name ? "" : " NULL");
	for(i=0; material->mat_name && i < material->n_mat; i++) {
		fprintf(fp, "%s\n", material->mat_name[i]);
	}

	fprintf(fp, "n_mat_item: %d\n", material->n_mat_item);

	fprintf(fp, "mat_item_index:%s\n", material->mat_item_index ? "" : " NULL");
	for(i=0; material->mat_item_index && i <= material->n_mat; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", material->mat_item_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_mat_subitem: %d\n", material->n_mat_subitem);

	fprintf(fp, "mat_subitem_index:%s\n",
			material->mat_subitem_index ? "" : " NULL");
	for(i=0; material->mat_subitem_index && i <= material->n_mat_item; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", material->mat_subitem_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_mat_table: %d\n", material->n_mat_table);

	fprintf(fp, "mat_table_index:%s\n",
			material->mat_table_index ? "" : " NULL");
	for(i=0; material->mat_table_index && i <= material->n_mat_subitem; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", material->mat_table_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "mat_val:%s\n", material->mat_val ? "" : " NULL");
	for(i=0; material->mat_table_index && material->mat_val &&
			i < material->mat_table_index[material->n_mat_subitem]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", material->mat_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "mat_temp:%s\n", material->mat_temp ? "" : " NULL");
	for(i=0; material->mat_temp && i < material->n_mat_table; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", material->mat_temp[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of MATERIAL\n");
}


void
HECMW_dist_print_mpc(const struct hecmwST_mpc *mpc, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "MPC:%s\n", mpc ? "" : " NULL");
	if(mpc == NULL) return;

	fprintf(fp, "n_mpc: %d\n", mpc->n_mpc);

	fprintf(fp, "mpc_index:%s\n", mpc->mpc_index ? "" : " NULL");
	for(i=0; mpc->mpc_index && i <= mpc->n_mpc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mpc->mpc_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "mpc_item:%s\n", mpc->mpc_item ? "" : " NULL");
	for(i=0; mpc->mpc_index && mpc->mpc_item &&
			i < mpc->mpc_index[mpc->n_mpc]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mpc->mpc_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "mpc_dof:%s\n", mpc->mpc_dof ? "" : " NULL");
	for(i=0; mpc->mpc_index && mpc->mpc_dof &&
			i < mpc->mpc_index[mpc->n_mpc]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", mpc->mpc_dof[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "mpc_val:%s\n", mpc->mpc_val ? "" : " NULL");
	for(i=0; mpc->mpc_index && mpc->mpc_val &&
			i < mpc->mpc_index[mpc->n_mpc]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", mpc->mpc_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of MPC\n");
}


void
HECMW_dist_print_amp(const struct hecmwST_amplitude *amp, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "AMPLITUDE:%s\n", amp ? "" : " NULL");
	if(amp == NULL) return;

	fprintf(fp, "n_amp: %d\n", amp->n_amp);

	fprintf(fp, "amp_name:%s\n", amp->amp_name ? "" : "NULL");
	for(i=0; amp->amp_name && i < amp->n_amp; i++) {
		fprintf(fp, "%s\n", amp->amp_name[i]);
	}

	fprintf(fp, "amp_type_definition:%s\n",
			amp->amp_type_definition ? "" : " NULL");
	for(i=0; amp->amp_type_definition && i < amp->n_amp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", amp->amp_type_definition[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "amp_type_time:%s\n", amp->amp_type_time ? "" : " NULL");
	for(i=0; amp->amp_type_time && i < amp->n_amp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", amp->amp_type_time[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "amp_type_value:%s\n", amp->amp_type_value ? "" : " NULL");
	for(i=0; amp->amp_type_value && i < amp->n_amp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", amp->amp_type_value[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "amp_index:%s\n", amp->amp_index ? "" : " NULL");
	for(i=0; amp->amp_index && i <= amp->n_amp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", amp->amp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "amp_val:%s\n", amp->amp_val ? "" : " NULL");
	for(i=0; amp->amp_index && amp->amp_val &&
			i < amp->amp_index[amp->n_amp]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", amp->amp_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "amp_table:%s\n", amp->amp_table ? "" : " NULL");
	for(i=0; amp->amp_index && amp->amp_table &&
			i < amp->amp_index[amp->n_amp]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", amp->amp_table[i]);
	}
	fprintf(fp, "\n");
}


void
HECMW_dist_print_ngrp(const struct hecmwST_node_grp *ngrp, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "NGROUP:%s\n", ngrp ? "" : " NULL");
	if(ngrp == NULL) return;

	fprintf(fp, "n_grp: %d\n", ngrp->n_grp);

	fprintf(fp, "grp_name:%s\n", ngrp->grp_name ? "" : "NULL");
	for(i=0; i < ngrp->n_grp; i++) {
		fprintf(fp, "%s\n", ngrp->grp_name[i]);
	}

	fprintf(fp, "grp_index:%s\n", ngrp->grp_index ? "" : " NULL");
	for(i=0; ngrp->grp_index && i <= ngrp->n_grp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", ngrp->grp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "grp_item:%s\n", ngrp->grp_item ? "" : " NULL");
	for(i=0; ngrp->grp_index && ngrp->grp_item &&
			i < ngrp->grp_index[ngrp->n_grp]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", ngrp->grp_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_bc: %d\n", ngrp->n_bc);

	fprintf(fp, "bc_grp_ID:%s\n", ngrp->bc_grp_ID ? "" : " NULL");
	for(i=0; ngrp->bc_grp_ID && i < ngrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", ngrp->bc_grp_ID[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_type:%s\n", ngrp->bc_grp_type ? "" : " NULL");
	for(i=0; ngrp->bc_grp_type && i < ngrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", ngrp->bc_grp_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_index:%s\n", ngrp->bc_grp_index ? "" : " NULL");
	for(i=0; ngrp->bc_grp_index && i <= ngrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", ngrp->bc_grp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_dof:%s\n", ngrp->bc_grp_dof ? "" : " NULL");
	for(i=0; ngrp->bc_grp_index && ngrp->bc_grp_dof &&
			i < ngrp->bc_grp_index[ngrp->n_bc]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", ngrp->bc_grp_dof[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_val:%s\n", ngrp->bc_grp_val ? "" : " NULL");
	for(i=0; ngrp->bc_grp_index && ngrp->bc_grp_val &&
			i < ngrp->bc_grp_index[ngrp->n_bc]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", ngrp->bc_grp_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of NGROUP\n");
}


void
HECMW_dist_print_egrp(const struct hecmwST_elem_grp *egrp, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "EGROUP:%s\n", egrp ? "" : " NULL");
	if(egrp == NULL) return;

	fprintf(fp, "n_grp: %d\n", egrp->n_grp);

	fprintf(fp, "grp_name:%s\n", egrp->grp_name ? "" : "NULL");
	for(i=0; i < egrp->n_grp; i++) {
		fprintf(fp, "%s\n", egrp->grp_name[i]);
	}

	fprintf(fp, "grp_index:%s\n", egrp->grp_index ? "" : " NULL");
	for(i=0; egrp->grp_index && i <= egrp->n_grp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", egrp->grp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "grp_item:%s\n", egrp->grp_item ? "" : " NULL");
	for(i=0; egrp->grp_index && egrp->grp_item &&
			i < egrp->grp_index[egrp->n_grp]; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", egrp->grp_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_bc: %d\n", egrp->n_bc);

	fprintf(fp, "bc_grp_ID:%s\n", egrp->bc_grp_ID ? "" : " NULL");
	for(i=0; egrp->bc_grp_ID && i < egrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", egrp->bc_grp_ID[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_type:%s\n", egrp->bc_grp_type ? "" : " NULL");
	for(i=0; egrp->bc_grp_type && i < egrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", egrp->bc_grp_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_index:%s\n", egrp->bc_grp_index ? "" : " NULL");
	for(i=0; egrp->bc_grp_index && i <= egrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", egrp->bc_grp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_val:%s\n", egrp->bc_grp_val ? "" : " NULL");
	for(i=0; egrp->bc_grp_index && egrp->bc_grp_val &&
			i < egrp->bc_grp_index[egrp->n_bc]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", egrp->bc_grp_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of EGROUP\n");
}


void
HECMW_dist_print_sgrp(const struct hecmwST_surf_grp *sgrp, FILE *fp)
{
	int i;
	const int NITEM = 10;

	if(fp == NULL) return;

	fprintf(fp, "SGROUP:%s\n", sgrp ? "" : " NULL");
	if(sgrp == NULL) return;

	fprintf(fp, "n_grp: %d\n", sgrp->n_grp);

	fprintf(fp, "grp_name:%s\n", sgrp->grp_name ? "" : "NULL");
	for(i=0; i < sgrp->n_grp; i++) {
		fprintf(fp, "%s\n", sgrp->grp_name[i]);
	}

	fprintf(fp, "grp_index:%s\n", sgrp->grp_index ? "" : " NULL");
	for(i=0; sgrp->grp_index && i <= sgrp->n_grp; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sgrp->grp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "grp_item:%s\n", sgrp->grp_item ? "" : " NULL");
	for(i=0; sgrp->grp_index && sgrp->grp_item &&
			i < sgrp->grp_index[sgrp->n_grp] * 2; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sgrp->grp_item[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "n_bc: %d\n", sgrp->n_bc);

	fprintf(fp, "bc_grp_ID:%s\n", sgrp->bc_grp_ID ? "" : " NULL");
	for(i=0; sgrp->bc_grp_ID && i < sgrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sgrp->bc_grp_ID[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_type:%s\n", sgrp->bc_grp_type ? "" : " NULL");
	for(i=0; sgrp->bc_grp_type && i < sgrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sgrp->bc_grp_type[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_index:%s\n", sgrp->bc_grp_index ? "" : " NULL");
	for(i=0; sgrp->bc_grp_index && i <= sgrp->n_bc; i++) {
		if(i != 0 && i % NITEM == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%d ", sgrp->bc_grp_index[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "bc_grp_val:%s\n", sgrp->bc_grp_val ? "" : " NULL");
	for(i=0; sgrp->bc_grp_index && sgrp->bc_grp_val &&
			i < sgrp->bc_grp_index[sgrp->n_bc]; i++) {
		if(i != 0 && i % (NITEM/2) == 0) {
			fprintf(fp, "\n");
		}
		fprintf(fp, "%E ", sgrp->bc_grp_val[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "END of SGROUP\n");
}


void
HECMW_dist_print(const struct hecmwST_local_mesh *mesh, FILE *fp)
{
	HECMW_dist_print_flags(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_header(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_gridfile(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_files(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_zero_temp(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_node(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_elem(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_pe(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_adapt(mesh, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_section(mesh->section, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_material(mesh->material, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_mpc(mesh->mpc, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_amp(mesh->amp, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_ngrp(mesh->node_group, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_egrp(mesh->elem_group, fp);
	fprintf(fp, "\n");
	HECMW_dist_print_sgrp(mesh->surf_group, fp);
	fprintf(fp, "\n");
}
