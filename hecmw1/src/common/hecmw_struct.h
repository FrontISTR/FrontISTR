/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_STRUCT_INCLUDED
#define HECMW_STRUCT_INCLUDED

#include "hecmw_util.h"

struct hecmwST_section {
  int n_sect;

  int *sect_type;
#define HECMW_SECT_TYPE_SOLID 1     /* 1:SOLID */
#define HECMW_SECT_TYPE_SHELL 2     /* 2:SHELL */
#define HECMW_SECT_TYPE_BEAM 3      /* 3:BEAM */
#define HECMW_SECT_TYPE_INTERFACE 4 /* 4:INTERFACE */

  int *sect_opt;
#define HECMW_SECT_OPT_PSTRESS 0       /* plane stress */
#define HECMW_SECT_OPT_PSTRAIN 1       /* plane strain */
#define HECMW_SECT_OPT_ASYMMETRY 2     /* axial symmetry */
#define HECMW_SECT_OPT_PSTRESS_RI 10   /* plane stress & reduced integral */
#define HECMW_SECT_OPT_PSTRAIN_RI 11   /* plane strain & reduced integral */
#define HECMW_SECT_OPT_ASYMMETRY_RI 12 /* axial symmetry & reduced integral */
  int *sect_mat_ID_index;
  int *sect_mat_ID_item;
  int *sect_I_index;
  int *sect_I_item;
  int *sect_R_index;
  double *sect_R_item;
};

struct hecmwST_material {
  int n_mat;
  int n_mat_item;
  int n_mat_subitem;
  int n_mat_table;
  char **mat_name;
  int *mat_item_index;
  int *mat_subitem_index;
  int *mat_table_index;
  double *mat_val;
  double *mat_temp;
};

struct hecmwST_mpc {
  int n_mpc;
  int *mpc_index;
  int *mpc_item;
  int *mpc_dof;
  double *mpc_val;
  double *mpc_const;
};

struct hecmwST_amplitude {
  int n_amp;
  char **amp_name;

  int *amp_type_definition;
#define HECMW_AMP_TYPEDEF_TABULAR 1 /* 1:TABULAR(default) */

  int *amp_type_time;
#define HECMW_AMP_TYPETIME_STEP 1 /* 1:STEP_TIME(default) */

  int *amp_type_value;
#define HECMW_AMP_TYPEVAL_RELATIVE 1 /* 1:RELATIVE(default) */
#define HECMW_AMP_TYPEVAL_ABSOLUTE 2 /* 2:ABSOLUTE */
  int *amp_index;
  double *amp_val;
  double *amp_table;
};

struct hecmwST_node_grp {
  int n_grp;
  char **grp_name;
  int *grp_index;
  int *grp_item;

  int n_bc;
  int *bc_grp_ID;

  int *bc_grp_type;
#define HECMW_BCGRPTYPE_DISPALCEMENT 1 /* 1:displacement */
#define HECMW_BCGRPTYPE_FLUX 2         /* 2:flux */
  int *bc_grp_index;
  int *bc_grp_dof;
  double *bc_grp_val;
};

struct hecmwST_elem_grp {
  int n_grp;
  char **grp_name;
  int *grp_index;
  int *grp_item;

  int n_bc;
  int *bc_grp_ID;

  int *bc_grp_type;
  int *bc_grp_index;
  double *bc_grp_val;
};

struct hecmwST_surf_grp {
  int n_grp;
  char **grp_name;
  int *grp_index;

  int *grp_item;
  int n_bc;
  int *bc_grp_ID;

  int *bc_grp_type;

  int *bc_grp_index;
  double *bc_grp_val;
};

struct hecmwST_contact_pair {
  int n_pair;
  char **name;
  int *type;
#define HECMW_CONTACT_TYPE_NODE_SURF 1 /* 1:NODE_SURF */
#define HECMW_CONTACT_TYPE_SURF_SURF 2 /* 2:SURF_SURF */
#define HECMW_CONTACT_TYPE_NODE_ELEM 3 /* 3:NODE_ELEM */
  int *slave_grp_id;
  int *slave_orisgrp_id;
  int *master_grp_id;
};

struct hecmwST_refine_origin {
  int *index;
  int *item_index;
  int *item_item;
};

struct hecmwST_local_mesh {
  int hecmw_flag_adapt;
  int hecmw_flag_initcon;
  int hecmw_flag_parttype;
#define HECMW_FLAG_PARTTYPE_UNKNOWN 0   /* 0:unknown */
#define HECMW_FLAG_PARTTYPE_NODEBASED 1 /* 1:Node-based */
#define HECMW_FLAG_PARTTYPE_ELEMBASED 2 /* 2:Element-based */
  int hecmw_flag_partdepth;
  int hecmw_flag_version;
  int hecmw_flag_partcontact;
#define HECMW_FLAG_PARTCONTACT_UNKNOWN 0    /* 0:unknown */
#define HECMW_FLAG_PARTCONTACT_AGGREGATE 1  /* 1:aggregate */
#define HECMW_FLAG_PARTCONTACT_DISTRIBUTE 2 /* 2:distribute */
#define HECMW_FLAG_PARTCONTACT_SIMPLE 3     /* 3:simple */

  char gridfile[HECMW_FILENAME_LEN + 1];
  int hecmw_n_file;
  char **files;
  char header[HECMW_HEADER_LEN + 1];
  double zero_temp;

  /* Node */
  int n_node;
  int n_node_gross;
  int nn_middle;
  int nn_internal;
  int *node_internal_list;

  int *node_ID;
  int *global_node_ID;

  double *node;
  int n_dof;
  int n_dof_grp;
  int n_dof_tot;
  int *node_dof_index;
  int *node_dof_item;

  int *node_val_index;
  double *node_val_item;

  int *node_init_val_index;
  double *node_init_val_item;

  /* Element */
  int n_elem;
  int n_elem_gross;
  int ne_internal;
  int *elem_internal_list;

  int *elem_ID;
  int *global_elem_ID;
  int *elem_type;
  int n_elem_type;
  int *elem_type_index;
  int *elem_type_item;
  int *elem_node_index;
  int *elem_node_item;
  int *section_ID;
  int *elem_mat_ID_index;
  int *elem_mat_ID_item;
  int n_elem_mat_ID;

  int *elem_mat_int_index;
  double *elem_mat_int_val;
  int *elem_val_index;
  double *elem_val_item;

  /* PE & Communication */
  int zero;
  HECMW_Comm HECMW_COMM;
  int PETOT;
  int PEsmpTOT;
  int my_rank;
  int errnof;
  int n_subdomain;
  int n_neighbor_pe;
  int *neighbor_pe;
  int *import_index;
  int *import_item;
  int *export_index;
  int *export_item;
  int *shared_index;
  int *shared_item;

  /* Adaptation */
  int coarse_grid_level;
  int n_adapt;
  int *when_i_was_refined_node;
  int *when_i_was_refined_elem;
  int *adapt_parent_type;
  int *adapt_type;
  int *adapt_level;

  int *adapt_parent;
  int *adapt_children_index;
  int *adapt_children_item;

  /* Refinement */
  int n_refine;
  int *node_old2new;
  int *node_new2old;
  int *elem_old2new;
  int *elem_new2old;
  int *n_node_refine_hist;

  struct hecmwST_section *section;
  struct hecmwST_material *material;
  struct hecmwST_mpc *mpc;
  struct hecmwST_amplitude *amp;
  struct hecmwST_node_grp *node_group;
  struct hecmwST_elem_grp *elem_group;
  struct hecmwST_surf_grp *surf_group;
  struct hecmwST_contact_pair *contact_pair;
  struct hecmwST_refine_origin *refine_origin;
};

#endif
