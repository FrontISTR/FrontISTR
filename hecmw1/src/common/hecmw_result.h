/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESULT_INCLUDED
#define HECMW_RESULT_INCLUDED

#include "hecmw_struct.h"

struct hecmwST_result_data {
  int ng_component;
  int nn_component;
  int ne_component;
  int *ng_dof;
  int *nn_dof;
  int *ne_dof;
  char **global_label;
  char **node_label;
  char **elem_label;
  double *global_val_item;
  double *node_val_item;
  double *elem_val_item;
};

extern void HECMW_result_free(struct hecmwST_result_data *result);

extern int HECMW_result_init(struct hecmwST_local_mesh *hecMESH,
                             int i_step, char *header, char *comment);
extern int HECMW_result_finalize(void);

extern int HECMW_result_write_by_name(char *name_ID);
extern int HECMW_result_write_by_addfname(char *name_ID, char *addfname);

extern int HECMW_result_write_ST_by_name(char *name_ID,
                                         struct hecmwST_result_data *result,
                                         int n_node, int n_elem, char *header, char *comment);

extern struct hecmwST_result_data *HECMW_result_read_by_name(char *name_ID,
                                                             int i_step);
extern struct hecmwST_result_data *HECMW_result_read_by_fname(char *filename);

extern int HECMW_result_get_nnode(void);
extern int HECMW_result_get_nelem(void);
extern char *HECMW_result_get_header(char *buff);
extern char *HECMW_result_get_comment(char *buff);
extern int *HECMW_result_get_nodeID(int *buff);
extern int *HECMW_result_get_elemID(int *buff);
extern void HECMW_result_free_nodeID(void);
extern void HECMW_result_free_elemID(void);

#endif
