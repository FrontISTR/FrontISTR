/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESULT_IO_INCLUDED
#define HECMW_RESULT_IO_INCLUDED

#include "hecmw_config.h"

enum HECMW_RESULT_DTYPE {
  HECMW_RESULT_DTYPE_MIN    = 1,
  HECMW_RESULT_DTYPE_NODE   = 1,
  HECMW_RESULT_DTYPE_ELEM   = 2,
  HECMW_RESULT_DTYPE_GLOBAL = 3,
  HECMW_RESULT_DTYPE_MAX    = 3
};

#ifdef OLD_RES_FORMAT
# define HECMW_RESULT_FILEVER_MAJOR 1
# define HECMW_RESULT_FILEVER_MINOR 0
#else
# define HECMW_RESULT_FILEVER_MAJOR 2
# define HECMW_RESULT_FILEVER_MINOR 0
#endif // OLD_RES_FORMAT

struct result_list {
  char *label;
  double *ptr;
  int n_dof;
  struct result_list *next;
};

struct hecmwST_result_io_data {
  int istep;
  int nnode;
  int nelem;
  char head[HECMW_HEADER_LEN + 1];
  char comment_line[HECMW_MSG_LEN + 1];

  struct result_list *global_list;
  struct result_list *node_list;
  struct result_list *elem_list;

  int *node_global_ID;
  int *elem_global_ID;

  int MPC_exist;
  int *eid_wo_MPC;
};

extern struct hecmwST_result_io_data ResIO;

extern void HECMW_result_io_finalize();
extern int HECMW_result_io_init(int n_node, int n_elem, int *nodeID, int *elemID,
                                int n_elem_type, int *elem_type_index, int *elem_type_item,
                                int i_step, char *header, char *comment);
extern int HECMW_result_io_add(int dtype, int n_dof, char *label,
                               double *ptr);

extern int HECMW_result_io_count_ng_comp(void);
extern int HECMW_result_io_count_nn_comp(void);
extern int HECMW_result_io_count_ne_comp(void);

#endif
