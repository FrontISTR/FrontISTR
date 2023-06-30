/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESULT_IO_INCLUDED
#define HECMW_RESULT_IO_INCLUDED

#include "hecmw_config.h"

#define COL_INT 10
#define COL_DOUBLE 5
#define LINEBUF_SIZE 1023

struct result_list {
  char *label;
  double *ptr;
  int n_dof;
  struct result_list *next;
};

extern int istep;
extern int nnode;
extern int nelem;
extern int filever_major;
extern int filever_minor;
extern char head[HECMW_HEADER_LEN + 1];
extern char comment_line[HECMW_MSG_LEN + 1];
extern char line_buf[LINEBUF_SIZE + 1];

extern struct result_list *global_list;
extern struct result_list *node_list;
extern struct result_list *elem_list;

extern int *node_global_ID;
extern int *elem_global_ID;


extern void HECMW_result_io_finalize();
extern int HECMW_result_io_init(int n_node, int n_elem, int *nodeID,
                                  int *elemID, int i_step,
                                  char *header, char *comment);
extern int HECMW_result_io_add(int dtype, int n_dof, char *label,
                            double *ptr);

extern int HECMW_result_io_count_ng_comp(void);
extern int HECMW_result_io_count_nn_comp(void);
extern int HECMW_result_io_count_ne_comp(void);

#endif
