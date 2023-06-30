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

enum HECMW_RESULT_DTYPE {
  HECMW_RESULT_DTYPE_MIN    = 1,
  HECMW_RESULT_DTYPE_NODE   = 1,
  HECMW_RESULT_DTYPE_ELEM   = 2,
  HECMW_RESULT_DTYPE_GLOBAL = 3,
  HECMW_RESULT_DTYPE_MAX    = 3
};

struct result_list {
  char *label;
  double *ptr;
  int n_dof;
  struct result_list *next;
};

extern int IStep;
extern int NNode;
extern int NElem;
extern int FileVer_Major;
extern int FileVer_Minor;
extern char Head[HECMW_HEADER_LEN + 1];
extern char Comment_Line[HECMW_MSG_LEN + 1];
extern char Line_Buf[LINEBUF_SIZE + 1];

extern struct result_list *Global_List;
extern struct result_list *Node_List;
extern struct result_list *Elem_List;

extern int *Node_Global_ID;
extern int *Elem_Global_ID;


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
