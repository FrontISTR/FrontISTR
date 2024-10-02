/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_PART_STRUCT
#define INC_HECMW_PART_STRUCT

#include "hecmw_part_define.h"

struct hecmw_part_edge_data {
  long long int n_edge;

  int *edge_node_item;
};

struct hecmw_part_node_data {
  int *node_elem_index;

  int *node_elem_item;
};

struct hecmw_part_cont_data {
  int n_domain;

  int depth;

  int type;

  int method;

  int n_rcb_div;

  int *rcb_axis;

  int is_print_ucd;

  char ucd_file_name[HECMW_FILENAME_LEN + 1]; /* ucd file name */

  int n_my_domain;

  int *my_domain;

  int contact;

  int is_print_part;

  char part_file_name[HECMW_FILENAME_LEN + 1];
};

#endif /* INC_HECMW_PART_STRUCT */
