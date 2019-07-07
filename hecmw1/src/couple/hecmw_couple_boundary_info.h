/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_BOUNDARY_INFO
#define INC_HECMW_COUPLE_BOUNDARY_INFO

#include "hecmw_struct.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_control.h"

struct hecmw_couple_boundary_item {
  int n;
  int *item;
};

struct hecmw_couple_boundary {
  int geom_type;
  int data_type;
  struct hecmw_couple_boundary_item *node;
  struct hecmw_couple_boundary_item *elem;
  struct hecmw_couple_boundary_item *surf;
  int *elem_node_index;
  int *elem_node_item;
};

extern void HECMW_couple_free_boundary_info(
    struct hecmw_couple_boundary *boundary);

extern struct hecmw_couple_boundary *HECMW_couple_alloc_boundary_info(void);

extern struct hecmw_couple_boundary *HECMW_couple_set_boundary_info(
    const char *boundary_id, int unit_specifier,
    const struct hecmwST_local_mesh *mesh);

#endif /* INC_HECMW_COUPLE_BOUNDARY_INFO */
