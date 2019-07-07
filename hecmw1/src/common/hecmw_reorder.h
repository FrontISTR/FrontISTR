/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_REORDER
#define INC_HECMW_REORDER

#include "hecmw_struct.h"

extern int HECMW_reorder_node_mpc(struct hecmwST_local_mesh *local_mesh);
extern int HECMW_reorder_node_dof(struct hecmwST_local_mesh *local_mesh);
extern int HECMW_reorder_elem_type(struct hecmwST_local_mesh *local_mesh);
extern int HECMW_reorder(struct hecmwST_local_mesh *local_mesh);

#endif
