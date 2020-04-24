/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_MESH_EDGE_INFO
#define INC_HECMW_MESH_EDGE_INFO

#include "hecmw_struct.h"
#include "hecmw_part_struct.h"

extern int HECMW_mesh_edge_info(struct hecmwST_local_mesh *mesh,
                                struct hecmw_part_edge_data *edge_data,
                                const int edge_create_type);

#endif /* INC_HECMW_MESH_EDGE_INFO */
