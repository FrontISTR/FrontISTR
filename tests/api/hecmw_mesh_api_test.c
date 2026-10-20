/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#undef NDEBUG
#include <stdio.h>
#include <assert.h>
#include "api/hecmw_api_local_mesh.h"

int main()
{
    printf("FrontISTR API Test : hecmwST_local_mesh\n");
    void* mesh1 = hecmw_api_mesh_new();
    double node1[] = {0.0,1.0,2.0};
    hecmw_api_mesh_set_node(mesh1,1,node1);
    assert(hecmw_api_mesh_n_node(mesh1)==1);
    hecmw_api_mesh_delete(mesh1);

    void* mesh2 = hecmw_api_mesh_new();
    double node2[] = {0.0,1.0,2.0,3.0,4.0,5.0};
    hecmw_api_mesh_set_node(mesh2,2,node2);
    assert(hecmw_api_mesh_n_node(mesh2)==2);
    hecmw_api_mesh_delete(mesh2);
    return 0;
}