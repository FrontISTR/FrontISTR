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
    double node2[] = {
        0.0,1.0,2.0,
        3.0,4.0,5.0,
        6.0,7.0,8.0,
    };
    hecmw_api_mesh_set_node(mesh2,3,node2);
    assert(hecmw_api_mesh_n_node(mesh2)==3);
    int elemtype[] = {
        111,
        231,
        111,
    };
    int elem[] = {
        1,2,
        1,2,3,
        1,3,
    };
    int sectID[] = {
        1,
        2,
        1,
    };
    hecmw_api_mesh_set_element(mesh2,3,elemtype,elem,sectID);
    const int* elemtype2 = hecmw_api_mesh_get_elem_type(mesh2);
    const int* elem2 = hecmw_api_mesh_get_elem_node_item(mesh2);
    const int* sectID2 = hecmw_api_mesh_get_section_id(mesh2);
    assert(elemtype2[0]==111);
    assert(elemtype2[1]==111);
    assert(elemtype2[2]==231);
    assert(elem2[0]==1);
    assert(elem2[1]==2);
    assert(elem2[2]==1);
    assert(elem2[3]==3);
    assert(elem2[4]==1);
    assert(elem2[5]==2);
    assert(elem2[6]==3);
    assert(sectID2[0]==1);
    assert(sectID2[1]==1);
    assert(sectID2[2]==2);
    hecmw_api_mesh_delete(mesh2);
    return 0;
}
