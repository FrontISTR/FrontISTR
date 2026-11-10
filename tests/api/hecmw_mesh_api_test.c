/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#undef NDEBUG
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "api/hecmw_api_local_mesh.h"

int mesh_test(){
    printf("FrontISTR API Test : hecmwST_local_mesh\n");
    void* mesh1 = hecmw_api_mesh_new();
    void* mesh2 = hecmw_api_mesh_new();

    assert(hecmw_api_mesh_n_node(mesh1)==0);
    assert(hecmw_api_mesh_n_node(mesh2)==0);
    assert(hecmw_api_mesh_n_elem(mesh1)==0);
    assert(hecmw_api_mesh_n_elem(mesh2)==0);
    assert(hecmw_api_mesh_n_elem_node_item(mesh1)==0);
    assert(hecmw_api_mesh_n_elem_node_item(mesh2)==0);
    assert(hecmw_api_mesh_n_ngrp(mesh1)==0);
    assert(hecmw_api_mesh_n_ngrp(mesh2)==0);
    assert(hecmw_api_mesh_n_sgrp(mesh1)==0);
    assert(hecmw_api_mesh_n_sgrp(mesh2)==0);
    assert(hecmw_api_mesh_n_egrp(mesh1)==0);
    assert(hecmw_api_mesh_n_egrp(mesh2)==0);

    double node1[] = {
        0.0,1.0,2.0
    };
    hecmw_api_mesh_set_node(mesh1,1,node1);

    double node2[] = {
        0.0,1.0,2.0,
        3.0,4.0,5.0,
        6.0,7.0,8.0,
    };
    hecmw_api_mesh_set_node(mesh2,3,node2);

    assert(hecmw_api_mesh_n_node(mesh1)==1);
    assert(hecmw_api_mesh_n_node(mesh2)==3);
    const double* node = hecmw_api_mesh_get_node(mesh2);
    assert(node[0]==0.0);
    assert(node[3]==3.0);
    assert(node[6]==6.0);

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
    hecmw_api_mesh_delete(mesh1);
    hecmw_api_mesh_delete(mesh2);
    return 0;
}

int surface_group_test()
{
    printf("FrontISTR API Test : surface group\n");
    void* mesh = hecmw_api_mesh_new();   
    double node[] = {
        0.0,0.0,0.0,
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        0.0,0.0,1.0,
    };
    hecmw_api_mesh_set_node(mesh,4,node);
    int elemtype[] = {
        341,
    };
    int elem[] = {
        1,2,3,4,
    };
    int sectID[] = {
        1,
    };
    hecmw_api_mesh_set_element(mesh,1,elemtype,elem,sectID);

    assert(hecmw_api_mesh_n_node(mesh)==4);
    assert(hecmw_api_mesh_n_elem(mesh)==1);
    assert(hecmw_api_mesh_n_elem_node_item(mesh)==4);
    assert(hecmw_api_mesh_n_ngrp(mesh)==0);
    assert(hecmw_api_mesh_n_sgrp(mesh)==0);
    assert(hecmw_api_mesh_n_egrp(mesh)==0);
    const int* elemtype2 = hecmw_api_mesh_get_elem_type(mesh);
    const int* elem2 = hecmw_api_mesh_get_elem_node_item(mesh);
    const int* sectID2 = hecmw_api_mesh_get_section_id(mesh);
    assert(elemtype2[0]==341);
    assert(elem2[0]==1);
    assert(elem2[1]==2);
    assert(elem2[2]==3);
    assert(elem2[3]==4);
    assert(sectID2[0]==1);

    int sgrp[] = {
        1,1,1,2
    };
    hecmw_api_mesh_append_sgrp(mesh,"S1",2,sgrp);
    assert(hecmw_api_mesh_n_sgrp(mesh)==1);
    char buf[64];
    hecmw_api_mesh_get_sgrp_name(mesh,1,buf,64);
    assert(strcmp("S1",buf)==0);
    int* array;
    int count=0;
    hecmw_api_mesh_get_sgrp(mesh,1,&array,&count);
    assert(count==2);
    assert(array[0]==1);
    assert(array[1]==1);
    assert(array[2]==1);
    assert(array[3]==2);
    hecmw_api_mesh_delete(mesh);
    return 0;
}

int node_group_test()
{
    printf("FrontISTR API Test : node group\n");
    void* mesh = hecmw_api_mesh_new();   
    double node[] = {
        0.0,0.0,0.0,
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        1.0,1.0,0.0,
        0.0,0.0,1.0,
        1.0,0.0,1.0,
        0.0,1.0,1.0,
        1.0,1.0,1.0,
    };
    hecmw_api_mesh_set_node(mesh,8,node);
    int elemtype[] = {
        361,
    };
    int elem[] = {
        1,2,4,3,5,6,8,7
    };
    int sectID[] = {
        1,
    };
    hecmw_api_mesh_set_element(mesh,1,elemtype,elem,sectID);

    assert(hecmw_api_mesh_n_node(mesh)==8);
    assert(hecmw_api_mesh_n_elem(mesh)==1);
    assert(hecmw_api_mesh_n_elem_node_item(mesh)==8);
    assert(hecmw_api_mesh_n_ngrp(mesh)==0);
    assert(hecmw_api_mesh_n_sgrp(mesh)==0);
    assert(hecmw_api_mesh_n_egrp(mesh)==0);

    int ngrp[] = {
        1,2,3,4,
    };
    hecmw_api_mesh_append_ngrp(mesh,"N1",4,ngrp);
    assert(hecmw_api_mesh_n_ngrp(mesh)==1);
    char buf[64];
    hecmw_api_mesh_get_ngrp_name(mesh,1,buf,64);
    assert(strcmp("N1",buf)==0);
    int* array;
    int count=0;
    hecmw_api_mesh_get_ngrp(mesh,1,&array,&count);
    assert(count==4);
    assert(array[0]==1);
    assert(array[1]==2);
    assert(array[2]==3);
    assert(array[3]==4);
    hecmw_api_mesh_delete(mesh);
    return 0;
}

int main()
{
    mesh_test();
    surface_group_test();
    node_group_test();

    return 0;
}
