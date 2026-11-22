/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#undef NDEBUG
#include <stdio.h>
#include <assert.h>
#include "api/hecmw_api_local_mesh.h"
#include "api/fstr_api_param.h"
#include "api/fstr_api_solid.h"

int main(void)
{
    printf("FrontISTR API Test : fstr_solid\n");
    void* mesh1 = hecmw_api_mesh_new();
    void* mesh2 = hecmw_api_mesh_new();
    void* param1 = fstr_api_param_new();
    void* param2 = fstr_api_param_new();
    fstr_api_param_init(param1,mesh1);
    fstr_api_param_init(param2,mesh2);
    void* solid1 = fstr_api_solid_new();
    void* solid2 = fstr_api_solid_new();
    fstr_api_solid_init(solid1,mesh1,param1);
    assert(fstr_api_solid_nmaterial_tot(solid1)==0);
    fstr_api_solid_init(solid2,mesh2,param2);
    fstr_api_solid_delete(solid1);
    fstr_api_solid_delete(solid2);
    fstr_api_param_delete(param1);
    fstr_api_param_delete(param2);
    hecmw_api_mesh_delete(mesh1);
    hecmw_api_mesh_delete(mesh2);
    return 0;
}