/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#undef NDEBUG
#include <stdio.h>
#include <assert.h>
#include "api/fstr_api_param.h"
#include "api/hecmw_api_local_mesh.h"

int main(void)
{
    printf("FrontISTR API Test : fstr_param\n");
    void* mesh1 = hecmw_api_mesh_new();
    void* mesh2 = hecmw_api_mesh_new();
    void* param1 = fstr_api_param_new();
    void* param2 = fstr_api_param_new();
    fstr_api_param_init(param1,mesh1);
    fstr_api_param_init(param2,mesh2);

    fstr_api_param_set_solution_type(param1,kstDYNAMIC);
    assert(fstr_api_param_solution_type(param1)==kstDYNAMIC);

    fstr_api_param_set_solution_type(param2,kstEIGEN);
    assert(fstr_api_param_solution_type(param2)==kstEIGEN);

    fstr_api_param_set_solver_method(param1,ksmGMRES);
    assert(fstr_api_param_solver_method(param1)==ksmGMRES);

    fstr_api_param_set_solver_method(param2,ksmCG);
    assert(fstr_api_param_solver_method(param2)==ksmCG);

    fstr_api_param_set_ref_temp(param2,123.0);
    assert(fstr_api_param_ref_temp(param2)==123.0);

    fstr_api_param_delete(param1);
    fstr_api_param_delete(param2);
    hecmw_api_mesh_delete(mesh1);
    hecmw_api_mesh_delete(mesh2);
    return 0;
}
