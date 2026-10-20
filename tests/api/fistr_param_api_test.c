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
    void* mesh = hecmw_api_mesh_new();
    void* param = fstr_api_param_new();
    fstr_api_param_init(param,mesh);

    fstr_api_param_set_solution_type(param,kstDYNAMIC);
    assert(fstr_api_param_solution_type(param)==kstDYNAMIC);

    fstr_api_param_set_solution_type(param,kstEIGEN);
    assert(fstr_api_param_solution_type(param)==kstEIGEN);

    fstr_api_param_set_solver_method(param,ksmGMRES);
    assert(fstr_api_param_solver_method(param)==ksmGMRES);

    fstr_api_param_delete(param);
    hecmw_api_mesh_delete(mesh);
    return 0;
}