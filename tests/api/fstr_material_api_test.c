/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#undef NDEBUG
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "api/hecmw_api_local_mesh.h"
#include "api/fstr_api_param.h"
#include "api/fstr_api_solid.h"
#include "api/fstr_api_material.h"

int main(void)
{
    printf("FrontISTR API Test : fstr_material\n");
    void* mesh1 = hecmw_api_mesh_new();
    void* mesh2 = hecmw_api_mesh_new();
    void* param1 = fstr_api_param_new();
    void* param2 = fstr_api_param_new();
    fstr_api_param_init(param1,mesh1);
    fstr_api_param_init(param2,mesh2);
    void* solid1 = fstr_api_solid_new();
    void* solid2 = fstr_api_solid_new();
    fstr_api_solid_init_material(solid1,4);
    // Fortran index (begin with 1)
    fstr_material_handle mat1 = fstr_api_solid_material(solid1,1);
    fstr_material_handle mat2 = fstr_api_solid_material(solid1,2);
    fstr_material_handle mat3 = fstr_api_solid_material(solid1,3);
    fstr_material_handle mat4 = fstr_api_solid_material(solid1,4);
    assert(fstr_api_material_type(mat1)==MATERIAL_NOT_DEFINED);
    assert(fstr_api_material_type(mat2)==MATERIAL_NOT_DEFINED);
    assert(fstr_api_material_type(mat3)==MATERIAL_NOT_DEFINED);
    assert(fstr_api_material_type(mat4)==MATERIAL_NOT_DEFINED);
    assert(fstr_api_material_nlgeom(mat1)==INFINITESIMAL);

    // elastic
    fstr_api_material_set_name(mat1,"MAT1");
    char buf[32];
    fstr_api_material_get_name(mat1,buf,32);
    assert(strcmp("MAT1",buf)==0);
    fstr_api_material_set_elastic(mat1,1000.0,0.3);
    assert(fstr_api_material_type(mat1)==ELASTIC);
    assert(fstr_api_material_nlgeom(mat1)==INFINITESIMAL);
    fstr_api_material_set_nlgeom(mat1,TOTALLAG);
    assert(fstr_api_material_nlgeom(mat1)==TOTALLAG);
    double youngs,poisson;
    fstr_api_material_get_elastic(mat1,&youngs,&poisson);
    assert(youngs==1000.0);
    assert(poisson==0.3);
    assert(fstr_api_material_isElastic(mat1));
    assert(fstr_api_material_isElastoplastic(mat1)==false);
    assert(fstr_api_material_isHyperelastic(mat1)==false);
    assert(fstr_api_material_isViscoelastic(mat1)==false);
    assert(fstr_api_material_isCreep(mat1)==false);

    // density & expansion
    fstr_api_material_set_density(mat1,1234.0);
    assert(fstr_api_material_density(mat1)==1234.0);
    fstr_api_material_set_expansion(mat1,0.567);
    assert(fstr_api_material_expansion(mat1)==0.567);

    // plastic
    fstr_api_material_set_name(mat2,"MAT2");
    double plas_var[] = {
        250.0, 0.0, 20.0,
        300.0, 0.1, 20.0,
        200.0, 0.0, 400.0,
        240.0, 0.1, 400.0,
    };
    fstr_api_material_set_plastic(mat2,"MISES","MULTILINEAR",plas_var,4,3);
    assert(fstr_api_material_nlgeom(mat2)==INFINITESIMAL);
    fstr_api_material_set_nlgeom(mat2,UPDATELAG);
    assert(fstr_api_material_nlgeom(mat2)==UPDATELAG);
    char yield[32], harden[32];
    int nvars, rows, cols;
    fstr_api_material_get_plastic_size(mat2,yield,harden,&nvars,&rows,&cols);
    assert(nvars==0);
    assert(rows==4);
    assert(cols==3);
    double var[100];
    fstr_api_material_get_plastic_vars(mat2,var);
    for(int i=0;i<rows*cols;i++){
        assert(var[i]==plas_var[i]);
    }
    assert(fstr_api_material_isElastoplastic(mat2));

    // viscoelastic
    double visco_var[] = {
        0.01, 0.12,
        0.1, 0.23,
        0.2, 0.34,
        0.4, 0.56,
        0.6, 0.78,
        0.7, 0.89
    };
    fstr_api_material_set_viscoelastic(mat1,visco_var,6,2);
    fstr_api_material_get_viscoelastic_size(mat1,&rows,&cols);
    assert(rows==6);
    assert(cols==2);
    fstr_api_material_get_viscoelastic_vars(mat1,var);
    for(int i=0;i<rows*cols;i++){
        assert(var[i]==visco_var[i]);
    }
    assert(fstr_api_material_isViscoelastic(mat1));

    // hyperelastic
    fstr_api_solid_init_material(solid2,1);
    // Fortran index (begin with 1)
    fstr_material_handle mat5 = fstr_api_solid_material(solid2,1);
    fstr_api_material_set_name(mat5,"MAT5");
    double hyper_var[] = {
        0.71, 1.7029, 0.1408
    };
    fstr_api_material_set_hyperelastic(mat5,"ARRUDA-BOYCE",hyper_var);
    assert(fstr_api_material_nlgeom(mat5)==INFINITESIMAL);
    fstr_api_material_set_nlgeom(mat5,TOTALLAG);
    assert(fstr_api_material_nlgeom(mat5)==TOTALLAG);
    char type[32];
    fstr_api_material_get_hyperelastic(mat5,type,var);
    assert(strcmp("ARRUDA-BOYCE",type)==0);
    for(int i=0;i<3;i++){
        assert(var[i]==hyper_var[i]);
    }

    fstr_api_solid_delete(solid1);
    fstr_api_solid_delete(solid2);
    fstr_api_param_delete(param1);
    fstr_api_param_delete(param2);
    hecmw_api_mesh_delete(mesh1);
    hecmw_api_mesh_delete(mesh2);
    return 0;
}
