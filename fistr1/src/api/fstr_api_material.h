/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
#ifndef fstr_api_materialH
#define fstr_api_materialH

#include <stdbool.h>
#include "api/fstr_api_solid.h"

typedef enum {
  INFINITESIMAL = 0,
  TOTALLAG  = 1,
  UPDATELAG = 2,
  CAUCHY = 2,
} fstr_material_nlgeom;

typedef enum {
  MATERIAL_NOT_DEFINED  = -1,
  USERMATERIAL          = 100000,
  ELASTIC               = 110000,
  MN_ORTHOELASTIC       = 111000,
  USERELASTIC           = 112000,
  EPLASTIC              = 120000,
  NEOHOOKE              = 130000,
  MOONEYRIVLIN          = 131000,
  ARRUDABOYCE           = 132000,
  USERHYPERELASTIC      = 133000,
  MOONEYRIVLIN_ANISO    = 134000,
  VISCOELASTIC          = 140000,
  NORTON                = 150000,
  INCOMP_NEWTONIAN      = 160000,
  CONNECTOR             = 170000,
} fstr_material_type;

typedef struct fstr_material_handle_t* fstr_material_handle;

fstr_material_handle fstr_api_solid_material(void* solid,int cid);
fstr_material_type fstr_api_material_type(fstr_material_handle material);
void fstr_api_material_get_name(fstr_material_handle material,char* name,int len);
void fstr_api_material_set_name(fstr_material_handle material,char* name);

bool fstr_api_material_isElastic(fstr_material_handle material);
bool fstr_api_material_isElastoplastic(fstr_material_handle material);
bool fstr_api_material_isHyperelastic(fstr_material_handle material);
bool fstr_api_material_isViscoelastic(fstr_material_handle material);
bool fstr_api_material_isCreep(fstr_material_handle material);

fstr_material_nlgeom fstr_api_material_nlgeom(void* material);
void fstr_api_material_set_nlgeom(fstr_material_handle material,fstr_material_nlgeom nlgeom);

void fstr_api_material_get_elastic(fstr_material_handle material,double* youngs,double* poisson);
void fstr_api_material_set_elastic(fstr_material_handle material,double youngs,double poisson);

double fstr_api_material_density(fstr_material_handle material);
void fstr_api_material_set_density(fstr_material_handle material,double density);

double fstr_api_material_expansion(fstr_material_handle material);
void fstr_api_material_set_expansion(fstr_material_handle material,double expansion);

void fstr_api_material_get_plastic_size(fstr_material_handle material,char* yield,char* harden,int* nvar,int* rows,int* cols);
void fstr_api_material_get_plastic_vars(fstr_material_handle material,double* variables);
void fstr_api_material_set_plastic(fstr_material_handle material,char* yield,char* harden,double* variables,int rows,int cols);

void fstr_api_material_get_hyperelastic(fstr_material_handle material,char* type,double* variables);
void fstr_api_material_set_hyperelastic(fstr_material_handle material,char* type,double* variables);

void fstr_api_material_get_viscoelastic_size(fstr_material_handle material,int* rows,int* cols);
void fstr_api_material_get_viscoelastic_vars(fstr_material_handle material,double* variables);
void fstr_api_material_set_viscoelastic(fstr_material_handle material,double* variables,int rows,int cols);

void fstr_api_material_get_creep_size(fstr_material_handle material,char* type,int* rows,int* cols);
void fstr_api_material_get_creep_vars(fstr_material_handle material,double* variables);
void fstr_api_material_set_creep(fstr_material_handle material,char* type,double* variables,int rows,int cols);

#endif
