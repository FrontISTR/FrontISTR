/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef fstr_api_solidH
#define fstr_api_solidH

typedef struct fstr_solid_handle_t* fstr_solid_handle;

fstr_solid_handle fstr_api_solid_new(void);
void fstr_api_solid_delete(fstr_solid_handle solid);
void fstr_api_solid_init(fstr_solid_handle solid,void* mesh,void* param);

int fstr_api_solid_nmaterial_tot(fstr_solid_handle solid);
void fstr_api_solid_init_material(fstr_solid_handle solid,int n);

#endif