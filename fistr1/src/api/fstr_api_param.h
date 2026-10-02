/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef hecmw_api_local_meshH
#define hecmw_api_local_meshH

#include <stdbool.h>

typedef enum {
  kstPRECHECK    = 0,
  kstSTATIC      = 1,
  kstEIGEN       = 2,
  kstHEAT        = 3,
  kstDYNAMIC     = 4,
  kstSTATICEIGEN = 6,
  kstNZPROF      = 7,
} fstr_solution_type;

typedef enum {
  ksmCG       = 1,
  ksmBiCGSTAB = 2,
  ksmGMRES    = 3,
  ksmGPBiCG   = 4,
  ksmGMRESR   = 5,
  ksmGMRESREN = 6,
  ksmPipeCG   = 8,
  ksmGroppCG  = 9,
  ksmDIRECT   = 101,
} fstr_solver_method;

void* fstr_api_param_new();
void fstr_api_param_delete(void* param);
void fstr_api_param_init(void* param,void* mesh);
fstr_solution_type fstr_api_param_solution_type(void* param);
fstr_solver_method fstr_api_param_solver_method(void* param);
bool fstr_api_param_nlgeom(void* param);
int fstr_api_param_fg_result(void* param);
int fstr_api_param_fg_visual(void* param);

#endif
