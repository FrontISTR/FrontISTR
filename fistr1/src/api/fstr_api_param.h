/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef hecmw_api_local_meshH
#define hecmw_api_local_meshH

void* fstr_api_param_new();
void fstr_api_param_delete(void* param);
void fstr_api_param_init(void* param,void* mesh);
int fstr_api_param_solutuin_type(void* param);

#endif
