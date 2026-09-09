/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
#ifndef hecmw_api_result_dataH
#define hecmw_api_result_dataH

void* hecmw_api_result_new();
void hecmw_api_result_delete(void* result);

int hecmw_api_result_ng_component(void* result);
void hecmw_api_result_global_val(void* result,int i,char* label,int label_len,double *value);

int hecmw_api_result_nn_component(void* result);
void hecmw_api_result_node_val(void* result,int i,int* dim,int* index,int* dof,char* label,int* label_len,double** value);

int hecmw_api_result_ne_component(void* result);
void hecmw_api_result_elem_val(void* result,int i,int* dim,int* index,int* dof,char* label,int* label_len,double** value);

#endif