#ifndef hecmw_local_mesh_apiH
#define hecmw_local_mesh_apiH

void* hecmw_api_mesh_new();
void hecmw_api_mesh_delete(void* mesh);

void hecmw_api_mesh_set_node(void* mesh,int nnode,const double* node);
int hecmw_api_mesh_n_node(const void* mesh);
const double* hecmw_api_mesh_get_node(const void* mesh);

void hecmw_api_mesh_set_n_dof(void* mesh,int ndof);
int hecmw_api_mesh_get_n_dof(const void* mesh);

void hecmw_api_mesh_set_element(void* mesh,int nelem,const int* elemtype,const int* element,const int* sectionID);
const int* hecmw_api_mesh_get_elem_type(const void* mesh);
const int* hecmw_api_mesh_get_elem_node_item(const void* mesh);
const int* hecmw_api_mesh_get_section_id(const void* mesh);
int hecmw_api_mesh_n_elem(const void* mesh);
int hecmw_api_mesh_n_elem_node_item(const void* mesh);

int hecmw_api_mesh_n_ngrp(const void* mesh);
void hecmw_api_mesh_append_ngrp(void* mesh,const char* grp_name,int count,int* list);
const char* hecmw_api_mesh_get_ngrp_name(const void* mesh,int i);
void hecmw_api_mesh_get_ngrp(const void* mesh,int i,int* array,int count);

int hecmw_api_mesh_n_sgrp(const void* mesh);
void hecmw_api_mesh_append_sgrp(void* mesh,const char* grp_name,int count,int* list);
const char* hecmw_api_mesh_get_sgrp_name(const void* mesh,int i);
void hecmw_api_mesh_get_sgrp(const void* mesh,int i,int* array,int count);

int hecmw_api_mesh_n_egrp(const void* mesh);
void hecmw_api_mesh_append_egrp(void* mesh,const char* grp_name,int count,int* list);
const char* hecmw_api_mesh_get_egrp_name(const void* mesh,int i);
void hecmw_api_mesh_get_egrp(const void* mesh,int i,int* array,int count);

#endif