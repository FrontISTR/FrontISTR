#ifndef hecmw_local_mesh_apiH
#define hecmw_local_mesh_apiH

void* hecmw_api_mesh_new();
void hecmw_api_mesh_delete(void* mesh);
void hecmw_api_mesh_set_node(void* mesh,int nnode,const double* node);
int hecmw_api_mesh_n_node(const void* mesh);
const double* hecmw_api_mesh_get_node(const void* mesh);

void hecmw_api_mesh_set_n_dof(void* mesh,int ndof);
int hecmw_api_mesh_get_n_dof(void* mesh);

#endif