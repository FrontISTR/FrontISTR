//----
// for linking HEC_MW3
//----

int mw_initialize_(int* argc, char** argv);
int mw_initialize_fstr_(int* argc, char** argv, char* cntrolfile);
int mw_finalize_();

int mw_file_read_();
int mw_file_read_fstr_();

int mw_refine_();
int mw_mg_construct_();
void mw_finalize_refine_();
void mw_finalize_mg_construct_();

int mw_get_num_of_assemble_model_();
void mw_select_assemble_model_(int* mglevel);
int mw_get_num_of_mesh_part_();
void mw_select_mesh_part_with_id_(int* mesh_id);
void mw_select_mesh_part_(int* index);

void mw_select_element_with_id_(int* elem_id);
void mw_select_element_(int* index);
int mw_get_element_type_();
int mw_get_num_of_element_vert_();
void mw_get_element_vert_node_id_(int v_node_id[]);
int mw_get_num_of_element_edge_();
void mw_get_element_edge_node_id_(int v_node_id[]);

void mw_get_node_coord_(int* node_id, double* x, double* y, double* z);
int mw_get_dof_(int* node_id);
int mw_get_dof_scalar_(int* node_id);
int mw_get_dof_vector_(int* node_id);

int mw_get_num_of_node_();
int mw_get_num_of_node_with_mesh_(int* imesh);
int mw_get_num_of_element_();
int mw_get_num_of_element_with_mesh_(int* imesh);

int mw_get_node_id_(int* index);
int mw_get_element_id_(int* index);
int mw_get_node_index_(int* id);
int mw_get_element_index_(int* id);

int mw_fistr_elemtype_to_mw3_elemtype_(int* fistr_elemtype);
int mw_mw3_elemtype_to_fistr_elemtype_(int* mw3_elemtype);

#ifdef HAVE_MPI
MPI_Comm mw_mpi_comm_();
#else
int mw_mpi_comm_();
#endif
int mw_get_rank_();
int mw_get_num_of_process_();

int mw_get_num_of_neibpe_(int* imesh);
int mw_get_transrank_(int* imesh, int* ipe);

int mw_get_num_of_comm_mesh_();
int mw_get_num_of_comm_node_(int* icmesh);
int mw_get_node_id_comm_node_(int* icmesh, int* icnode);
