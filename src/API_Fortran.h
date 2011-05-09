//
//  API_Fortran.h ( Fortran & C )
//
//
//                     2009.4.20
//                     2009.3.26
//                     k.Takeda

#ifndef API_Fortran_AF05492D_2BFF_49fd_AEE3_818B0B426793
#define API_Fortran_AF05492D_2BFF_49fd_AEE3_818B0B426793

#ifdef __cplusplus
extern "C" {
#endif

//----
// HEC_MW3 construct & destruct
//----
int mw_initialize_(int* argc, char** argv, char* path);// for C
int mw_initialize_1_(char* argv1, int* argv1_len, char* path, int* path_len);
int mw_initialize_2_(char* argv1, int* argv1_len, char* argv2, int* argv2_len, char* path, int* path_len);

int mw_finalize_();

//----
// file i/o API
//----
int mw_file_read_();
int mw_file_write_();

//----
// linear solver API
//----
//int mw_initialize_matrix_();
//int mw_initialize_vector_();
void mw_gene_linear_algebra_(int* num_of_algebra, int dof[]);
void mw_select_algebra_(int* ieq);

int mw_matrix_add_elem_(int* imesh, int* ielem,  double elem_matrix[]);// standard
int mw_matrix_add_node_(int* imesh, int* i_nid, int* j_nid, double nodal_matrix[]);

void mw_matrix_clear_(int* imesh);// matrix 0 clear
void mw_vector_clear_(int* imesh);// matrix 0 clear

int mw_matrix_add_elem_24_(int* imesh, int* ielem, double elem_matrix[][24]);//Hexa   8Node * 3DOF, Quad 8Node * 3DOF, Quad 4Node * 6DOF
int mw_matrix_add_elem_60_(int* imesh, int* ielem, double elem_matrix[][60]);//Hexa  20Node * 3DOF
int mw_matrix_add_elem_12_(int* imesh, int* ielem, double elem_matirx[][12]);//Tetra  4Node * 3DOF, Quad 4Node * 3DOF, Beam 2Node * 6DOF
int mw_matrix_add_elem_30_(int* imesh, int* ielem, double elem_matirx[][30]);//Tetra 10Node * 3DOF, Tri  6Node * 5DOF
int mw_matrix_add_elem_18_(int* imesh, int* ielem, double elem_matirx[][18]);//Prism  6Node * 3DOF, Tri  6Node * 3DOF, Beam 3Node * 6DOF
int mw_matirx_add_elem_45_(int* imesh, int* ielem, double elem_matirx[][45]);//Prism 15Node * 3DOF
int mw_matirx_add_elem_20_(int* imesh, int* ielem, double elem_matirx[][20]);//Quad   4Node * 5DOF
int mw_matrix_add_elem_40_(int* imesh, int* ielem, double elem_matirx[][40]);//Quad   8Node * 5DOF
int mw_matrix_add_elem_15_(int* imesh, int* ielem, double elem_matirx[][15]);//Tri    3Node * 5DOF, Beam 3Node * 5DOF
int mw_matirx_add_elem_9_(int* imesh, int* ielem, double elem_matirx[][9]);  //Tri    3Node * 3DOF, Beam 3Node * 3DOF
int mw_matirx_add_elem_48_(int* imesh, int* ielem, double elem_matirx[][48]);//Quad   8Node * 6DOF
int mw_matirx_add_elem_6_(int* imesh, int* ielem, double elem_matirx[][6]);  //Beam   2Node * 3DOF
int mw_matirx_add_elem_10_(int* imesh, int* ielem, double elem_matirx[][10]);//Beam   2Node * 5DOF

int mw_matrix_set_bc_(int* imesh, int* inode, int* idof, double* value1, double* value2);// matrix-D and solution_vector
int mw_matrix_rhs_set_bc_(int* imesh, int* inode, int* idof, double* value1, double* value2);// matirix-D and rhs_vector
int mw_rhs_set_bc_(int* imesh, int* inode, int* idof, double* value);

int mw_solve_(int* iter_max, double* tolerance, int* method, int* pre_condition);

//void mw_store_matrix_();
//void mw_load_matrix_();


//----
// MG construct (refine)
//----
int mw_refine_();      // refine == mg_construct
int mw_mg_construct_();// mg_construct == refine
void mw_finalize_refine_();      // release memory (final proc for mesh construct)
void mw_finalize_mg_construct_();// release memory (final proc for mesh construct) == finalize_refine

//----
// model
//----
// assemble model
int mw_get_num_of_assemble_model_();
void mw_select_assemble_model_(int* mglevel);
// mesh part
int mw_get_num_of_mesh_part_();
void mw_select_mesh_part_with_id_(int* mesh_id);
void mw_select_mesh_part_(int* index);
// element
void mw_select_element_with_id_(int* elem_id);
void mw_select_element_(int* index);
int mw_get_element_type_();
int mw_get_num_of_element_vert_();

void mw_get_element_vert_node_id_(int v_node_id[]);
int mw_get_num_of_element_edge_();
void mw_get_element_edge_node_id_(int v_node_id[]);

// node
void mw_get_node_coord_(int* node_id, double* x, double* y, double* z);
int mw_get_dof_(int* node_id);
int mw_get_dof_scalar_(int* node_id);
int mw_get_dof_vector_(int* node_id);

void mw_set_node_value_(int* node_id, double value[]);
void mw_set_node_value_with_dof_(int* node_id, int* idof, double* value);
void mw_get_node_value_(int* node_id, double value[]);
void mw_get_node_value_with_dof_(int* node_id, int* idof, double* value);

void mw_set_sv_node_value_(int* node_id, double v_value[], double s_value[]);
void mw_set_sv_node_value_with_dof_(int* node_id, int* v_dof, double* v_value, int* s_dof, double* s_value);
void mw_get_sv_node_value_(int* node_id, double v_value[], double s_value[]);
void mw_get_sv_node_value_with_dof_(int* node_id, int* v_dof, double* v_value, int* s_dof, double* s_value);

// node size, element size
int mw_get_num_of_node_();// in select mesh_part
int mw_get_num_of_node_with_mesh_(int* imesh);
int mw_get_num_of_element_();// in select mesh_part
int mw_get_num_of_element_with_mesh_(int* imesh);

// id
int mw_get_element_id_(int* index);
int mw_get_node_id_(int* index);

//----
// node type
//----
int mw_nodetype_s_();
int mw_nodetype_v_();
int mw_nodetype_sv_();
//----
// element type
//----
int mw_elemtype_hexa_();
int mw_elemtype_hexa2_();
int mw_elemtype_tetra_();
int mw_elemtype_tetra2_();
int mw_elemtype_prism_();
int mw_elemtype_prism2_();
int mw_elemtype_quad_();
int mw_elemtype_quad2_();
int mw_elemtype_triangle_();
int mw_elemtype_triangle2_();
int mw_elemtype_line_();
int mw_elemtype_line2_();
//----
// frontISTR element type
//----
int mw_fistr_elemtype_hexa_();
int mw_fistr_elemtype_hexa2_();
int mw_fistr_elemtype_tetra_();
int mw_fistr_elemtype_tetra2_();
int mw_fistr_elemtype_prism_();
int mw_fistr_elemtype_prism2_();
int mw_fistr_elemtype_quad_();
int mw_fistr_elemtype_quad2_();
int mw_fistr_elemtype_triangle_();
int mw_fistr_elemtype_triangle2_();
int mw_fistr_elemtype_line_();
int mw_fistr_elemtype_line2_();
//----
// frontISTR element type => MW3 element type
//----
int mw_fistr_elemtype_to_mw3_elemtype_(int* fistr_elemtype);
//----
// MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
//----
int mw_mw3_elemtype_to_fistr_elemtype_(int* mw3_elemtype);


//----
// shape function
//----
int mw_get_num_of_integ_point_(int* shape_type);
void mw_shape_function_on_pt_(int* shape_type, int* igauss, double N[]);

void mw_shape_function_hexa81_(int* igauss, int* ishape, double* N);
void mw_shape_function_hexa82_(int* igauss, int* ishape, double* N);
void mw_shape_function_hexa201_(int* igauss, int* ishape, double* N);
void mw_shape_function_hexa202_(int* igauss, int* ishape, double* N);
void mw_shape_function_hexa203_(int* igauss, int* ishape, double* N);
void mw_shape_function_tetra41_(int* igauss, int* ishape, double* N);
void mw_shape_function_tetra101_(int* igauss, int* ishape, double* N);
void mw_shape_function_tetra104_(int* igauss, int* ishape, double* N);
void mw_shape_function_tetra1015_(int* igauss, int* ishape, double* N);
void mw_shape_function_prism62_(int* igauss, int* ishape, double* N);
void mw_shape_function_prism156_(int* igauss, int* ishape, double* N);
void mw_shape_function_prism159_(int* igauss, int* ishape, double* N);
void mw_shape_function_prism1518_(int* igauss, int* ishape, double* N);
void mw_shape_function_quad41_(int* igauss, int* ishape, double* N);
void mw_shape_function_quad84_(int* igauss, int* ishape, double* N);
void mw_shape_function_quad89_(int* igauss, int* ishape, double* N);
void mw_shape_function_tri31_(int* igauss, int* ishape, double* N);
void mw_shape_function_tri63_(int* igauss, int* ishape, double* N);
void mw_shape_function_line21_(int* igauss, int* ishape, double* N);
void mw_shape_function_line32_(int* igauss, int* ishape, double* N);

//----
// shape function deriv (rst coord)
//----
void mw_dndr_(int* shape_type, double dndr[]);
void mw_dndr_hexa81_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_hexa82_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_hexa201_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_hexa202_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_hexa203_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_tetra41_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_tetra101_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_tetra104_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_tetra1015_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_prism62_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_prism156_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_prism159_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_prism1518_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_quad41_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_quad84_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_quad89_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_tri31_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_tri63_(int* igauss, int* ishape, int* iaxis, double* dndr);
void mw_dndr_line21_(int* igauss, int* ishape, double* dndr);
void mw_dndr_line32_(int* igauss, int* ishape, double* dndr);
//----
// shape function deriv (xyz coord)
//----
void mw_dndx_(int* elem_type, int* num_of_integ, int* ielem, double dndx[]);
void mw_det_jacobian_(int* elem_type, int* num_of_integ, int* igauss, double* det_j);
void mw_weight_(int* elem_type, int* num_of_integ, int* igauss, double* w);

//----
// shape function type
//----
int mw_shapetype_hexa81_();
int mw_shapetype_hexa82_();
int mw_shapetype_hexa201_();
int mw_shapetype_hexa202_();
int mw_shapetype_hexa203_();
int mw_shapetype_tetra41_();
int mw_shapetype_tetra101_();
int mw_shapetype_tetra104_();
int mw_shapetype_tetra1015_();
int mw_shapetype_prism62_();
int mw_shapetype_prism156_();
int mw_shapetype_prism159_();
int mw_shapetype_prism1518_();
int mw_shapetype_quad41_();
int mw_shapetype_quad84_();
int mw_shapetype_quad89_();
int mw_shapetype_tri31_();
int mw_shapetype_tri63_();
int mw_shapetype_line21_();
int mw_shapetype_line32_();


//--
// boundary mesh
//--
int mw_get_num_of_boundary_bnode_mesh_();
int mw_get_num_of_boundary_bface_mesh_();
int mw_get_num_of_boundary_bedge_mesh_();
int mw_get_num_of_boundary_bvolume_mesh_();
int mw_get_num_of_bnode_in_bnode_mesh_(int* ibmesh);
int mw_get_num_of_bnode_in_bface_mesh_(int* ibmesh);
int mw_get_num_of_bnode_in_bedge_mesh_(int* ibmesh);
int mw_get_num_of_bnode_in_bvolume_mesh_(int* ibmesh);
int mw_get_num_of_dof_in_bnode_mesh_(int* ibmesh, int* ibnode);
int mw_get_num_of_dof_in_bface_mesh_(int* ibmesh);
int mw_get_num_of_dof_in_bedge_mesh_(int* ibmesh);
int mw_get_num_of_dof_in_bvolume_mesh_(int* ibmesh);
//--
// value of boundary node
//--
double mw_get_bnode_value_in_bnode_mesh_(int* ibmesh, int* ibnode, int* idof);
double mw_get_bnode_value_in_bface_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel);
double mw_get_bnode_value_in_bedge_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel);
double mw_get_bnode_value_in_bvolume_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel);
int mw_get_node_id_in_bnode_mesh_(int* ibmesh, int* ibnode);
int mw_get_node_id_in_bface_mesh_(int* ibmesh, int* ibnode);
int mw_get_node_id_in_bedge_mesh_(int* ibmesh, int* ibnode);
int mw_get_node_id_in_bvolume_mesh_(int* ibmesh, int* ibnode);
//--
// value of boundary face, edge, volume
//--
int mw_get_num_of_bface_(int* ibmesh);
double mw_get_bface_value_(int* ibmesh, int* ibface, int* idof);
int mw_get_num_of_bedge_(int* ibmesh);
double mw_get_bedge_value_(int* ibmesh, int* ibedge, int* idof);
int mw_get_num_of_bvolume_(int* ibmesh);
double mw_get_bvolume_value_(int* ibmesh, int* ibvol, int* idof);
//--
// boundary_mesh name
//--
int mw_get_bnode_mesh_namelength_(int* ibmesh);
void mw_get_bnode_mesh_name_(int* ibmesh, char* name, int* name_len);
int mw_get_bface_mesh_namelength_(int* ibmesh);
void mw_get_bface_mesh_name_(int* ibmesh, char* name, int* name_len);
int mw_get_bvolume_mesh_namelength_(int* ibmesh);
void mw_get_bvolume_mesh_name_(int* ibmesh, char* name, int* name_len);
int mw_get_bedge_mesh_namelength_(int* ibmesh);
void mw_get_bedge_mesh_name_(int* ibmesh, char* name, int* name_len);

//--
// mpi
//--
int mw_mpi_int_();   // MPI_INT
int mw_mpi_double_();// MPI_DOUBLE
int mw_mpi_comm_();  // MPI_COMM_WORLD

int mw_mpi_sum_();// op  ,use allreduce_r argument
int mw_mpi_max_();// op  ,use allreduce_r argument
int mw_mpi_min_();// op  ,use allreduce_r argument

void mw_allreduce_r_(double val[], int* val_size, int* op);
void mw_allreduce_i_(int val[], int* val_size, int* op);
int mw_barrier_();
int mw_abort_(int* error);
int mw_allgather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt);
int mw_allgather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt);
int mw_gather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root);
int mw_gather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root);
int mw_scatter_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root);
int mw_scatter_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root);
int mw_get_rank_();
int mw_get_num_of_process_();
void mw_send_recv_r2_(double buf[], int* dof_size);// bufの値を送信, 受信値をNodeとbufに代入. bufのサイズ == NumOfCommNode * dof_size
void mw_send_recv_r_();// 通信Nodeの値を入れ替えて更新


//--
// Element_Group { select AssyModel, select Mesh }
//--
int mw_get_num_of_elementgroup_();
int mw_get_num_of_element_id_(int* igrp);
int mw_get_element_id_with_elementgroup_(int* igrp, int* index);
int mw_get_elementgroup_name_length_(int* igrp);
void mw_get_elementgroup_name_(int* igrp, char* name, int* name_len);

//----
// logger
//----
void mw_logger_set_mode_(int* mode);
void mw_logger_set_device_(int* mode, int* device);
void mw_logger_info_ (int* mode, char* message, int* str_len);
//----
// logger parameter
//----
int mw_get_error_mode_();
int mw_get_warn_mode_();
int mw_get_info_mode_();
int mw_get_debug_mode_();
int mw_get_disk_device_();
int mw_get_display_device_();

#ifdef __cplusplus
}
#endif

#endif // API_Fortran.h

