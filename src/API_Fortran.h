/* 
 * File:   API_Fortran.hxx
 * Author: ktakeda
 *
 * Created on 2010/10/04, 19:18
 */

#ifndef API_FORTRAN_HXX_VISUAL_CPP
#define	API_FORTRAN_HXX_VISUAL_CPP


#ifdef __cplusplus
extern "C" {
#endif

//----
// HEC_MW3 construct & destruct
//----
__declspec(dllexport) int mw_initialize_(int* argc, char** argv, char* path);// for C
__declspec(dllexport) int mw_initialize_1_(char* argv1, int* argv1_len, char* path, int* path_len);
__declspec(dllexport) int mw_initialize_2_(char* argv1, int* argv1_len, char* argv2, int* argv2_len, char* path, int* path_len);

__declspec(dllexport) int mw_finalize_();

//----
// file i/o API
//----
__declspec(dllexport) int mw_file_read_();
__declspec(dllexport) int mw_file_write_();

//----
// linear solver API
//----
__declspec(dllexport) int mw_initialize_matrix_();
__declspec(dllexport) int mw_initialize_vector_();

__declspec(dllexport) int mw_matrix_add_elem_(int* imesh,  int* ielem,  double elem_matrix[]);// standard

__declspec(dllexport) int mw_matrix_add_elem_24_(int* imesh, int* ielem, double elem_matrix[][24]);//Hexa   8Node * 3DOF, Quad 8Node * 3DOF, Quad 4Node * 6DOF
__declspec(dllexport) int mw_matrix_add_elem_60_(int* imesh, int* ielem, double elem_matrix[][60]);//Hexa  20Node * 3DOF
__declspec(dllexport) int mw_matrix_add_elem_12_(int* imesh, int* ielem, double elem_matirx[][12]);//Tetra  4Node * 3DOF, Quad 4Node * 3DOF, Beam 2Node * 6DOF
__declspec(dllexport) int mw_matrix_add_elem_30_(int* imesh, int* ielem, double elem_matirx[][30]);//Tetra 10Node * 3DOF, Tri  6Node * 5DOF
__declspec(dllexport) int mw_matrix_add_elem_18_(int* imesh, int* ielem, double elem_matirx[][18]);//Prism  6Node * 3DOF, Tri  6Node * 3DOF, Beam 3Node * 6DOF
__declspec(dllexport) int mw_matirx_add_elem_45_(int* imesh, int* ielem, double elem_matirx[][45]);//Prism 15Node * 3DOF
__declspec(dllexport) int mw_matirx_add_elem_20_(int* imesh, int* ielem, double elem_matirx[][20]);//Quad   4Node * 5DOF
__declspec(dllexport) int mw_matrix_add_elem_40_(int* imesh, int* ielem, double elem_matirx[][40]);//Quad   8Node * 5DOF
__declspec(dllexport) int mw_matrix_add_elem_15_(int* imesh, int* ielem, double elem_matirx[][15]);//Tri    3Node * 5DOF, Beam 3Node * 5DOF
__declspec(dllexport) int mw_matirx_add_elem_9_(int* imesh, int* ielem, double elem_matirx[][9]);  //Tri    3Node * 3DOF, Beam 3Node * 3DOF
__declspec(dllexport) int mw_matirx_add_elem_48_(int* imesh, int* ielem, double elem_matirx[][48]);//Quad   8Node * 6DOF
__declspec(dllexport) int mw_matirx_add_elem_6_(int* imesh, int* ielem, double elem_matirx[][6]);  //Beam   2Node * 3DOF
__declspec(dllexport) int mw_matirx_add_elem_10_(int* imesh, int* ielem, double elem_matirx[][10]);//Beam   2Node * 5DOF

__declspec(dllexport) int mw_matrix_set_bc_(int* imesh, int* inode, int* idof, double* value1, double* value2);
__declspec(dllexport) int mw_rhs_set_bc_(int* imesh, int* inode, int* idof, double* value);

__declspec(dllexport) int mw_solve_(int* iter_max, double* tolerance, int* method, int* pre_condition);

__declspec(dllexport) void mw_store_matrix_();
__declspec(dllexport) void mw_load_matrix_();


//----
// MG construct (refine)
//----
__declspec(dllexport) int mw_refine_();      // refine == mg_construct
__declspec(dllexport) int mw_mg_construct_();// mg_construct == refine
__declspec(dllexport) void mw_finalize_refine_();      // release memory (final proc for mesh construct)
__declspec(dllexport) void mw_finalize_mg_construct_();// release memory (final proc for mesh construct) == finalize_refine

//----
// model
//----
// assemble model
__declspec(dllexport) int mw_get_num_of_assemble_model_();
__declspec(dllexport) void mw_select_assemble_model_(int* mglevel);
// mesh part
__declspec(dllexport) int mw_get_num_of_mesh_part_();
__declspec(dllexport) void mw_select_mesh_part_with_id_(int* mesh_id);
__declspec(dllexport) void mw_select_mesh_part_(int* index);
// element
__declspec(dllexport) void mw_select_element_with_id_(int* elem_id);
__declspec(dllexport) void mw_select_element_(int* index);
__declspec(dllexport) int mw_get_element_type_();
__declspec(dllexport) int mw_get_num_of_element_vert_();

__declspec(dllexport) void mw_get_element_vert_node_id_(int v_node_id[]);
__declspec(dllexport) int mw_get_num_of_element_edge_();
__declspec(dllexport) void mw_get_element_edge_node_id_(int v_node_id[]);

// node
__declspec(dllexport) void mw_get_node_coord_(int* node_id, double* x, double* y, double* z);
__declspec(dllexport) int mw_get_dof_(int* node_id);
__declspec(dllexport) int mw_get_dof_scalar_(int* node_id);
__declspec(dllexport) int mw_get_dof_vector_(int* node_id);

__declspec(dllexport) void mw_set_node_value_(int* node_id, double value[]);
__declspec(dllexport) void mw_set_node_value_with_dof_(int* node_id, int* idof, double* value);
__declspec(dllexport) void mw_get_node_value_(int* node_id, double value[]);
__declspec(dllexport) void mw_get_node_value_with_dof_(int* node_id, int* idof, double* value);

__declspec(dllexport) void mw_set_sv_node_value_(int* node_id, double v_value[], double s_value[]);
__declspec(dllexport) void mw_set_sv_node_value_with_dof_(int* node_id, int* v_dof, double* v_value, int* s_dof, double* s_value);
__declspec(dllexport) void mw_get_sv_node_value_(int* node_id, double v_value[], double s_value[]);
__declspec(dllexport) void mw_get_sv_node_value_with_dof_(int* node_id, int* v_dof, double* v_value, int* s_dof, double* s_value);

// node size, element size
__declspec(dllexport) int mw_get_num_of_node_();// in select mesh_part
__declspec(dllexport) int mw_get_num_of_node_with_mesh_(int* imesh);
__declspec(dllexport) int mw_get_num_of_element_();// in select mesh_part
__declspec(dllexport) int mw_get_num_of_element_with_mesh_(int* imesh);

//----
// node type
//----
__declspec(dllexport) int mw_nodetype_s_();
__declspec(dllexport) int mw_nodetype_v_();
__declspec(dllexport) int mw_nodetype_sv_();
//----
// element type
//----
__declspec(dllexport) int mw_elemtype_hexa_();
__declspec(dllexport) int mw_elemtype_hexa2_();
__declspec(dllexport) int mw_elemtype_tetra_();
__declspec(dllexport) int mw_elemtype_tetra2_();
__declspec(dllexport) int mw_elemtype_prism_();
__declspec(dllexport) int mw_elemtype_prism2_();
__declspec(dllexport) int mw_elemtype_quad_();
__declspec(dllexport) int mw_elemtype_quad2_();
__declspec(dllexport) int mw_elemtype_triangle_();
__declspec(dllexport) int mw_elemtype_triangle2_();
__declspec(dllexport) int mw_elemtype_line_();
__declspec(dllexport) int mw_elemtype_line2_();





//----
// shape function
//----
__declspec(dllexport) int mw_get_num_of_integ_point_(int* shape_type);
__declspec(dllexport) void mw_shape_function_on_pt_(int* shape_type, int* igauss, double N[]);

__declspec(dllexport) void mw_shape_function_hexa81_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa82_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa201_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa202_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa203_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra41_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra101_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra104_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra1015_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism62_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism156_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism159_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism1518_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_quad41_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_quad84_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_quad89_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tri31_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tri63_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_line21_(int* igauss, int* ishape, double* N);
__declspec(dllexport) void mw_shape_function_line32_(int* igauss, int* ishape, double* N);

//----
// shape function deriv (rst coord)
//----
__declspec(dllexport) void mw_dndr_(int* shape_type, double dndr[]);
__declspec(dllexport) void mw_dndr_hexa81_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa82_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa201_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa202_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa203_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra41_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra101_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra104_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra1015_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism62_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism156_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism159_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism1518_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_quad41_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_quad84_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_quad89_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tri31_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tri63_(int* igauss, int* ishape, int* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_line21_(int* igauss, int* ishape, double* dndr);
__declspec(dllexport) void mw_dndr_line32_(int* igauss, int* ishape, double* dndr);
//----
// shape function deriv (xyz coord)
//----
__declspec(dllexport) void mw_dndx_(int* elem_type, int* num_of_integ, int* ielem, double dndx[]);
__declspec(dllexport) void mw_det_jacobian_(int* elem_type, int* num_of_integ, int* igauss, double* det_j);
__declspec(dllexport) void mw_weight_(int* elem_type, int* num_of_integ, int* igauss, double* w);

//----
// shape function type
//----
__declspec(dllexport) int mw_shapetype_hexa81_();
__declspec(dllexport) int mw_shapetype_hexa82_();
__declspec(dllexport) int mw_shapetype_hexa201_();
__declspec(dllexport) int mw_shapetype_hexa202_();
__declspec(dllexport) int mw_shapetype_hexa203_();
__declspec(dllexport) int mw_shapetype_tetra41_();
__declspec(dllexport) int mw_shapetype_tetra101_();
__declspec(dllexport) int mw_shapetype_tetra104_();
__declspec(dllexport) int mw_shapetype_tetra1015_();
__declspec(dllexport) int mw_shapetype_prism62_();
__declspec(dllexport) int mw_shapetype_prism156_();
__declspec(dllexport) int mw_shapetype_prism159_();
__declspec(dllexport) int mw_shapetype_prism1518_();
__declspec(dllexport) int mw_shapetype_quad41_();
__declspec(dllexport) int mw_shapetype_quad84_();
__declspec(dllexport) int mw_shapetype_quad89_();
__declspec(dllexport) int mw_shapetype_tri31_();
__declspec(dllexport) int mw_shapetype_tri63_();
__declspec(dllexport) int mw_shapetype_line21_();
__declspec(dllexport) int mw_shapetype_line32_();


//--
// boundary mesh
//--
__declspec(dllexport) int mw_get_num_of_boundary_bnode_mesh_();
__declspec(dllexport) int mw_get_num_of_boundary_bface_mesh_();
__declspec(dllexport) int mw_get_num_of_boundary_bedge_mesh_();
__declspec(dllexport) int mw_get_num_of_boundary_bvolume_mesh_();
__declspec(dllexport) int mw_get_num_of_bnode_in_bnode_mesh_(int* ibmesh);
__declspec(dllexport) int mw_get_num_of_bnode_in_bface_mesh_(int* ibmesh);
__declspec(dllexport) int mw_get_num_of_bnode_in_bedge_mesh_(int* ibmesh);
__declspec(dllexport) int mw_get_num_of_bnode_in_bvolume_mesh_(int* ibmesh);
__declspec(dllexport) int mw_get_num_of_dof_in_bnode_mesh_(int* ibmesh, int* ibnode);
__declspec(dllexport) int mw_get_num_of_dof_in_bface_mesh_(int* ibmesh);
__declspec(dllexport) int mw_get_num_of_dof_in_bedge_mesh_(int* ibmesh);
__declspec(dllexport) int mw_get_num_of_dof_in_bvolume_mesh_(int* ibmesh);
//--
// value of boundary node
//--
__declspec(dllexport) double mw_get_bnode_value_in_bnode_mesh_(int* ibmesh, int* ibnode, int* idof);
__declspec(dllexport) double mw_get_bnode_value_in_bface_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel);
__declspec(dllexport) double mw_get_bnode_value_in_bedge_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel);
__declspec(dllexport) double mw_get_bnode_value_in_bvolume_mesh_(int* ibmesh, int* ibnode, int* idof, int* mglevel);
__declspec(dllexport) int mw_get_node_id_in_bnode_mesh_(int* ibmesh, int* ibnode);
__declspec(dllexport) int mw_get_node_id_in_bface_mesh_(int* ibmesh, int* ibnode);
__declspec(dllexport) int mw_get_node_id_in_bedge_mesh_(int* ibmesh, int* ibnode);
__declspec(dllexport) int mw_get_node_id_in_bvolume_mesh_(int* ibmesh, int* ibnode);
//--
// value of boundary face, edge, volume
//--
__declspec(dllexport) int mw_get_num_of_bface_(int* ibmesh);
__declspec(dllexport) double mw_get_bface_value_(int* ibmesh, int* ibface, int* idof);
__declspec(dllexport) int mw_get_num_of_bedge_(int* ibmesh);
__declspec(dllexport) double mw_get_bedge_value_(int* ibmesh, int* ibedge, int* idof);
__declspec(dllexport) int mw_get_num_of_bvolume_(int* ibmesh);
__declspec(dllexport) double mw_get_bvolume_value_(int* ibmesh, int* ibvol, int* idof);


//--
// mpi
//--
__declspec(dllexport) int mw_mpi_sum_();// op  ,use allreduce_r argument
__declspec(dllexport) int mw_mpi_max_();// op  ,use allreduce_r argument
__declspec(dllexport) int mw_mpi_min_();// op  ,use allreduce_r argument
__declspec(dllexport) void mw_allreduce_r_(double val[], int* val_size, int* op);
__declspec(dllexport) void mw_send_recv_r2_(double buf[], int* dof_size);// bufの値を送信, 受信値をNodeとbufに代入. bufのサイズ == NumOfCommNode * dof_size
__declspec(dllexport) void mw_send_recv_r_();// 通信Nodeの値を入れ替えて更新

//----
// logger
//----
__declspec(dllexport) void mw_logger_set_mode_(int* mode);
__declspec(dllexport) void mw_logger_set_device_(int* mode, int* device);
__declspec(dllexport) void mw_logger_info_ (int* mode, char* message, int* str_len);
//----
// logger parameter
//----
__declspec(dllexport) int mw_get_error_mode_();
__declspec(dllexport) int mw_get_warn_mode_();
__declspec(dllexport) int mw_get_info_mode_();
__declspec(dllexport) int mw_get_debug_mode_();
__declspec(dllexport) int mw_get_disk_device_();
__declspec(dllexport) int mw_get_display_device_();

#ifdef __cplusplus
}
#endif


#endif	/* API_FORTRAN_HXX */

