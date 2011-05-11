//
//  API_Fortran.h ( Fortran & C )
//
//
//                     2011.4.04
//                     2009.3.26
//                     k.Takeda

#ifndef API_Fortran_AF05492D_2BFF_49fd_AEE3_818B0B426793
#define API_Fortran_AF05492D_2BFF_49fd_AEE3_818B0B426793

#include "TypeDef.h"
#include <cstring>
#include <cmath>
#include <cstdlib>//getenv
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

#define ERROR 0
#define SUCCESS 1


//----
// MW3 construct & destruct
//----
iint mw_initialize_(int* argc, char** argv);//MW3 standard
iint mw_initialize_fstr_(int* argc, char** argv, char* ctrlname);// fstr style
iint mw_finalize_();
//----
// banner
//----
void mw_banner_();

//----
// file i/o API
//----
iint mw_file_read_(char* basename);    // mesh_file read         : General
//iint mw_file_read_bin_(char* basename);// mesh_file(binary) read : General
iint mw_file_read_fstr_();    // mesh_file read              : FrontISTR
//iint mw_file_read_bin_fstr_();// mesh_file read (binary file): FrontISTR
iint mw_file_debug_write_();  // data_check file write
//----
// fstr file_name (hecmw_ctrl)
//----
iint mw_get_fstr_filename_length_mesh_();
iint mw_get_fstr_filename_length_control_();
iint mw_get_fstr_filename_length_result_();
iint mw_get_fstr_filename_length_restart_();
iint mw_get_fstr_filename_length_part_in_();
iint mw_get_fstr_filename_length_part_out_();
iint mw_get_fstr_filename_length_vis_mesh_();
iint mw_get_fstr_filename_length_vis_in_();
iint mw_get_fstr_filename_length_vis_out_();

void mw_get_fstr_filename_mesh_(char name[], iint* len);
void mw_get_fstr_filename_control_(char name[], iint* len);
void mw_get_fstr_filename_result_(char name[], iint* len);
void mw_get_fstr_filename_restart_(char name[], iint* len);
void mw_get_fstr_filename_part_in_(char name[], iint* len);
void mw_get_fstr_filename_part_out_(char name[], iint* len);
void mw_get_fstr_filename_vis_mesh_(char name[], iint* len);
void mw_get_fstr_filename_vis_in_(char name[], iint* len);
void mw_get_fstr_filename_vis_out_(char name[], iint* len);

//----
// result, format = %d:iint*, %f:double*(fixed), %e:double*(scientific), %s:const char*
//----
void mw_rlt_start_(iint* step);
void mw_rlt_start_bin_(iint* step);
iint mw_rlt_print_(iint* width, char* format, ... );
void mw_rlt_end_();
// output *.inp
void mw_print_avs_basis_();// basis variable output
void mw_rec_avs_label_(iint* imesh, char* label, char* unit, iint* ndof);
void mw_rec_avs_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
void mw_print_avs_fem_();  // record variable output

//----
// restart, linear_algebra_equation info
//----
iint mw_file_write_res_(iint* step);// restart file write
iint mw_set_restart_(iint* step);// restart file read, linear_algebra_equation info
iint mw_file_write_res_bin_(iint* step);// restart file write (binary file)
iint mw_set_restart_bin_(iint* step);// restart file read, linear_algebra_equation info (binary file)

//----
// linear solver API
//----
void mw_gene_linear_algebra_(iint* num_of_algebra, iint dof[]);
void mw_select_algebra_(iint* ieq);

iint mw_matrix_add_elem_(iint* imesh, iint* ielem,  double elem_matrix[]);// standard
iint mw_matrix_add_node_(iint* imesh, iint* i_index, iint* j_index, double nodal_matrix[]);

void mw_matrix_clear_(iint* imesh);// matrix 0 clear
void mw_vector_clear_(iint* imesh);// matrix 0 clear

iint mw_matrix_add_elem_24_(iint* imesh, iint* ielem, double elem_matrix[][24]);//Hexa   8Node * 3DOF, Quad 8Node * 3DOF, Quad 4Node * 6DOF
iint mw_matrix_add_elem_60_(iint* imesh, iint* ielem, double elem_matrix[][60]);//Hexa  20Node * 3DOF
iint mw_matrix_add_elem_12_(iint* imesh, iint* ielem, double elem_matirx[][12]);//Tetra  4Node * 3DOF, Quad 4Node * 3DOF, Beam 2Node * 6DOF
iint mw_matrix_add_elem_30_(iint* imesh, iint* ielem, double elem_matirx[][30]);//Tetra 10Node * 3DOF, Tri  6Node * 5DOF
iint mw_matrix_add_elem_18_(iint* imesh, iint* ielem, double elem_matirx[][18]);//Prism  6Node * 3DOF, Tri  6Node * 3DOF, Beam 3Node * 6DOF
iint mw_matirx_add_elem_45_(iint* imesh, iint* ielem, double elem_matirx[][45]);//Prism 15Node * 3DOF
iint mw_matirx_add_elem_20_(iint* imesh, iint* ielem, double elem_matirx[][20]);//Quad   4Node * 5DOF
iint mw_matrix_add_elem_40_(iint* imesh, iint* ielem, double elem_matirx[][40]);//Quad   8Node * 5DOF
iint mw_matrix_add_elem_15_(iint* imesh, iint* ielem, double elem_matirx[][15]);//Tri    3Node * 5DOF, Beam 3Node * 5DOF
iint mw_matirx_add_elem_9_(iint* imesh, iint* ielem, double elem_matirx[][9]);  //Tri    3Node * 3DOF, Beam 3Node * 3DOF
iint mw_matirx_add_elem_48_(iint* imesh, iint* ielem, double elem_matirx[][48]);//Quad   8Node * 6DOF
iint mw_matirx_add_elem_6_(iint* imesh, iint* ielem, double elem_matirx[][6]);  //Beam   2Node * 3DOF
iint mw_matirx_add_elem_10_(iint* imesh, iint* ielem, double elem_matirx[][10]);//Beam   2Node * 5DOF


iint mw_matrix_rhs_set_bc2_(iint* imesh, iint* inode, iint* dof, double* diagval, double* solval);// diagonal=1, off_diag=0 and X
iint mw_matrix_rhs_set_bc_(iint* imesh, iint* inode, iint* dof, double* diagval, double* rhsval); // penalty method, diagonal and rhs_vector
iint mw_rhs_set_bc_(iint* imesh, iint* inode, iint* dof, double* value);
iint mw_rhs_add_bc_(iint* imesh, iint* inode, iint* dof, double* value);

iint mw_solve_(iint* iter_max, double* tolerance, iint* method, iint* pre_condition);

//--
// solution_vector copy,  at select MG-Level && select Equation
//--
void mw_get_solution_vector_(double buf[], iint* imesh);
void mw_get_solution_assy_vector_(double buf[]);

//--
// rhs_vector copy,  at select MG-Level && select Equation
//--
void mw_get_rhs_vector_(double buf[], iint* imesh);
void mw_get_rhs_assy_vector_(double buf[]);

//--
// solution vector value, rhs vector value
//--
double mw_get_solution_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
double mw_get_rhs_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
//--
// solution vector dof, rhs vector dof
//--
iint mw_get_solution_assy_vector_dof_();
iint mw_get_rhs_assy_vector_dof_();

//--
// AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
//--
void mw_mult_vector_(double vX[], double vB[]);




//----
// MG construct (refine)
//----
iint mw_refine_(int* num_of_refine);      // refine == mg_construct
iint mw_mg_construct_(int* num_of_refine);// mg_construct == refine
void mw_finalize_refine_();      // release memory (final proc for mesh construct)
void mw_finalize_mg_construct_();// release memory (final proc for mesh construct) == finalize_refine

//----
// model
//----
// assemble model
iint mw_get_num_of_assemble_model_();
void mw_select_assemble_model_(iint* mglevel);
// mesh part
iint mw_get_num_of_mesh_part_();
void mw_select_mesh_part_with_id_(iint* mesh_id);
void mw_select_mesh_part_(iint* index);
// element
void mw_select_element_with_id_(iint* elem_id);
void mw_select_element_(iint* index);
iint mw_get_element_type_();
iint mw_get_num_of_element_vert_();

void mw_get_element_vert_node_id_(iint v_node_id[]);
iint mw_get_num_of_element_edge_();
void mw_get_element_edge_node_id_(iint v_node_id[]);

// node
void mw_get_node_coord_(iint* node_id, double* x, double* y, double* z);
iint mw_get_dof_(iint* node_id);
iint mw_get_dof_scalar_(iint* node_id);
iint mw_get_dof_vector_(iint* node_id);

////void mw_set_node_value_(iint* node_id, double value[]);
////void mw_set_node_value_with_dof_(iint* node_id, iint* idof, double* value);
////void mw_get_node_value_(iint* node_id, double value[]);
////void mw_get_node_value_with_dof_(iint* node_id, iint* idof, double* value);
////
////void mw_set_sv_node_value_(iint* node_id, double v_value[], double s_value[]);
////void mw_set_sv_node_value_with_dof_(iint* node_id, iint* v_dof, double* v_value, iint* s_dof, double* s_value);
////void mw_get_sv_node_value_(iint* node_id, double v_value[], double s_value[]);
////void mw_get_sv_node_value_with_dof_(iint* node_id, iint* v_dof, double* v_value, iint* s_dof, double* s_value);

// node size, element size
iint mw_get_num_of_node_();// in select mesh_part
iint mw_get_num_of_node_with_mesh_(iint* imesh);
iint mw_get_num_of_element_();// in select mesh_part
iint mw_get_num_of_element_with_mesh_(iint* imesh);

// id && index
iint mw_get_node_id_(iint* index);
iint mw_get_element_id_(iint* index);
iint mw_get_node_index_(iint* id);
iint mw_get_element_index_(iint* id);

//--
// node connectivity,  itemU[]:upper triangle, itemL[]:lower triangle
//--
void mw_construct_node_connect_fem_(iint* node_id);
void mw_get_node_connect_fem_size_(iint* num_itemu, iint* num_iteml);
void mw_get_node_connect_fem_item_(iint itemU[], iint itemL[]);
//--
// node around element_id
//--
iint mw_get_num_of_aggregate_element_(iint* node_id);
iint mw_get_aggregate_element_id_(iint* node_id, iint* index);

//----
// node type
//----
iint mw_nodetype_s_();
iint mw_nodetype_v_();
iint mw_nodetype_sv_();
//----
// element type
//----
iint mw_elemtype_hexa_();
iint mw_elemtype_hexa2_();
iint mw_elemtype_tetra_();
iint mw_elemtype_tetra2_();
iint mw_elemtype_prism_();
iint mw_elemtype_prism2_();
iint mw_elemtype_quad_();
iint mw_elemtype_quad2_();
iint mw_elemtype_triangle_();
iint mw_elemtype_triangle2_();
iint mw_elemtype_line_();
iint mw_elemtype_line2_();
//----
// frontISTR element type
//----
iint mw_fistr_elemtype_hexa_();
iint mw_fistr_elemtype_hexa2_();
iint mw_fistr_elemtype_tetra_();
iint mw_fistr_elemtype_tetra2_();
iint mw_fistr_elemtype_prism_();
iint mw_fistr_elemtype_prism2_();
iint mw_fistr_elemtype_quad_();
iint mw_fistr_elemtype_quad2_();
iint mw_fistr_elemtype_triangle_();
iint mw_fistr_elemtype_triangle2_();
iint mw_fistr_elemtype_line_();
iint mw_fistr_elemtype_line2_();
//----
// frontISTR element type => MW3 element type
//----
iint mw_fistr_elemtype_to_mw3_elemtype_(iint* fistr_elemtype);
//----
// MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
//----
iint mw_mw3_elemtype_to_fistr_elemtype_(iint* mw3_elemtype);


//----
// shape function
//----
iint mw_get_num_of_integ_point_(iint* shape_type);
void mw_shape_function_on_pt_(iint* shape_type, iint* igauss, double N[]);

void mw_shape_function_hexa81_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_hexa82_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_hexa201_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_hexa202_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_hexa203_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_tetra41_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_tetra101_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_tetra104_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_tetra1015_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_prism62_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_prism156_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_prism159_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_prism1518_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_quad41_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_quad84_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_quad89_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_tri31_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_tri63_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_line21_(iint* igauss, iint* ishape, double* N);
void mw_shape_function_line32_(iint* igauss, iint* ishape, double* N);

//----
// shape function deriv (rst coord)
//----
void mw_dndr_(iint* shape_type, double dndr[]);
void mw_dndr_hexa81_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_hexa82_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_hexa201_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_hexa202_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_hexa203_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_tetra41_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_tetra101_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_tetra104_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_tetra1015_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_prism62_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_prism156_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_prism159_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_prism1518_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_quad41_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_quad84_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_quad89_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_tri31_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_tri63_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
void mw_dndr_line21_(iint* igauss, iint* ishape, double* dndr);
void mw_dndr_line32_(iint* igauss, iint* ishape, double* dndr);
//----
// shape function deriv (xyz coord)
//----
void mw_dndx_(iint* elem_type, iint* num_of_integ, iint* ielem, double dndx[]);
void mw_det_jacobian_(iint* elem_type, iint* num_of_integ, iint* igauss, double* det_j);
void mw_weight_(iint* elem_type, iint* num_of_integ, iint* igauss, double* w);

//----
// shape function type
//----
iint mw_shapetype_hexa81_();
iint mw_shapetype_hexa82_();
iint mw_shapetype_hexa201_();
iint mw_shapetype_hexa202_();
iint mw_shapetype_hexa203_();
iint mw_shapetype_tetra41_();
iint mw_shapetype_tetra101_();
iint mw_shapetype_tetra104_();
iint mw_shapetype_tetra1015_();
iint mw_shapetype_prism62_();
iint mw_shapetype_prism156_();
iint mw_shapetype_prism159_();
iint mw_shapetype_prism1518_();
iint mw_shapetype_quad41_();
iint mw_shapetype_quad84_();
iint mw_shapetype_quad89_();
iint mw_shapetype_tri31_();
iint mw_shapetype_tri63_();
iint mw_shapetype_line21_();
iint mw_shapetype_line32_();


//--
// boundary mesh
//--
// number of boudary_mesh
iint mw_get_num_of_boundary_bnode_mesh_();
iint mw_get_num_of_boundary_bface_mesh_();
iint mw_get_num_of_boundary_bedge_mesh_();
iint mw_get_num_of_boundary_bvolume_mesh_();
// BND type for each boundary_mesh { Neumann || Dirichlet }
iint mw_get_bnd_type_bnode_mesh_(iint* ibmesh);
iint mw_get_bnd_type_bface_mesh_(iint* ibmesh);
iint mw_get_bnd_type_bedge_mesh_(iint* ibmesh);
iint mw_get_bnd_type_bvolume_mesh_(iint* ibmesh);
// BND type number
iint mw_get_neumann_type_();
iint mw_get_dirichlet_type_();
// number of bnode for each boundary_mesh
iint mw_get_num_of_bnode_in_bnode_mesh_(iint* ibmesh);
iint mw_get_num_of_bnode_in_bface_mesh_(iint* ibmesh);
iint mw_get_num_of_bnode_in_bedge_mesh_(iint* ibmesh);
iint mw_get_num_of_bnode_in_bvolume_mesh_(iint* ibmesh);
// number of DOF for each boundary_mesh
iint mw_get_num_of_dof_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
iint mw_get_num_of_dof_in_bface_mesh_(iint* ibmesh);
iint mw_get_num_of_dof_in_bedge_mesh_(iint* ibmesh);
iint mw_get_num_of_dof_in_bvolume_mesh_(iint* ibmesh);
// DOF number for each boundary_mesh ( DOF_index => DOF Number )
iint mw_get_dof_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* idof);
iint mw_get_dof_bface_mesh_(iint* ibmesh, iint* idof);
iint mw_get_dof_bedge_mesh_(iint* ibmesh, iint* idof);
iint mw_get_dof_bvolume_mesh_(iint* ibmesh, iint* idof);
//--
// value of boundary node
//--
double mw_get_bnode_value_in_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* dof);
double mw_get_bnode_value_in_bface_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
double mw_get_bnode_value_in_bedge_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
double mw_get_bnode_value_in_bvolume_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
iint mw_get_node_id_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
iint mw_get_node_id_in_bface_mesh_(iint* ibmesh, iint* ibnode);
iint mw_get_node_id_in_bedge_mesh_(iint* ibmesh, iint* ibnode);
iint mw_get_node_id_in_bvolume_mesh_(iint* ibmesh, iint* ibnode);
//--
// value of boundary face, edge, volume
//--
iint mw_get_num_of_bface_(iint* ibmesh);
double mw_get_bface_value_(iint* ibmesh, iint* ibface, iint* dof);
iint mw_get_num_of_bedge_(iint* ibmesh);
double mw_get_bedge_value_(iint* ibmesh, iint* ibedge, iint* dof);
iint mw_get_num_of_bvolume_(iint* ibmesh);
double mw_get_bvolume_value_(iint* ibmesh, iint* ibvol, iint* dof);
//--
// node_id, face, edge, volume
//--
iint mw_get_num_of_node_bface_(iint* ibmesh, iint* ibface);
iint mw_get_node_id_bface_(iint* ibmesh, iint* ibface, iint* ibnode);
iint mw_get_num_of_node_bedge_(iint* ibmesh, iint* ibedge);
iint mw_get_node_id_bedge_(iint* ibmesh, iint* ibedge, iint* ibnode);
iint mw_get_num_of_node_bvolume_(iint* ibmesh, iint* ibvol);
iint mw_get_node_id_bvolume_(iint* ibmesh, iint* ibvol, iint* ibnode);
//--
// boundary_mesh name
//--
iint mw_get_bnode_mesh_namelength_(iint* ibmesh);
void mw_get_bnode_mesh_name_(iint* ibmesh, char* name, iint* name_len);
iint mw_get_bface_mesh_namelength_(iint* ibmesh);
void mw_get_bface_mesh_name_(iint* ibmesh, char* name, iint* name_len);
iint mw_get_bvolume_mesh_namelength_(iint* ibmesh);
void mw_get_bvolume_mesh_name_(iint* ibmesh, char* name, iint* name_len);
iint mw_get_bedge_mesh_namelength_(iint* ibmesh);
void mw_get_bedge_mesh_name_(iint* ibmesh, char* name, iint* name_len);
//--
// entity_id of boundary_mesh (for FrontISTR)
//--
iint mw_get_edge_id_bedge_(iint* ibmesh, iint* ibedge);
iint mw_get_elem_id_bedge_(iint* ibmesh, iint* ibedge);
iint mw_get_face_id_bface_(iint* ibmesh, iint* ibface);
iint mw_get_elem_id_bface_(iint* ibmesh, iint* ibface);
iint mw_get_elem_id_bvolume_(iint* ibmesh, iint* ibvol);




//--
// mpi
//--
int mw_mpi_int_();   // MPI_INT
int mw_mpi_double_();// MPI_DOUBLE
int mw_mpi_comm_();  // MPI_COMM_WORLD

int mw_mpi_sum_();// op  ,use allreduce_r argument
int mw_mpi_max_();// op  ,use allreduce_r argument
int mw_mpi_min_();// op  ,use allreduce_r argument

int mw_get_rank_();
int mw_get_num_of_process_();

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
int mw_bcast_i_(int buf[], int* cnt, int* root);
int mw_bcast_r_(double buf[], int* cnt, int* root);
int mw_bcast_s_(char buf[], int* cnt, int* root);

// 以下の4メソッドは、ペア
int mw_get_num_of_neibpe_(int* imesh);//Meshパーツが通信する相手の数
int mw_get_transrank_(int* imesh, int* ipe);//通信Mesh毎のランク番号
void mw_send_recv_r_(double buf[], int* num_of_node, int* dof_size, int* trans_rank);//bufの値を送信-受信
void mw_send_recv_i_(int buf[], int* num_of_node, int* dof_size, int* trans_rank);
////void mw_send_recv_r2_(double buf[], iint* dof_size);// bufの値を送信, 受信値をNodeとbufに代入. bufのサイズ == NumOfCommNode * dof_size
////void mw_send_recv_r_();// 通信Nodeの値を入れ替えて更新


//--
// CommMesh2 (通信Mesh) { select された AssyModel,Meshを対象 }
// CommNode (通信Node) :  Visualizer 用途
//--
iint mw_get_num_of_comm_mesh_();
iint mw_get_num_of_comm_node_(iint* icmesh);
iint mw_get_node_id_comm_node_(iint* icmesh, iint* icnode);//MeshのNodeID


//--
// Element_Group { select AssyModel, select Mesh }
//--
iint mw_get_num_of_elementgroup_();
iint mw_get_num_of_element_id_(iint* igrp);
iint mw_get_element_id_with_elementgroup_(iint* igrp, iint* index);
iint mw_get_elementgroup_name_length_(iint* igrp);
void mw_get_elementgroup_name_(iint* igrp, char* name, iint* name_len);

//----
// logger
//----
void mw_logger_set_mode_(iint* mode);
void mw_logger_set_device_(iint* mode, iint* device);
void mw_logger_info_ (iint* mode, char* message, iint* str_len);
//----
// logger parameter
//----
iint mw_get_error_mode_();
iint mw_get_warn_mode_();
iint mw_get_info_mode_();
iint mw_get_debug_mode_();
iint mw_get_disk_device_();
iint mw_get_display_device_();

#ifdef __cplusplus
}
#endif

#endif // API_Fortran.h

