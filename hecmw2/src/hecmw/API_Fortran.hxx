/* 
 * File:   API_Fortran.hxx
 * Author: ktakeda
 *
 * Created on 2010/10/04, 19:18
 */

#ifndef API_FORTRAN_HXX_VISUAL_CPP
#define	API_FORTRAN_HXX_VISUAL_CPP

#include "TypeDef.h"
#include <cstring>
#include <cmath>
#include <cstdlib>

#ifdef __cplusplus
extern "C" {
#endif

#define ERROR 0
#define SUCCESS 1


//----
// MW3 construct & destruct
//----
__declspec(dllexport) iint mw_initialize_(int* argc, char** argv);//MW3 standard
__declspec(dllexport) iint mw_initialize_fstr_(int* argc, char** argv, char* ctrlname);//fstr style
__declspec(dllexport) iint mw_finalize_();

//----
// banner
//----
__declspec(dllexport) void mw_banner_();

//----
// file i/o API
//----
__declspec(dllexport) iint mw_file_read_(char* basename);    // mesh_file read         : General
//__declspec(dllexport) iint mw_file_read_bin_(char* basename);// mesh_file(binary) read : General
__declspec(dllexport) iint mw_file_read_fstr_();    // mesh_file read              : FrontISTR
//__declspec(dllexport) iint mw_file_read_bin_fstr_();// mesh_file read (binary file): FrontISTR
__declspec(dllexport) iint mw_file_debug_write_();  // data_check file write
//----
// fstr file_name (hecmw_ctrl)
//----
__declspec(dllexport) iint mw_get_fstr_filename_length_mesh_();
__declspec(dllexport) iint mw_get_fstr_filename_length_control_();
__declspec(dllexport) iint mw_get_fstr_filename_length_result_();
__declspec(dllexport) iint mw_get_fstr_filename_length_restart_();
__declspec(dllexport) iint mw_get_fstr_filename_length_part_in_();
__declspec(dllexport) iint mw_get_fstr_filename_length_part_out_();
__declspec(dllexport) iint mw_get_fstr_filename_length_vis_mesh_();
__declspec(dllexport) iint mw_get_fstr_filename_length_vis_in_();
__declspec(dllexport) iint mw_get_fstr_filename_length_vis_out_();

__declspec(dllexport) void mw_get_fstr_filename_mesh_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_control_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_result_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_restart_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_part_in_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_part_out_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_vis_mesh_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_vis_in_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_vis_out_(char name[], iint* len);

//----
// result, format = %d:iint*, %f:double*(fixed), %e:double*(scientific), %s:const char*
//----
__declspec(dllexport) void mw_rlt_start_(iint* step);
__declspec(dllexport) void mw_rlt_start_bin_(iint* step);
__declspec(dllexport) iint mw_rlt_print_(iint* width, char* format, ... );
__declspec(dllexport) void mw_rlt_end_();
// output *.inp
__declspec(dllexport) void mw_print_avs_basis_();// basis variable output
__declspec(dllexport) void mw_rec_avs_label_(iint* imesh, char* label, char* unit, iint* ndof);
__declspec(dllexport) void mw_rec_avs_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
__declspec(dllexport) void mw_print_avs_fem_();  // record variable output

//----
// restart, linear_algebra_equation info
//----
__declspec(dllexport) iint mw_file_write_res_(iint* step);// restart file write
__declspec(dllexport) iint mw_set_restart_(iint* step);// restart file read, linear_algebra_equation info
__declspec(dllexport) iint mw_file_write_res_bin_(iint* step);// restart file write (binary file)
__declspec(dllexport) iint mw_set_restart_bin_(iint* step);// restart file read, linear_algebra_equation info (binary file)

//----
// linear solver API
//----
__declspec(dllexport) void mw_gene_linear_algebra_(iint* num_of_algebra, iint dof[]);
__declspec(dllexport) void mw_select_algebra_(iint* ieq);

__declspec(dllexport) iint mw_matrix_add_elem_(iint* imesh,  iint* ielem,  double elem_matrix[]);// standard
__declspec(dllexport) iint mw_matrix_add_node_(iint* imesh, iint* i_index, iint* j_index, double nodal_matrix[]);

__declspec(dllexport) void mw_matrix_clear_(iint* imesh);// matrix 0 clear
__declspec(dllexport) void mw_vector_clear_(iint* imesh);// matrix 0 clear

__declspec(dllexport) iint mw_matrix_add_elem_24_(iint* imesh, iint* ielem, double elem_matrix[][24]);//Hexa   8Node * 3DOF, Quad 8Node * 3DOF, Quad 4Node * 6DOF
__declspec(dllexport) iint mw_matrix_add_elem_60_(iint* imesh, iint* ielem, double elem_matrix[][60]);//Hexa  20Node * 3DOF
__declspec(dllexport) iint mw_matrix_add_elem_12_(iint* imesh, iint* ielem, double elem_matirx[][12]);//Tetra  4Node * 3DOF, Quad 4Node * 3DOF, Beam 2Node * 6DOF
__declspec(dllexport) iint mw_matrix_add_elem_30_(iint* imesh, iint* ielem, double elem_matirx[][30]);//Tetra 10Node * 3DOF, Tri  6Node * 5DOF
__declspec(dllexport) iint mw_matrix_add_elem_18_(iint* imesh, iint* ielem, double elem_matirx[][18]);//Prism  6Node * 3DOF, Tri  6Node * 3DOF, Beam 3Node * 6DOF
__declspec(dllexport) iint mw_matirx_add_elem_45_(iint* imesh, iint* ielem, double elem_matirx[][45]);//Prism 15Node * 3DOF
__declspec(dllexport) iint mw_matirx_add_elem_20_(iint* imesh, iint* ielem, double elem_matirx[][20]);//Quad   4Node * 5DOF
__declspec(dllexport) iint mw_matrix_add_elem_40_(iint* imesh, iint* ielem, double elem_matirx[][40]);//Quad   8Node * 5DOF
__declspec(dllexport) iint mw_matrix_add_elem_15_(iint* imesh, iint* ielem, double elem_matirx[][15]);//Tri    3Node * 5DOF, Beam 3Node * 5DOF
__declspec(dllexport) iint mw_matirx_add_elem_9_(iint* imesh, iint* ielem, double elem_matirx[][9]);  //Tri    3Node * 3DOF, Beam 3Node * 3DOF
__declspec(dllexport) iint mw_matirx_add_elem_48_(iint* imesh, iint* ielem, double elem_matirx[][48]);//Quad   8Node * 6DOF
__declspec(dllexport) iint mw_matirx_add_elem_6_(iint* imesh, iint* ielem, double elem_matirx[][6]);  //Beam   2Node * 3DOF
__declspec(dllexport) iint mw_matirx_add_elem_10_(iint* imesh, iint* ielem, double elem_matirx[][10]);//Beam   2Node * 5DOF


__declspec(dllexport) iint mw_matrix_rhs_set_bc2_(iint* imesh, iint* inode, iint* dof, double* diagval, double* solval);// diagonal=1, off_diag=0 and X
__declspec(dllexport) iint mw_matrix_rhs_set_bc_(iint* imesh, iint* inode, iint* dof, double* diagval, double* rhsval); // penalty method, diagonal and rhs_vector
__declspec(dllexport) iint mw_rhs_set_bc_(iint* imesh, iint* inode, iint* dof, double* value);// rhs_vector
__declspec(dllexport) iint mw_rhs_add_bc_(iint* imesh, iint* inode, iint* dof, double* value);

__declspec(dllexport) iint mw_solve_(iint* iter_max, double* tolerance, iint* method, iint* pre_condition);

//--
// solution_vector copy,  at select MG-Level && select Equation
//--
__declspec(dllexport) void mw_get_solution_vector_(double buf[], iint* imesh);
__declspec(dllexport) void mw_get_solution_assy_vector_(double buf[]);
//--
// rhs_vector copy,  at select MG-Level && select Equation
//--
__declspec(dllexport) void mw_get_rhs_vector_(double buf[], iint* imesh);
__declspec(dllexport) void mw_get_rhs_assy_vector_(double buf[]);

//--
// solution vector value, rhs vector value
//--
__declspec(dllexport) double mw_get_solution_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
__declspec(dllexport) double mw_get_rhs_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
//--
// solution vector dof, rhs vector dof
//--
__declspec(dllexport) iint mw_get_solution_assy_vector_dof_();
__declspec(dllexport) iint mw_get_rhs_assy_vector_dof_();

//--
// AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
//--
__declspec(dllexport) void mw_mult_vector_(double vX[], double vB[]);




//----
// MG construct (refine)
//----
__declspec(dllexport) iint mw_refine_(int* num_of_refine);      // refine == mg_construct
__declspec(dllexport) iint mw_mg_construct_(int* num_of_refine);// mg_construct == refine
__declspec(dllexport) void mw_finalize_refine_();      // release memory (final proc for mesh construct)
__declspec(dllexport) void mw_finalize_mg_construct_();// release memory (final proc for mesh construct) == finalize_refine

//----
// model
//----
// assemble model
__declspec(dllexport) iint mw_get_num_of_assemble_model_();
__declspec(dllexport) void mw_select_assemble_model_(iint* mglevel);
// mesh part
__declspec(dllexport) iint mw_get_num_of_mesh_part_();
__declspec(dllexport) void mw_select_mesh_part_with_id_(iint* mesh_id);
__declspec(dllexport) void mw_select_mesh_part_(iint* index);
// element
__declspec(dllexport) void mw_select_element_with_id_(iint* elem_id);
__declspec(dllexport) void mw_select_element_(iint* index);
__declspec(dllexport) iint mw_get_element_type_();
__declspec(dllexport) iint mw_get_num_of_element_vert_();

__declspec(dllexport) void mw_get_element_vert_node_id_(iint v_node_id[]);
__declspec(dllexport) iint mw_get_num_of_element_edge_();
__declspec(dllexport) void mw_get_element_edge_node_id_(iint v_node_id[]);

// node
__declspec(dllexport) void mw_get_node_coord_(iint* node_id, double* x, double* y, double* z);
__declspec(dllexport) iint mw_get_dof_(iint* node_id);
__declspec(dllexport) iint mw_get_dof_scalar_(iint* node_id);
__declspec(dllexport) iint mw_get_dof_vector_(iint* node_id);

////__declspec(dllexport) void mw_set_node_value_(iint* node_id, double value[]);
////__declspec(dllexport) void mw_set_node_value_with_dof_(iint* node_id, iint* idof, double* value);
////__declspec(dllexport) void mw_get_node_value_(iint* node_id, double value[]);
////__declspec(dllexport) void mw_get_node_value_with_dof_(iint* node_id, iint* idof, double* value);

////__declspec(dllexport) void mw_set_sv_node_value_(iint* node_id, double v_value[], double s_value[]);
////__declspec(dllexport) void mw_set_sv_node_value_with_dof_(iint* node_id, iint* v_dof, double* v_value, iint* s_dof, double* s_value);
////__declspec(dllexport) void mw_get_sv_node_value_(iint* node_id, double v_value[], double s_value[]);
////__declspec(dllexport) void mw_get_sv_node_value_with_dof_(iint* node_id, iint* v_dof, double* v_value, iint* s_dof, double* s_value);

// node size, element size
__declspec(dllexport) iint mw_get_num_of_node_();// in select mesh_part
__declspec(dllexport) iint mw_get_num_of_node_with_mesh_(iint* imesh);
__declspec(dllexport) iint mw_get_num_of_element_();// in select mesh_part
__declspec(dllexport) iint mw_get_num_of_element_with_mesh_(iint* imesh);

// id && index
__declspec(dllexport) iint mw_get_node_id_(iint* index);
__declspec(dllexport) iint mw_get_element_id_(iint* index);
__declspec(dllexport) iint mw_get_node_index_(iint* id);
__declspec(dllexport) iint mw_get_element_index_(iint* id);

//--
// node connectivity,  itemU[]:upper triangle, itemL[]:lower triangle
//--
__declspec(dllexport) void mw_construct_node_connect_fem_(iint* node_id);
__declspec(dllexport) void mw_get_node_connect_fem_size_(iint* num_itemu, iint* num_iteml);
__declspec(dllexport) void mw_get_node_connect_fem_item_(iint itemU[], iint itemL[]);
//--
// node around element_id
//--
__declspec(dllexport) iint mw_get_num_of_aggregate_element_(iint* node_id);
__declspec(dllexport) iint mw_get_aggregate_element_id_(iint* node_id, iint* index);


//----
// node type
//----
__declspec(dllexport) iint mw_nodetype_s_();
__declspec(dllexport) iint mw_nodetype_v_();
__declspec(dllexport) iint mw_nodetype_sv_();
//----
// element type
//----
__declspec(dllexport) iint mw_elemtype_hexa_();
__declspec(dllexport) iint mw_elemtype_hexa2_();
__declspec(dllexport) iint mw_elemtype_tetra_();
__declspec(dllexport) iint mw_elemtype_tetra2_();
__declspec(dllexport) iint mw_elemtype_prism_();
__declspec(dllexport) iint mw_elemtype_prism2_();
__declspec(dllexport) iint mw_elemtype_quad_();
__declspec(dllexport) iint mw_elemtype_quad2_();
__declspec(dllexport) iint mw_elemtype_triangle_();
__declspec(dllexport) iint mw_elemtype_triangle2_();
__declspec(dllexport) iint mw_elemtype_line_();
__declspec(dllexport) iint mw_elemtype_line2_();
//----
// frontISTR element type
//----
__declspec(dllexport) iint mw_fistr_elemtype_hexa_();
__declspec(dllexport) iint mw_fistr_elemtype_hexa2_();
__declspec(dllexport) iint mw_fistr_elemtype_tetra_();
__declspec(dllexport) iint mw_fistr_elemtype_tetra2_();
__declspec(dllexport) iint mw_fistr_elemtype_prism_();
__declspec(dllexport) iint mw_fistr_elemtype_prism2_();
__declspec(dllexport) iint mw_fistr_elemtype_quad_();
__declspec(dllexport) iint mw_fistr_elemtype_quad2_();
__declspec(dllexport) iint mw_fistr_elemtype_triangle_();
__declspec(dllexport) iint mw_fistr_elemtype_triangle2_();
__declspec(dllexport) iint mw_fistr_elemtype_line_();
__declspec(dllexport) iint mw_fistr_elemtype_line2_();
//----
// frontISTR element type => MW3 element type
//----
__declspec(dllexport) iint mw_fistr_elemtype_to_mw3_elemtype_(iint* fistr_elemtype);
//----
// MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
//----
__declspec(dllexport) iint mw_mw3_elemtype_to_fistr_elemtype_(iint* mw3_elemtype);


//----
// shape function
//----
__declspec(dllexport) iint mw_get_num_of_integ_point_(iint* shape_type);
__declspec(dllexport) void mw_shape_function_on_pt_(iint* shape_type, iint* igauss, double N[]);

__declspec(dllexport) void mw_shape_function_hexa81_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa82_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa201_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa202_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_hexa203_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra41_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra101_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra104_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tetra1015_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism62_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism156_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism159_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_prism1518_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_quad41_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_quad84_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_quad89_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tri31_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_tri63_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_line21_(iint* igauss, iint* ishape, double* N);
__declspec(dllexport) void mw_shape_function_line32_(iint* igauss, iint* ishape, double* N);

//----
// shape function deriv (rst coord)
//----
__declspec(dllexport) void mw_dndr_(iint* shape_type, double dndr[]);
__declspec(dllexport) void mw_dndr_hexa81_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa82_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa201_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa202_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_hexa203_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra41_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra101_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra104_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tetra1015_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism62_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism156_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism159_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_prism1518_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_quad41_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_quad84_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_quad89_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tri31_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_tri63_(iint* igauss, iint* ishape, iint* iaxis, double* dndr);
__declspec(dllexport) void mw_dndr_line21_(iint* igauss, iint* ishape, double* dndr);
__declspec(dllexport) void mw_dndr_line32_(iint* igauss, iint* ishape, double* dndr);
//----
// shape function deriv (xyz coord)
//----
__declspec(dllexport) void mw_dndx_(iint* elem_type, iint* num_of_integ, iint* ielem, double dndx[]);
__declspec(dllexport) void mw_det_jacobian_(iint* elem_type, iint* num_of_integ, iint* igauss, double* det_j);
__declspec(dllexport) void mw_weight_(iint* elem_type, iint* num_of_integ, iint* igauss, double* w);

//----
// shape function type
//----
__declspec(dllexport) iint mw_shapetype_hexa81_();
__declspec(dllexport) iint mw_shapetype_hexa82_();
__declspec(dllexport) iint mw_shapetype_hexa201_();
__declspec(dllexport) iint mw_shapetype_hexa202_();
__declspec(dllexport) iint mw_shapetype_hexa203_();
__declspec(dllexport) iint mw_shapetype_tetra41_();
__declspec(dllexport) iint mw_shapetype_tetra101_();
__declspec(dllexport) iint mw_shapetype_tetra104_();
__declspec(dllexport) iint mw_shapetype_tetra1015_();
__declspec(dllexport) iint mw_shapetype_prism62_();
__declspec(dllexport) iint mw_shapetype_prism156_();
__declspec(dllexport) iint mw_shapetype_prism159_();
__declspec(dllexport) iint mw_shapetype_prism1518_();
__declspec(dllexport) iint mw_shapetype_quad41_();
__declspec(dllexport) iint mw_shapetype_quad84_();
__declspec(dllexport) iint mw_shapetype_quad89_();
__declspec(dllexport) iint mw_shapetype_tri31_();
__declspec(dllexport) iint mw_shapetype_tri63_();
__declspec(dllexport) iint mw_shapetype_line21_();
__declspec(dllexport) iint mw_shapetype_line32_();


//--
// boundary mesh
//--
// number of boudary_mesh
__declspec(dllexport) iint mw_get_num_of_boundary_bnode_mesh_();
__declspec(dllexport) iint mw_get_num_of_boundary_bface_mesh_();
__declspec(dllexport) iint mw_get_num_of_boundary_bedge_mesh_();
__declspec(dllexport) iint mw_get_num_of_boundary_bvolume_mesh_();
// BND type for each boundary_mesh { Neumann || Dirichlet }
__declspec(dllexport) iint mw_get_bnd_type_bnode_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_bnd_type_bface_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_bnd_type_bedge_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_bnd_type_bvolume_mesh_(iint* ibmesh);
// BND type number
__declspec(dllexport) iint mw_get_neumann_type_();
__declspec(dllexport) iint mw_get_dirichlet_type_();
// number of bnode for each boundary_mesh
__declspec(dllexport) iint mw_get_num_of_bnode_in_bnode_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_bnode_in_bface_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_bnode_in_bedge_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_bnode_in_bvolume_mesh_(iint* ibmesh);
// number of DOF for each boundary_mesh
__declspec(dllexport) iint mw_get_num_of_dof_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_dof_in_bface_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_dof_in_bedge_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_dof_in_bvolume_mesh_(iint* ibmesh);
// DOF number for each boundary_mesh ( DOF_index => DOF Number )
__declspec(dllexport) iint mw_get_dof_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* idof);
__declspec(dllexport) iint mw_get_dof_bface_mesh_(iint* ibmesh, iint* idof);
__declspec(dllexport) iint mw_get_dof_bedge_mesh_(iint* ibmesh, iint* idof);
__declspec(dllexport) iint mw_get_dof_bvolume_mesh_(iint* ibmesh, iint* idof);
//--
// value of boundary node
//--
__declspec(dllexport) double mw_get_bnode_value_in_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* dof);
__declspec(dllexport) double mw_get_bnode_value_in_bface_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
__declspec(dllexport) double mw_get_bnode_value_in_bedge_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
__declspec(dllexport) double mw_get_bnode_value_in_bvolume_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
__declspec(dllexport) iint mw_get_node_id_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_node_id_in_bface_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_node_id_in_bedge_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_node_id_in_bvolume_mesh_(iint* ibmesh, iint* ibnode);
//--
// value of boundary face, edge, volume
//--
__declspec(dllexport) iint mw_get_num_of_bface_(iint* ibmesh);
__declspec(dllexport) double mw_get_bface_value_(iint* ibmesh, iint* ibface, iint* dof);
__declspec(dllexport) iint mw_get_num_of_bedge_(iint* ibmesh);
__declspec(dllexport) double mw_get_bedge_value_(iint* ibmesh, iint* ibedge, iint* dof);
__declspec(dllexport) iint mw_get_num_of_bvolume_(iint* ibmesh);
__declspec(dllexport) double mw_get_bvolume_value_(iint* ibmesh, iint* ibvol, iint* dof);
//--
// node_id, face, edge, volume
//--
__declspec(dllexport) iint mw_get_num_of_node_bface_(iint* ibmesh, iint* ibface);
__declspec(dllexport) iint mw_get_node_id_bface_(iint* ibmesh, iint* ibface, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_node_bedge_(iint* ibmesh, iint* ibedge);
__declspec(dllexport) iint mw_get_node_id_bedge_(iint* ibmesh, iint* ibedge, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_node_bvolume_(iint* ibmesh, iint* ibvol);
__declspec(dllexport) iint mw_get_node_id_bvolume_(iint* ibmesh, iint* ibvol, iint* ibnode);
//--
// boundary_mesh name
//--
__declspec(dllexport) iint mw_get_bnode_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bnode_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_bface_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bface_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_bvolume_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bvolume_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_bedge_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bedge_mesh_name_(iint* ibmesh, char* name, iint* name_len);
//--
// entity_id of boundary_mesh (for FrontISTR)
//--
__declspec(dllexport) iint mw_get_edge_id_bedge_(iint*ibmesh, iint* ibedge);
__declspec(dllexport) iint mw_get_elem_id_bedge_(iint* ibmesh, iint* ibedge);
__declspec(dllexport) iint mw_get_face_id_bface_(iint* ibmesh, iint* ibface);
__declspec(dllexport) iint mw_get_elem_id_bface_(iint* ibmesh, iint* ibface);
__declspec(dllexport) iint mw_get_elem_id_bvolume_(iint* ibmesh, iint* ibvol);





//--
// mpi
//--
__declspec(dllexport) int mw_mpi_int_();   // MPI_INT
__declspec(dllexport) int mw_mpi_double_();// MPI_DOUBLE
__declspec(dllexport) int mw_mpi_comm_();  // MPI_COMM_WORLD

__declspec(dllexport) int mw_mpi_sum_();// op  ,use allreduce_r argument
__declspec(dllexport) int mw_mpi_max_();// op  ,use allreduce_r argument
__declspec(dllexport) int mw_mpi_min_();// op  ,use allreduce_r argument

__declspec(dllexport) int mw_get_rank_();
__declspec(dllexport) int mw_get_num_of_process_();

__declspec(dllexport) void mw_allreduce_r_(double val[], int* val_size, int* op);
__declspec(dllexport) void mw_allreduce_i_(int val[], int* val_size, int* op);
__declspec(dllexport) int mw_barrier_();
__declspec(dllexport) int mw_abort_(int* error);
__declspec(dllexport) int mw_allgather_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt);
__declspec(dllexport) int mw_allgather_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt);
__declspec(dllexport) int mw_gather_r_(double sendbuf[], int* sendcnt, double recvbuf[], iint* recvcnt, iint* root);
__declspec(dllexport) int mw_gather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], int* recvcnt, int* root);
__declspec(dllexport) int mw_scatter_r_(double sendbuf[], int* sendcnt, double recvbuf[], int* recvcnt, int* root);
__declspec(dllexport) int mw_scatter_i_(int sendbuf[], int* sendcnt, int recvbuf[], int* recvcnt, int* root);
__declspec(dllexport) int mw_bcast_i_(int buf[], int* cnt, int* root);
__declspec(dllexport) int mw_bcast_r_(double buf[], int* cnt, int* root);
__declspec(dllexport) int mw_bcast_s_(char buf[], int* cnt, int* root);

// 以下の4メソッドは、ペア
__declspec(dllexport) int mw_get_num_of_neibpe_(int* imesh);//Meshパーツが通信する相手の数
__declspec(dllexport) int mw_get_transrank_(int* imesh, int* ipe);//通信Mesh毎のランク番号
__declspec(dllexport) void mw_send_recv_r_(double buf[], int* num_of_node, int* dof_size, int* trans_rank);//bufの値を送信-受信
__declspec(dllexport) void mw_send_recv_i_(int buf[], int* num_of_node, int* dof_size, int* trans_rank);
////__declspec(dllexport) void mw_send_recv_r2_(double buf[], iint* dof_size);// bufの値を送信, 受信値をNodeとbufに代入. bufのサイズ == NumOfCommNode * dof_size
////__declspec(dllexport) void mw_send_recv_r_();// 通信Nodeの値を入れ替えて更新


//--
// CommMesh2 (通信Mesh) { select された AssyModel,Meshを対象 }
// CommNode (通信Node) :  Visualizer 用途
//--
__declspec(dllexport) iint mw_get_num_of_comm_mesh_();
__declspec(dllexport) iint mw_get_num_of_comm_node_(iint* icmesh);
__declspec(dllexport) iint mw_get_node_id_comm_node_(iint* icmesh, iint* icnode);//MeshのNodeID


//--
// Element_Group { select AssyModel, select Mesh }
//--
__declspec(dllexport) iint mw_get_num_of_elementgroup_();
__declspec(dllexport) iint mw_get_num_of_element_id_(iint* igrp);
__declspec(dllexport) iint mw_get_element_id_with_elementgroup_(iint* igrp, iint* index);
__declspec(dllexport) iint mw_get_elementgroup_name_length_(iint* igrp);
__declspec(dllexport) void mw_get_elementgroup_name_(iint* igrp, char* name, iint* name_len);


//----
// logger
//----
__declspec(dllexport) void mw_logger_set_mode_(iint* mode);
__declspec(dllexport) void mw_logger_set_device_(iint* mode, iint* device);
__declspec(dllexport) void mw_logger_info_ (iint* mode, char* message, iint* str_len);
//----
// logger parameter
//----
__declspec(dllexport) iint mw_get_error_mode_();
__declspec(dllexport) iint mw_get_warn_mode_();
__declspec(dllexport) iint mw_get_info_mode_();
__declspec(dllexport) iint mw_get_debug_mode_();
__declspec(dllexport) iint mw_get_disk_device_();
__declspec(dllexport) iint mw_get_display_device_();

#ifdef __cplusplus
}
#endif


#endif	/* API_FORTRAN_HXX */

