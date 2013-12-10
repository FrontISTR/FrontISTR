/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.2 beta
|
|   ./src/API_Fortran.hxx
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
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

__declspec(dllexport) iint mw_initialize_(int* argc, char** argv);
__declspec(dllexport) iint mw_initialize_fstr_(int* argc, char** argv, char* ctrlname);
__declspec(dllexport) iint mw_finalize_();
__declspec(dllexport) void mw_banner_();
__declspec(dllexport) iint mw_revocap_refine_(char* filename, iint* n_refine);

__declspec(dllexport) iint mw_file_read_(char* basename);    
__declspec(dllexport) iint mw_file_read_fstr_();    
__declspec(dllexport) iint mw_file_debug_write_();  
__declspec(dllexport) iint mw_get_fstr_filename_length_mesh_();
__declspec(dllexport) iint mw_get_fstr_filename_length_control_();
__declspec(dllexport) iint mw_get_fstr_filename_length_result_();
__declspec(dllexport) iint mw_get_fstr_filename_length_restart_();
__declspec(dllexport) iint mw_get_fstr_filename_length_part_in_();
__declspec(dllexport) iint mw_get_fstr_filename_length_part_out_();
__declspec(dllexport) iint mw_get_fstr_filename_length_vis_mesh_();
__declspec(dllexport) iint mw_get_fstr_filename_length_vis_in_();
__declspec(dllexport) iint mw_get_fstr_filename_length_vis_out_();
__declspec(dllexport) iint mw_get_fstr_filename_length_cadfit_();

__declspec(dllexport) void mw_get_fstr_filename_mesh_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_control_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_result_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_restart_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_part_in_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_part_out_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_vis_mesh_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_vis_in_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_vis_out_(char name[], iint* len);
__declspec(dllexport) void mw_get_fstr_filename_cadfit_(char name[], iint* len);
__declspec(dllexport) iint mw_get_fstr_refine_num_();
__declspec(dllexport) iint mw_get_fstr_refine_type_();

__declspec(dllexport) void mw_rlt_start_(iint* step);
__declspec(dllexport) void mw_rlt_start_bin_(iint* step);
__declspec(dllexport) iint mw_rlt_print_(iint* width, const char* format, ... );// format= %d:int32, %f:double(fixed), %e:double(scientific), %s:const char*
__declspec(dllexport) void mw_rlt_end_();

__declspec(dllexport) void mw_print_avs_basis_(iint* ieq);
__declspec(dllexport) void mw_rec_avs_label_(iint* imesh, char* label, char* unit, iint* ndof);
__declspec(dllexport) void mw_rec_avs_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
__declspec(dllexport) void mw_print_avs_fem_();

__declspec(dllexport) void mw_rec_vtk_label_(iint* imesh, char* label, char* unit, iint* ndof);
__declspec(dllexport) void mw_rec_vtk_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
__declspec(dllexport) void mw_print_vtk_fem_();

__declspec(dllexport) void mw_rec_uns_label_(iint* imesh, char* label, char* unit, iint* ndof);
__declspec(dllexport) void mw_rec_uns_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
__declspec(dllexport) void mw_print_uns_fem_();


__declspec(dllexport) iint mw_file_write_res_(iint* step);
__declspec(dllexport) iint mw_set_restart_(iint* step, iint* num_of_level, iint* num_of_algebra);
__declspec(dllexport) iint mw_file_write_res_bin_(iint* step);
__declspec(dllexport) iint mw_set_restart_bin_(iint* step, iint* num_of_level, iint* num_of_algebra);

__declspec(dllexport) void mw_gene_linear_algebra_(iint* num_of_algebra, iint* global_num_of_mesh, iint* dof);
__declspec(dllexport) void mw_gene_linear_algebra_assymodel_(iint* num_of_algebra, iint* global_num_of_mesh, iint* dof, iint* num_of_contact, double* transCoeff);
__declspec(dllexport) void mw_select_algebra_(iint* ieq);

__declspec(dllexport) iint mw_matrix_add_elem_(iint* imesh,  iint* ielem,  double elem_matrix[]);
__declspec(dllexport) iint mw_matrix_add_node_(iint* imesh, iint* i_index, iint* j_index, double nodal_matrix[]);
__declspec(dllexport) void mw_matrix_clear_(iint* imesh);
__declspec(dllexport) void mw_vector_clear_(iint* imesh);
__declspec(dllexport) void mw_assy_matrix_clear_();
__declspec(dllexport) void mw_assy_vector_clear_();

__declspec(dllexport) iint mw_matrix_add_elem_24_(iint* imesh, iint* ielem, double elem_matrix[][24]);
__declspec(dllexport) iint mw_matrix_add_elem_60_(iint* imesh, iint* ielem, double elem_matrix[][60]);
__declspec(dllexport) iint mw_matrix_add_elem_12_(iint* imesh, iint* ielem, double elem_matirx[][12]);
__declspec(dllexport) iint mw_matrix_add_elem_30_(iint* imesh, iint* ielem, double elem_matirx[][30]);
__declspec(dllexport) iint mw_matrix_add_elem_18_(iint* imesh, iint* ielem, double elem_matirx[][18]);
__declspec(dllexport) iint mw_matirx_add_elem_45_(iint* imesh, iint* ielem, double elem_matirx[][45]);
__declspec(dllexport) iint mw_matirx_add_elem_20_(iint* imesh, iint* ielem, double elem_matirx[][20]);
__declspec(dllexport) iint mw_matrix_add_elem_40_(iint* imesh, iint* ielem, double elem_matirx[][40]);
__declspec(dllexport) iint mw_matrix_add_elem_15_(iint* imesh, iint* ielem, double elem_matirx[][15]);
__declspec(dllexport) iint mw_matirx_add_elem_9_(iint* imesh, iint* ielem, double elem_matirx[][9]);  
__declspec(dllexport) iint mw_matirx_add_elem_48_(iint* imesh, iint* ielem, double elem_matirx[][48]);
__declspec(dllexport) iint mw_matirx_add_elem_6_(iint* imesh, iint* ielem, double elem_matirx[][6]);  
__declspec(dllexport) iint mw_matirx_add_elem_10_(iint* imesh, iint* ielem, double elem_matirx[][10]);
__declspec(dllexport) iint mw_matrix_rhs_set_bc2_(iint* imesh, iint* node_id, iint* dof, double* diagval, double* solval);
__declspec(dllexport) iint mw_matrix_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* diagval, double* rhsval);
__declspec(dllexport) iint mw_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* value);
__declspec(dllexport) iint mw_rhs_add_bc_(iint* imesh, iint* node_id, iint* dof, double* value);
__declspec(dllexport) iint mw_nl_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* value);
__declspec(dllexport) iint mw_nl_rhs_add_bc_(iint* imesh, iint* node_id, iint* dof, double* value);

//--
// solver
//--
__declspec(dllexport) iint mw_solve_(iint* iter_max, double* tolerance, iint* method, iint* pre_condition);

__declspec(dllexport) void mw_get_solution_vector_(double buf[], iint* imesh);
__declspec(dllexport) void mw_get_solution_assy_vector_(double buf[]);
__declspec(dllexport) void mw_get_rhs_vector_(double buf[], iint* imesh);
__declspec(dllexport) void mw_get_rhs_assy_vector_(double buf[]);
__declspec(dllexport) void mw_get_rhs_load_(double buf[], iint* imesh);
__declspec(dllexport) void mw_get_rhs_assy_load_(double buf[]);
__declspec(dllexport) double mw_get_solution_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
__declspec(dllexport) double mw_get_rhs_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
__declspec(dllexport) iint mw_get_solution_assy_vector_dof_(iint* imesh);
__declspec(dllexport) iint mw_get_rhs_assy_vector_dof_(iint* imesh);

__declspec(dllexport) void mw_dump_assy_matrix_();// Debug : AssyMatrix display
__declspec(dllexport) void mw_dump_rhs_assy_vector_();// Debug : RHSAssyVector display

__declspec(dllexport) void mw_matvec_assy_(iint* mglevel, iint* ieq, double assy_x[], double assy_y[]);
__declspec(dllexport) void mw_matvec_(iint* mglevel, iint* imesh, iint* ieq, double x[], double y[]);



__declspec(dllexport) iint mw_refine_(iint* num_of_refine);
__declspec(dllexport) iint mw_mg_construct_(iint* num_of_refine);
__declspec(dllexport) void mw_finalize_refine_();      
__declspec(dllexport) void mw_finalize_mg_construct_();

__declspec(dllexport) iint mw_get_num_of_assemble_model_();
__declspec(dllexport) void mw_select_assemble_model_(iint* mglevel);
__declspec(dllexport) iint mw_get_num_of_mesh_part_();
__declspec(dllexport) void mw_select_mesh_part_with_id_(iint* mesh_id);
__declspec(dllexport) void mw_select_mesh_part_(iint* index);
__declspec(dllexport) iint mw_get_mesh_part_id_(iint* mesh_index);//-- mesh_index to  mesh_id
__declspec(dllexport) iint mw_get_mesh_part_index_(iint* mesh_id);//-- mesh_id    to  mesh_index
__declspec(dllexport) void mw_select_element_with_id_(iint* elem_id);
__declspec(dllexport) void mw_select_element_(iint* index);
__declspec(dllexport) iint mw_get_element_type_();
__declspec(dllexport) iint mw_get_num_of_element_vert_();
__declspec(dllexport) void mw_get_element_vert_node_id_(iint v_node_id[]);
__declspec(dllexport) void mw_get_element_vert_node_index_(iint v_node_index[]);
__declspec(dllexport) iint mw_get_num_of_element_edge_();
__declspec(dllexport) iint mw_get_num_of_element_face_();
__declspec(dllexport) void mw_get_element_edge_node_id_(iint v_node_id[]);

__declspec(dllexport) iint mw_get_element_face_element_id_(iint iface);
__declspec(dllexport) iint mw_get_num_Of_element_edge_element_(iint iedge);
__declspec(dllexport) iint mw_get_element_edge_element_id_(iint iedge, iint i);
__declspec(dllexport) void mw_setup_neighbors_();

__declspec(dllexport) void mw_get_node_coord_(iint* node_id, double* x, double* y, double* z);
__declspec(dllexport) iint mw_get_dof_(iint* node_id);
__declspec(dllexport) iint mw_get_dof_scalar_(iint* node_id);
__declspec(dllexport) iint mw_get_dof_vector_(iint* node_id);
__declspec(dllexport) iint mw_get_num_of_node_();
__declspec(dllexport) iint mw_get_num_of_node_with_mesh_(iint* imesh);
__declspec(dllexport) iint mw_get_num_of_element_();
__declspec(dllexport) iint mw_get_num_of_element_with_mesh_(iint* imesh);
__declspec(dllexport) iint mw_get_node_id_(iint* index);
__declspec(dllexport) iint mw_get_element_id_(iint* index);
__declspec(dllexport) iint mw_get_node_index_(iint* id);
__declspec(dllexport) iint mw_get_element_index_(iint* id);

__declspec(dllexport) iint mw_get_num_of_parent_node_(iint* id, iint* mglevel);
__declspec(dllexport) iint mw_get_parent_node_id_(iint* id, iint* mglevel, iint* index);

__declspec(dllexport) void mw_construct_node_connect_fem_(iint* node_id);
__declspec(dllexport) void mw_get_node_connect_fem_size_(iint* num_itemu, iint* num_iteml);
__declspec(dllexport) void mw_get_node_connect_fem_item_(iint itemU[], iint itemL[]);
__declspec(dllexport) iint mw_get_num_of_aggregate_element_(iint* node_id);
__declspec(dllexport) iint mw_get_aggregate_element_id_(iint* node_id, iint* index);
__declspec(dllexport) iint mw_nodetype_s_();
__declspec(dllexport) iint mw_nodetype_v_();
__declspec(dllexport) iint mw_nodetype_sv_();

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

__declspec(dllexport) iint mw_fistr_elemtype_to_mw3_elemtype_(iint* fistr_elemtype);
__declspec(dllexport) iint mw_mw3_elemtype_to_fistr_elemtype_(iint* mw3_elemtype);

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
__declspec(dllexport) void mw_dndx_(iint* elem_type, iint* num_of_integ, iint* ielem, double dndx[]);
__declspec(dllexport) void mw_det_jacobian_(iint* elem_type, iint* num_of_integ, iint* igauss, double* det_j);
__declspec(dllexport) void mw_weight_(iint* elem_type, iint* num_of_integ, iint* igauss, double* w);

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
// 境界条件
//--
__declspec(dllexport) iint mw_get_num_of_boundary_bnode_mesh_();
__declspec(dllexport) iint mw_get_num_of_boundary_bface_mesh_();
__declspec(dllexport) iint mw_get_num_of_boundary_bedge_mesh_();
__declspec(dllexport) iint mw_get_num_of_boundary_bvolume_mesh_();
__declspec(dllexport) iint mw_get_bnd_type_bnode_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_bnd_type_bface_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_bnd_type_bedge_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_bnd_type_bvolume_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_neumann_type_();
__declspec(dllexport) iint mw_get_dirichlet_type_();
__declspec(dllexport) iint mw_get_num_of_bnode_in_bnode_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_bnode_in_bface_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_bnode_in_bedge_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_bnode_in_bvolume_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_dof_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_dof_in_bface_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_dof_in_bedge_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_num_of_dof_in_bvolume_mesh_(iint* ibmesh);
__declspec(dllexport) iint mw_get_dof_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* idof);
__declspec(dllexport) iint mw_get_dof_bface_mesh_(iint* ibmesh, iint* idof);
__declspec(dllexport) iint mw_get_dof_bedge_mesh_(iint* ibmesh, iint* idof);
__declspec(dllexport) iint mw_get_dof_bvolume_mesh_(iint* ibmesh, iint* idof);
__declspec(dllexport) double mw_get_bnode_value_in_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* dof);
__declspec(dllexport) double mw_get_bnode_value_in_bface_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
__declspec(dllexport) double mw_get_bnode_value_in_bedge_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
__declspec(dllexport) double mw_get_bnode_value_in_bvolume_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
__declspec(dllexport) iint mw_get_node_id_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_node_id_in_bface_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_node_id_in_bedge_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_node_id_in_bvolume_mesh_(iint* ibmesh, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_bface_(iint* ibmesh);
__declspec(dllexport) double mw_get_bface_value_(iint* ibmesh, iint* ibface, iint* dof);
__declspec(dllexport) iint mw_get_num_of_bedge_(iint* ibmesh);
__declspec(dllexport) double mw_get_bedge_value_(iint* ibmesh, iint* ibedge, iint* dof);
__declspec(dllexport) iint mw_get_num_of_bvolume_(iint* ibmesh);
__declspec(dllexport) double mw_get_bvolume_value_(iint* ibmesh, iint* ibvol, iint* dof);
__declspec(dllexport) iint mw_get_num_of_node_bface_(iint* ibmesh, iint* ibface);
__declspec(dllexport) iint mw_get_node_id_bface_(iint* ibmesh, iint* ibface, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_node_bedge_(iint* ibmesh, iint* ibedge);
__declspec(dllexport) iint mw_get_node_id_bedge_(iint* ibmesh, iint* ibedge, iint* ibnode);
__declspec(dllexport) iint mw_get_num_of_node_bvolume_(iint* ibmesh, iint* ibvol);
__declspec(dllexport) iint mw_get_node_id_bvolume_(iint* ibmesh, iint* ibvol, iint* ibnode);
__declspec(dllexport) iint mw_get_bnode_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bnode_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_bface_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bface_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_bvolume_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bvolume_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_bedge_mesh_namelength_(iint* ibmesh);
__declspec(dllexport) void mw_get_bedge_mesh_name_(iint* ibmesh, char* name, iint* name_len);
__declspec(dllexport) iint mw_get_edge_id_bedge_(iint*ibmesh, iint* ibedge);
__declspec(dllexport) iint mw_get_elem_id_bedge_(iint* ibmesh, iint* ibedge);
__declspec(dllexport) iint mw_get_face_id_bface_(iint* ibmesh, iint* ibface);
__declspec(dllexport) iint mw_get_elem_id_bface_(iint* ibmesh, iint* ibface);
__declspec(dllexport) iint mw_get_elem_id_bvolume_(iint* ibmesh, iint* ibvol);

//--
// MPI
//--
__declspec(dllexport) MPI_Datatype mw_mpi_int_();
__declspec(dllexport) MPI_Datatype mw_mpi_iint_();
__declspec(dllexport) MPI_Datatype mw_mpi_uiint_();
__declspec(dllexport) MPI_Datatype mw_mpi_double_();
__declspec(dllexport) MPI_Comm mw_mpi_comm_();
__declspec(dllexport) MPI_Op mw_mpi_sum_();
__declspec(dllexport) MPI_Op mw_mpi_max_();
__declspec(dllexport) MPI_Op mw_mpi_min_();

__declspec(dllexport) int mw_get_rank_();
__declspec(dllexport) int mw_get_num_of_process_();
__declspec(dllexport) void mw_allreduce_r_(double val[], iint* val_size, MPI_Op* op);
__declspec(dllexport) void mw_allreduce_i_(iint val[], iint* val_size, MPI_Op* op);

__declspec(dllexport) iint mw_barrier_();
__declspec(dllexport) iint mw_abort_(int* error);
__declspec(dllexport) iint mw_allgather_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt);
__declspec(dllexport) iint mw_allgather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt);

__declspec(dllexport) iint mw_gather_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt, int* root);
__declspec(dllexport) iint mw_gather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt, int* root);

__declspec(dllexport) iint mw_scatter_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt, int* root);
__declspec(dllexport) iint mw_scatter_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt, int* root);

__declspec(dllexport) iint mw_bcast_i_(iint buf[], iint* cnt, int* root);
__declspec(dllexport) iint mw_bcast_r_(double buf[], iint* cnt, int* root);
__declspec(dllexport) iint mw_bcast_s_(char buf[], iint* cnt, int* root);

__declspec(dllexport) iint mw_get_num_of_neibpe_(iint* imesh);
__declspec(dllexport) iint mw_get_transrank_(iint* imesh, iint* ipe);

__declspec(dllexport) iint mw_send_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank);
__declspec(dllexport) iint mw_send_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank);

__declspec(dllexport) iint mw_recv_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank);
__declspec(dllexport) iint mw_recv_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank);

__declspec(dllexport) void mw_send_recv_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank);
__declspec(dllexport) void mw_send_recv_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank);

//--
// 通信テーブル:CommMesh2
//--
__declspec(dllexport) iint mw_get_num_of_comm_mesh_();
__declspec(dllexport) iint mw_get_num_of_comm_node_(iint* icmesh);
__declspec(dllexport) iint mw_get_node_id_comm_node_(iint* icmesh, iint* icnode);

//--
// 接合面:ContactMesh
//--
__declspec(dllexport) iint mw_get_num_of_contact_();
__declspec(dllexport) iint mw_get_contact_id_(iint* icont);


__declspec(dllexport) iint mw_get_num_of_elementgroup_();
__declspec(dllexport) iint mw_get_num_of_element_id_(iint* igrp);
__declspec(dllexport) iint mw_get_element_id_with_elementgroup_(iint* igrp, iint* index);
__declspec(dllexport) iint mw_get_elementgroup_name_length_(iint* igrp);
__declspec(dllexport) void mw_get_elementgroup_name_(iint* igrp, char* name, iint* name_len);


__declspec(dllexport) void mw_logger_set_mode_(iint* mode);
__declspec(dllexport) void mw_logger_set_device_(iint* mode, iint* device);
__declspec(dllexport) void mw_logger_info_mssg_ (iint* mode, const char* message);
__declspec(dllexport) void mw_logger_info_(iint* mode, const char* format, ...);
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
