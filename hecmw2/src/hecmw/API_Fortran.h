/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/API_Fortran.h
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
#ifndef API_Fortran_AF05492D_2BFF_49fd_AEE3_818B0B426793
#define API_Fortran_AF05492D_2BFF_49fd_AEE3_818B0B426793
#include "TypeDef.h"

#include <cstring>
#include <cmath>
#include <cstdlib>
using namespace std;
#ifdef __cplusplus
extern "C" {
#endif

    iint mw_initialize_(int* argc, char** argv);
    iint mw_initialize_fstr_(int* argc, char** argv, char* ctrlname);
    iint mw_finalize_();

    void mw_banner_();
    iint mw_revocap_refine_(char* filename, iint* n_refine);

    iint mw_file_read_(char* basename);
    iint mw_file_read_fstr_();
    iint mw_file_debug_write_();

    iint mw_get_fstr_filename_length_mesh_();
    iint mw_get_fstr_filename_length_control_();
    iint mw_get_fstr_filename_length_result_();
    iint mw_get_fstr_filename_length_restart_();
    iint mw_get_fstr_filename_length_part_in_();
    iint mw_get_fstr_filename_length_part_out_();
    iint mw_get_fstr_filename_length_vis_mesh_();
    iint mw_get_fstr_filename_length_vis_in_();
    iint mw_get_fstr_filename_length_vis_out_();
    iint mw_get_fstr_filename_length_cadfit_();

    void mw_get_fstr_filename_mesh_(char name[], iint* len);
    void mw_get_fstr_filename_control_(char name[], iint* len);
    void mw_get_fstr_filename_result_(char name[], iint* len);
    void mw_get_fstr_filename_restart_(char name[], iint* len);
    void mw_get_fstr_filename_part_in_(char name[], iint* len);
    void mw_get_fstr_filename_part_out_(char name[], iint* len);
    void mw_get_fstr_filename_vis_mesh_(char name[], iint* len);
    void mw_get_fstr_filename_vis_in_(char name[], iint* len);
    void mw_get_fstr_filename_vis_out_(char name[], iint* len);
    void mw_get_fstr_filename_cadfit_(char name[], iint* len);
    iint mw_get_fstr_refine_num_();
    iint mw_get_fstr_refine_type_();

    void mw_rlt_start_(iint* step);
    void mw_rlt_start_bin_(iint* step);
    iint mw_rlt_print_(iint* width, const char* format, ... );//format= %d:int32, %f:double(fixed), %e:double(scientific), %s:const char*
    void mw_rlt_end_();

    void mw_print_avs_basis_(iint* ieq);
    void mw_rec_avs_label_(iint* imesh, char* label, char* unit, iint* ndof);
    void mw_rec_avs_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
    void mw_print_avs_fem_();

    void mw_rec_vtk_label_(iint* imesh, char* label, char* unit, iint* ndof);
    void mw_rec_vtk_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
    void mw_print_vtk_fem_();

    void mw_rec_uns_label_(iint* imesh, char* label, char* unit, iint* ndof);
    void mw_rec_uns_variable_(iint* imesh, iint* num_of_node, char* label, double* value);
    void mw_print_uns_fem_();


    iint mw_file_write_res_(iint* step);
    iint mw_set_restart_(iint* step, iint* num_of_level, iint* num_of_algebra);
    iint mw_file_write_res_bin_(iint* step);
    iint mw_set_restart_bin_(iint* step, iint* num_of_level, iint* num_of_algebra);

    void mw_gene_linear_algebra_(iint* num_of_algebra, iint* global_num_of_mesh, iint* dof);
    void mw_gene_linear_algebra_assymodel_(iint* num_of_algebra, iint* global_num_of_mesh, iint* dof, iint* num_of_contact, double* transCoeff);
    void mw_select_algebra_(iint* ieq);

    iint mw_matrix_add_elem_(iint* imesh, iint* ielem,  double elem_matrix[]);
    iint mw_matrix_add_node_(iint* imesh, iint* i_index, iint* j_index, double nodal_matrix[]);
    void mw_matrix_clear_(iint* imesh);
    void mw_vector_clear_(iint* imesh);
    void mw_assy_matrix_clear_();
    void mw_assy_vector_clear_();

    iint mw_matrix_add_elem_24_(iint* imesh, iint* ielem, double elem_matrix[][24]);
    iint mw_matrix_add_elem_60_(iint* imesh, iint* ielem, double elem_matrix[][60]);
    iint mw_matrix_add_elem_12_(iint* imesh, iint* ielem, double elem_matirx[][12]);
    iint mw_matrix_add_elem_30_(iint* imesh, iint* ielem, double elem_matirx[][30]);
    iint mw_matrix_add_elem_18_(iint* imesh, iint* ielem, double elem_matirx[][18]);
    iint mw_matirx_add_elem_45_(iint* imesh, iint* ielem, double elem_matirx[][45]);
    iint mw_matirx_add_elem_20_(iint* imesh, iint* ielem, double elem_matirx[][20]);
    iint mw_matrix_add_elem_40_(iint* imesh, iint* ielem, double elem_matirx[][40]);
    iint mw_matrix_add_elem_15_(iint* imesh, iint* ielem, double elem_matirx[][15]);
    iint mw_matirx_add_elem_9_(iint* imesh, iint* ielem, double elem_matirx[][9]);
    iint mw_matirx_add_elem_48_(iint* imesh, iint* ielem, double elem_matirx[][48]);
    iint mw_matirx_add_elem_6_(iint* imesh, iint* ielem, double elem_matirx[][6]);
    iint mw_matirx_add_elem_10_(iint* imesh, iint* ielem, double elem_matirx[][10]);
    iint mw_matrix_rhs_set_bc2_(iint* imesh, iint* node_id, iint* dof, double* diagval, double* solval);
    iint mw_matrix_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* diagval, double* rhsval);
    iint mw_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* value);
    iint mw_rhs_add_bc_(iint* imesh, iint* node_id, iint* dof, double* value);
    iint mw_nl_rhs_set_bc_(iint* imesh, iint* node_id, iint* dof, double* value);//---- 節点集中荷重(Nodal_Load)として扱う:Rank大には値は入らない.
    iint mw_nl_rhs_add_bc_(iint* imesh, iint* node_id, iint* dof, double* value);//---- 節点集中荷重(Nodal_Load)として扱う:Rank大には値は入らない.

//--
// solver
//--
    iint mw_solve_(iint* iter_max, double* tolerance, iint* method, iint* pre_condition);

    void mw_get_solution_vector_(double buf[], iint* imesh);
    void mw_get_solution_assy_vector_(double buf[]);
    void mw_get_rhs_vector_(double buf[], iint* imesh);
    void mw_get_rhs_assy_vector_(double buf[]);
    void mw_get_rhs_load_(double buf[], iint* imesh);
    void mw_get_rhs_assy_load_(double buf[]);
    double mw_get_solution_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
    double mw_get_rhs_assy_vector_val_(iint* imesh, iint* inode, iint* idof);
    iint mw_get_solution_assy_vector_dof_(iint* imesh);
    iint mw_get_rhs_assy_vector_dof_(iint* imesh);

    void mw_dump_assy_matrix_();// Debug : AssyMatrix display
    void mw_dump_rhs_assy_vector_();// Debug : RHSAssyVector display

//--
// y = A*x : x update, y sumup
//--
    void mw_matvec_assy_(iint* mglevel, iint* ieq, double assy_x[], double assy_y[]);
    void mw_matvec_(iint* mglevel, iint* imesh, iint* ieq, double x[], double y[]);


    iint mw_refine_(iint* num_of_refine);
    iint mw_mg_construct_(iint* num_of_refine);
    void mw_finalize_refine_();
    void mw_finalize_mg_construct_();


    iint mw_get_num_of_assemble_model_();
    void mw_select_assemble_model_(iint* mglevel);
    iint mw_get_num_of_mesh_part_();
    void mw_select_mesh_part_with_id_(iint* mesh_id);
    void mw_select_mesh_part_(iint* index);
    iint mw_get_mesh_part_id_(iint* mesh_index);//-- mesh_index to  mesh_id
    iint mw_get_mesh_part_index_(iint* mesh_id);//-- mesh_id    to  mesh_index
    void mw_select_element_with_id_(iint* elem_id);
    void mw_select_element_(iint* index);
    iint mw_get_element_type_();
    iint mw_get_num_of_element_vert_();
    void mw_get_element_vert_node_id_(iint v_node_id[]);
    void mw_get_element_vert_node_index_(iint v_node_index[]);
    iint mw_get_num_of_element_edge_();
    iint mw_get_num_of_element_face_();
    void mw_get_element_edge_node_id_(iint v_node_id[]);

    iint mw_get_element_face_element_id_(iint iface);
    iint mw_get_num_Of_element_edge_element_(iint iedge);
    iint mw_get_element_edge_element_id_(iint iedge, iint i);
    void mw_setup_neighbors_();

    void mw_get_node_coord_(iint* node_id, double* x, double* y, double* z);
    iint mw_get_dof_(iint* node_id);
    iint mw_get_dof_scalar_(iint* node_id);
    iint mw_get_dof_vector_(iint* node_id);
    iint mw_get_num_of_node_();
    iint mw_get_num_of_node_with_mesh_(iint* imesh);
    iint mw_get_num_of_element_();
    iint mw_get_num_of_element_with_mesh_(iint* imesh);
    iint mw_get_node_id_(iint* index);
    iint mw_get_element_id_(iint* index);
    iint mw_get_node_index_(iint* id);
    iint mw_get_element_index_(iint* id);

    iint mw_get_num_of_parent_node_(iint* id, iint* mglevel);
    iint mw_get_parent_node_id_(iint* id, iint* mglevel, iint* index);

    void mw_construct_node_connect_fem_(iint* node_id);
    void mw_get_node_connect_fem_size_(iint* num_itemu, iint* num_iteml);
    void mw_get_node_connect_fem_item_(iint itemU[], iint itemL[]);
    iint mw_get_num_of_aggregate_element_(iint* node_id);
    iint mw_get_aggregate_element_id_(iint* node_id, iint* index);
    iint mw_nodetype_s_();
    iint mw_nodetype_v_();
    iint mw_nodetype_sv_();

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

    iint mw_fistr_elemtype_to_mw3_elemtype_(iint* fistr_elemtype);
    iint mw_mw3_elemtype_to_fistr_elemtype_(iint* mw3_elemtype);

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
    void mw_dndx_(iint* elem_type, iint* num_of_integ, iint* ielem, double dndx[]);
    void mw_det_jacobian_(iint* elem_type, iint* num_of_integ, iint* igauss, double* det_j);
    void mw_weight_(iint* elem_type, iint* num_of_integ, iint* igauss, double* w);

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
// 境界条件
//--
    iint mw_get_num_of_boundary_bnode_mesh_();
    iint mw_get_num_of_boundary_bface_mesh_();
    iint mw_get_num_of_boundary_bedge_mesh_();
    iint mw_get_num_of_boundary_bvolume_mesh_();
    iint mw_get_bnd_type_bnode_mesh_(iint* ibmesh);
    iint mw_get_bnd_type_bface_mesh_(iint* ibmesh);
    iint mw_get_bnd_type_bedge_mesh_(iint* ibmesh);
    iint mw_get_bnd_type_bvolume_mesh_(iint* ibmesh);
    iint mw_get_neumann_type_();
    iint mw_get_dirichlet_type_();
    iint mw_get_num_of_bnode_in_bnode_mesh_(iint* ibmesh);
    iint mw_get_num_of_bnode_in_bface_mesh_(iint* ibmesh);
    iint mw_get_num_of_bnode_in_bedge_mesh_(iint* ibmesh);
    iint mw_get_num_of_bnode_in_bvolume_mesh_(iint* ibmesh);
    iint mw_get_num_of_dof_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
    iint mw_get_num_of_dof_in_bface_mesh_(iint* ibmesh);
    iint mw_get_num_of_dof_in_bedge_mesh_(iint* ibmesh);
    iint mw_get_num_of_dof_in_bvolume_mesh_(iint* ibmesh);
    iint mw_get_dof_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* idof);
    iint mw_get_dof_bface_mesh_(iint* ibmesh, iint* idof);
    iint mw_get_dof_bedge_mesh_(iint* ibmesh, iint* idof);
    iint mw_get_dof_bvolume_mesh_(iint* ibmesh, iint* idof);
    double mw_get_bnode_value_in_bnode_mesh_(iint* ibmesh, iint* ibnode, iint* dof);
    double mw_get_bnode_value_in_bface_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
    double mw_get_bnode_value_in_bedge_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
    double mw_get_bnode_value_in_bvolume_mesh_(iint* ibmesh, iint* ibnode, iint* dof, iint* mglevel);
    iint mw_get_node_id_in_bnode_mesh_(iint* ibmesh, iint* ibnode);
    iint mw_get_node_id_in_bface_mesh_(iint* ibmesh, iint* ibnode);
    iint mw_get_node_id_in_bedge_mesh_(iint* ibmesh, iint* ibnode);
    iint mw_get_node_id_in_bvolume_mesh_(iint* ibmesh, iint* ibnode);
    iint mw_get_num_of_bface_(iint* ibmesh);
    double mw_get_bface_value_(iint* ibmesh, iint* ibface, iint* dof);
    iint mw_get_num_of_bedge_(iint* ibmesh);
    double mw_get_bedge_value_(iint* ibmesh, iint* ibedge, iint* dof);
    iint mw_get_num_of_bvolume_(iint* ibmesh);
    double mw_get_bvolume_value_(iint* ibmesh, iint* ibvol, iint* dof);
    iint mw_get_num_of_node_bface_(iint* ibmesh, iint* ibface);
    iint mw_get_node_id_bface_(iint* ibmesh, iint* ibface, iint* ibnode);
    iint mw_get_num_of_node_bedge_(iint* ibmesh, iint* ibedge);
    iint mw_get_node_id_bedge_(iint* ibmesh, iint* ibedge, iint* ibnode);
    iint mw_get_num_of_node_bvolume_(iint* ibmesh, iint* ibvol);
    iint mw_get_node_id_bvolume_(iint* ibmesh, iint* ibvol, iint* ibnode);
    iint mw_get_bnode_mesh_namelength_(iint* ibmesh);
    void mw_get_bnode_mesh_name_(iint* ibmesh, char* name, iint* name_len);
    iint mw_get_bface_mesh_namelength_(iint* ibmesh);
    void mw_get_bface_mesh_name_(iint* ibmesh, char* name, iint* name_len);
    iint mw_get_bvolume_mesh_namelength_(iint* ibmesh);
    void mw_get_bvolume_mesh_name_(iint* ibmesh, char* name, iint* name_len);
    iint mw_get_bedge_mesh_namelength_(iint* ibmesh);
    void mw_get_bedge_mesh_name_(iint* ibmesh, char* name, iint* name_len);
    iint mw_get_edge_id_bedge_(iint* ibmesh, iint* ibedge);
    iint mw_get_elem_id_bedge_(iint* ibmesh, iint* ibedge);
    iint mw_get_face_id_bface_(iint* ibmesh, iint* ibface);
    iint mw_get_elem_id_bface_(iint* ibmesh, iint* ibface);
    iint mw_get_elem_id_bvolume_(iint* ibmesh, iint* ibvol);
//--
// MPI
//--
    MPI_Datatype mw_mpi_int_();
    MPI_Datatype mw_mpi_iint_();
    MPI_Datatype mw_mpi_uiint_();
    MPI_Datatype mw_mpi_double_();
    MPI_Comm mw_mpi_comm_();
    MPI_Op mw_mpi_sum_();
    MPI_Op mw_mpi_max_();
    MPI_Op mw_mpi_min_();

    int mw_get_rank_();
    int mw_get_num_of_process_();
    void mw_allreduce_r_(double val[], iint* val_size, MPI_Op* op);
    void mw_allreduce_i_(iint val[], iint* val_size, MPI_Op* op);

    iint mw_barrier_();
    iint mw_abort_(int* error);
    iint mw_allgather_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt);
    iint mw_allgather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt);

    iint mw_gather_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt, int* root);
    iint mw_gather_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt, int* root);

    iint mw_scatter_r_(double sendbuf[], iint* sendcnt, double recvbuf[], iint* recvcnt, int* root);
    iint mw_scatter_i_(iint sendbuf[], iint* sendcnt, iint recvbuf[], iint* recvcnt, int* root);

    iint mw_bcast_i_(iint buf[], iint* cnt, int* root);
    iint mw_bcast_r_(double buf[], iint* cnt, int* root);
    iint mw_bcast_s_(char buf[], iint* cnt, int* root);

    iint mw_get_num_of_neibpe_(iint* imesh);
    iint mw_get_transrank_(iint* imesh, iint* ipe);

    iint mw_send_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank);
    iint mw_send_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank);

    iint mw_recv_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank);
    iint mw_recv_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank);

    void mw_send_recv_r_(double buf[], iint* num_of_node, iint* dof_size, int* trans_rank);
    void mw_send_recv_i_(iint buf[], iint* num_of_node, iint* dof_size, int* trans_rank);

    void mw_sumup_(iint* mglevel, iint* imesh, double vec[], iint* num_of_dof);
//--
// 通信テーブル:CommMesh2
//--
    iint mw_get_num_of_comm_mesh_();
    iint mw_get_num_of_comm_node_(iint* icmesh);
    iint mw_get_node_id_comm_node_(iint* icmesh, iint* icnode);

//--
// 接合面:ContactMesh
//--
    iint mw_get_num_of_contact_();
    iint mw_get_contact_id_(iint* icont);


    iint mw_get_num_of_elementgroup_();
    iint mw_get_num_of_element_id_(iint* igrp);
    iint mw_get_element_id_with_elementgroup_(iint* igrp, iint* index);
    iint mw_get_elementgroup_name_length_(iint* igrp);
    void mw_get_elementgroup_name_(iint* igrp, char* name, iint* name_len);

    void mw_logger_set_mode_(iint* mode);
    void mw_logger_set_device_(iint* mode, iint* device);
    void mw_logger_info_mssg_(iint* mode, const char* message);
    void mw_logger_info_(iint* mode, const char* format, ...);
    iint mw_get_error_mode_();
    iint mw_get_warn_mode_();
    iint mw_get_info_mode_();
    iint mw_get_debug_mode_();
    iint mw_get_disk_device_();
    iint mw_get_display_device_();
#ifdef __cplusplus
}
#endif
#endif
