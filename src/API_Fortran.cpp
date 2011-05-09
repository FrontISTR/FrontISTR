//
//  API_Fortran.cpp
//
//                  2009.4.20
//                  2009.2.27
//                  k.Takeda
#include "API_Fortran.h"
#include "HEC_MW3.h"

pmw::CMWMain *pMW;

//----
// 1. HEC_MW3 construct & destruct
//----
void mw_initialize_(){
    pMW = pmw::CMWMain::Instance();
}
void mw_finalize_(){
    ;
}


//----
// 2. file i/o API
//----
void mw_file_import_syscontrol_(int* modetype ){
    ;//
}
void mw_file_export_(int* modetype){
    ;//
}


//----
// 3.VISUALIZER API
//----
void mw_visualize_(int* assy_id, int* result_type){
    ;
}


//----
// 4. shape function API
//----
// gaussian integral point API
void mw_gauss_(int* num_of_points, double *gzi_1d, double* weight){
    ;
}
void mw_gauss_tetra_on_pt_(int* order, double* tetra_coord, double* weight){
    ;
}
void mw_gauss_tetra_(int* order, double** tetra_coord, double* weight){
    ;
}
void mw_gauss_tri_on_pt_(int* order, int* index, double* tri_coord, double* weight){
    ;
}
void mw_gauss_tri_(int* order, double** tri_cood, double* weight){
    ;
}


//----
//
//---
void mw_elem_shapefunc_g_on_pt_(int* elem_type, double* gzi_coord, double* n){
    ;
}
void mw_elem_shapefunc_g_(int* elem_type, int* integ_order, double** n, int* num_of_point){
    ;
}


//----
//
//---
void mw_elem_dndr_on_pt_(int* elemtype, double* gzi_coord, double** dNdr){
    ;
}
void mw_elem_dndr_(int* elemtype, int* integ_order, double*** dNdr){
    ;
}


//----
//
//---
void mw_elem_shapefunc_g2_on_pt_(int* elem_type, double* gzi_coord, double* N , double **dNdr){
    ;
}
void mw_elem_shapefunc_g2_(int* elem_type, int* integ_order, double** N, double*** dNdr){
    ;
}


//----
//
//---
void mw_elem_jacobian1_on_pt_(int* assy_id, int* elem_id, double* gzi_coord, double** jacobian){
    ;
}
void mw_elem_jacobian1_(int* assy_id, int* elem_id, int* integ_order, double*** jacobian){
    ;
}


//----
//
//---
void mw_elem_jacobian2_on_pt_(int* assy_id, int* elem_id, double** dNdr, double** jacobian){
    ;
}
void mw_elem_jacobian2_(int* assy_id, int* elem_id, int* integ_order, double*** dNdr, double*** jacobian){
    ;
}


//----
//
//---
void mw_matrix33_inv_(double** mat, double** invmat){
    ;
}

//----
//
//----
void mw_matrix33_det_(double** mat, double* det){
    ;
}

//----
// Element Volume
//----
void mw_element_volume_(int* assy_id, int* elem_id, double* volume){
    ;
}

//----
// Shape Function grad  (x,y,z)
//----
void mw_elem_dndx1_on_pt_(int* assy_id, int* elem_id, double* gzi_coord, double** dNdx){
    ;
}
void mw_elem_dndx1_(int* assy_id, int* elem_id, int* integ_order, double*** dNdx){
    ;
}
void mw_elem_dndx2_on_pt_( int* assy_id,  int* elem_id,  double** dNdr,  double** jacobian_inv,  double** dNdx){
    ;
}
void mw_elem_dndx2_( int* assy_id,  int* elem_id,  int* integ_order,  double*** dNdr,  double*** jacobian_inv,  double*** dNdx){
    ;
}
void mw_elem_shapefunc_x_on_pt_( int* assy_id, int* elem_id,  double* gzi_coord,  double* N,  double** dNdx){
    ;
}
void mw_elem_shapefunc_x_( int* assy_id, int* elem_id, int* integ_order,  double** N, double*** dNdx){
    ;
}

//----
// Extrapolation
//----
void mw_value_extrapolation_( int* assy_id,  int* node_id,  int* value_type){
    ;
}

//----
// B matrix API
//----
void mw_elem_b_matrix1_on_pt_( int* assy_id,  int* elem_id,  int* integ_order,  int* integ_index,  double* weight,  double* det_jacobi, double** b_mat){
    ;
}
void mw_elem_b_matrix1_( int* assy_id, int* elem_id,  int* integ_order,  double* weight,  double* det_jacobi,  double*** b_mat){
    ;
}
void mw_elem_b_matrix2_on_pt_( int* assy_id, int*  elem_id, double** dNdx,  double** b_mat){
    ;
}
void mw_elem_b_matrix2_( int* assy_id, int* elem_id,  int* integ_order, double***  dNdx,  double*** b_mat){
    ;
}

//----
// disp grad tensor
//----
void mw_elem_f_vector_( int* assy_id,  int* elem_id,  int* local_node_num,  double* f_vec ){
    ;
}

//----
// velocity grad tensor L,stretch tensor D, spin tensor W API
//----
void mw_elem_l_matrix_(  int* assy_id,  int* elem_id,  int* local_node_num, double** l_mat){
    ;
}
void mw_elem_d_matrix_(  int* assy_id,  int* elem_id,  int* local_node_num, double** d_mat ){
    ;
}
void mw_elem_w_matrix_(  int* assy_id,  int* elem_id,  int* local_node_num, double** w_mat ){
    ;
}

//----
//([B] = [Z1][Z2])
//----
void mw_elem_z1_matrix_ (  int* assy_id,  int* elem_id,  int* local_node_num,   double** z1_mat ){
    ;
}
void mw_elem_z2_matrix_ (   int* assy_id,  int* elem_id,  int* local_node_num,   double** z2_mat ){
    ;
}


//----
// 5.matrix API
//----
void mw_initialize_matrix_(  int* number_of_matrix ){
    ;
}
void mw_get_num_of_matrix_( int*  number_of_matrix ){
    ;
}
void mw_finalize_matrix_( int* matrix_id ){
    ;
}

//----
// non-zero matrix construct API
//----
void mw_compress_matrix_(int* matrix_id ){
    ;
}

//----
// elem_matrix -> global matrix API
//----
void mw_matrix_add_elem_(int* matrix_id, int* assy_id,  int* elem_id,  double** elem_matrix){
    ;
}

//----
// matrix transpose
//----
void mw_matrix_transposed_( double* b_mat, double* bt_mat){
    ;
}

//----
// matrix product
//----
void mw_matrix_product_(double* b_mat, double* bt_mat, double* btb_mat){
    ;
}
void mw_matrix_product_s_(double* b_mat, double* bt_mat, double* scalar_val, double* btb_mat){
    ;
}


//----
// 6.vector API
//----
void mw_initialize_vector_ ( int* number_of_vector ){
    ;
}
void mw_get_num_of_vector_ (int*  number_of_vector ){
    ;
}
void mw_finalize_vector_ ( int* vector_id ){
    ;
}

//----
// matrix vector product API
//----
void mw_multiply_matrix_vector_ ( double** matrix,  double* vector, double* result_vector){
    ;
}

//----
// vector product, norm API
//----
void mw_vector_inner_product_ (double* vector1,  double* vector2,  double* result){
    ;
}
void mw_vector_norm_2_ ( double* vector, double* norm){
    ;
}
void mw_vector_norm_inf_ (double* vector, double* norm){
    ;
}


//----
// 7. boundary API
//----
// element load -> add for global vector
//
void mw_vector_add_elem_ ( int* vector_id,  int* assy_id,  int* elem_id, double* elem_vector){
    ;
}

// boundary set API
void mw_set_bc_ ( double** matrix,  double* rhs_vector, int* node_id, int* dof_id, double* val){
    ;
}


//----
// 8. MPC API
//----
void mw_set_mpc_ (double** matrix,  double* rhs_vector,  int* num_terms,  int* node_id,  int* dof_id,  double* coef){
    ;
}


//----
// 9. linear solver API
//----
// iterative API
//
void mw_solve_cg_ ( ){
    ;
}
void mw_solve_gmres_ ( ){
    ;
}
void mw_solve_bicgstab_ ( ){
    ;
}
void mw_solve_gpbicg_ ( ){
    ;
}

//----
// multigrid API
//----
void mw_solve_gmg_( int* cycle_type,  int* solve_type){
    ;
}

//----
// solver pre-process API
//----
void mw_pre_solve_ ( int*  pretype ){
    ;
}

//----
// direct solver API
//----
void mw_solve_direct_ ( ){
    ;
}


//----
// 10.essention API
//----
// Mesh data API
void mw_get_number_of_assy_ ( int* num_of_assy){
    ;
}
void mw_get_number_of_node_ ( int* assy_id,  int* num_of_node){
    ;
}
void mw_get_number_of_element_ ( int* assy_id,  int* num_of_element){
    ;
}

//----
// Element Type API
//----
void mw_get_type_of_element_ (int* assy_id,  int* elem_id, int* element_type){
    ;
}
void mw_get_number_of_local_node_(int* assy_id,  int* elem_id,  int* num_of_local_node){
    ;
}

//----
// Node DOF Type API
//----
void mw_get_type_of_node_ ( int* assy_id, int* elem_id,  int* local_node_id,  int* node_type){
    ;
}

//----
// Group API
//----
void mw_get_num_of_group_ ( int* group_type,  int* num_of_group){
    ;
}
void mw_get_num_of_node_nodegroup_ (int*  group_id, int*  num_of_node){
    ;
}
void mw_get_num_of_elem_elemgroup_ ( int* group_id, int*  num_of_element){
    ;
}
void mw_get_node_index_ ( int* group_id,  int* node_index){
    ;
}
void mw_get_elem_index_ ( int* group_id,  int* element_index){
    ;
}

//----
// Material API
//----
void mw_get_num_of_material_ ( int* num_of_material ){
    ;
}
void mw_get_poisson_ (int* matrial_id,  double* poisson_val){
    ;
}
void mw_get_young_coeff_ ( int* material_id,  double* val_e ){
    ;
}
void mw_get_thermal_conduct_ (int*  material_id,   double* val_conductivity_k){
    ;
}
void mw_get_thermal_expansion_ ( int* material_id,   double* val_expansion){
    ;
}

//----
// global info API
//----
void mw_get_max_ (int* val_type,  double* max_val){
    ;
}
void mw_get_min_ (int* val_type, double* min_val){
    ;
}
void mw_get_average_ (int* val_type,  double* ave_val){
    ;
}
void mw_get_mw_sys_ (int* num_of_precess){
    ;
}

//----
// Section API
//----

void mw_get_section_prop_count_ ( int* section_tag_type, int*  num_of_section_property){
    ;
}
void mw_get_section_prop_ (int* section_tag_type, int* sectin_prop_no,  double*  property){
    ;
}


//----
// 11.Utility PI
//----
void mw_get_memory_size_ (int* assy_id,  int* num_of_byte){
    ;
}

//----
// Logger API
//----
void mw_logger_mode_ ( int* state_type){
    ;
}
void mw_logger_property_ ( int* state_type,  int* output_type){
    ;
}
void mw_logger_monitor_ ( int* state_type, int* id,  double* value, char* message, int* str_len){
    ;
}
void mw_logger_info_ (int* state_type, char* message, int* str_len){
    ;
}
