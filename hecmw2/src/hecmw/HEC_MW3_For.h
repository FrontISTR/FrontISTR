!
! ----------------------------------------------------------
!|
!| Software Name :HEC-MW Ver 4.3 beta
!|
!|   ../src/HEC_MW3_For.h
!|
!|                     Written by T.Takeda,    2012/06/01
!|                                Y.Sato,      2012/06/01
!|                                K.Goto,      2010/01/12
!|                                K.Matsubara, 2010/06/01
!|
!|   Contact address : IIS, The University of Tokyo CISS
!|
! ----------------------------------------------------------
!
!----
! HEC_MW3 construct & destruct
!----
integer mw_initialize
integer mw_initialize_fstr
integer mw_finalize
!----
! banner
!----
external mw_banner
!----
! revocap_refiner
!----
integer mw_revocap_refine
!----
! file i/o 
!----
integer mw_file_read
integer mw_file_read_bin
integer mw_file_read_fstr
integer mw_file_read_bin_fstr
integer mw_file_debug_write
!----
! fstr (hecmw_ctrl)
!----
integer mw_get_fstr_filename_length_mesh
integer mw_get_fstr_filename_length_control
integer mw_get_fstr_filename_length_result
integer mw_get_fstr_filename_length_restart
integer mw_get_fstr_filename_length_part_in
integer mw_get_fstr_filename_length_part_out
integer mw_get_fstr_filename_length_vis_mesh
integer mw_get_fstr_filename_length_vis_in
integer mw_get_fstr_filename_length_vis_out
integer mw_get_fstr_filename_length_cadfit

external mw_get_fstr_filename_mesh
external mw_get_fstr_filename_control
external mw_get_fstr_filename_result
external mw_get_fstr_filename_restart
external mw_get_fstr_filename_part_in
external mw_get_fstr_filename_part_out
external mw_get_fstr_filename_vis_mesh
external mw_get_fstr_filename_vis_in
external mw_get_fstr_filename_vis_out
external mw_get_fstr_filename_cadfit
integer mw_get_fstr_refine_num
integer mw_get_fstr_refine_type
!----
! result, format= %d:int32*, %f:double*(fixed), %e:double*(scientific), %s:const char*
!----
external mw_rlt_start
external mw_rlt_start_bin
integer  mw_rlt_print
external mw_rlt_end
!----
! output *.inp
!----
external mw_print_avs
external mw_rec_avs_label
external mw_rec_avs_variable
external mw_print_avs_fem
!----
! output *.vtu
!----
external mw_rec_vtk_label
external mw_rec_vtk_variable
external mw_print_vtk_fem
!----
! output *.uns
!----
external mw_rec_uns_label
external mw_rec_uns_variable
external mw_print_uns_fem
!----
! restart, linear_algebra_equation info
!----
integer mw_file_write_res
integer mw_set_restart
integer mw_file_write_res_bin
integer mw_set_restart_bin
!----
! linear solver API
!----
external mw_gene_linear_algebra
external mw_gene_linear_algebra_assymodel
external mw_select_algebra
integer mw_matrix_add_elem
integer mw_matrix_add_node
external mw_matrix_clear
external mw_vector_clear
external mw_assy_matrix_clear
external mw_assy_vector_clear

integer mw_matrix_add_elem_24
integer mw_matrix_add_elem_60
integer mw_matrix_add_elem_12
integer mw_matrix_add_elem_30
integer mw_matrix_add_elem_18
integer mw_matirx_add_elem_45
integer mw_matirx_add_elem_20
integer mw_matrix_add_elem_40
integer mw_matrix_add_elem_15
integer mw_matirx_add_elem_9
integer mw_matirx_add_elem_48
integer mw_matirx_add_elem_6
integer mw_matirx_add_elem_10
integer mw_matrix_rhs_set_bc2
integer mw_matrix_rhs_set_bc
integer mw_rhs_set_bc
integer mw_rhs_add_bc

integer mw_nl_rhs_set_bc
integer mw_nl_rhs_add_bc
!--
! solver
!--
integer mw_solve
!--
! solution_vector copy, at select MG-Level && select Equation
!--
external mw_get_solution_vector
external mw_get_solution_assy_vector
!--
! rhs_vector copy,  at select MG-Level && select Equation
!--
external mw_get_rhs_vector
external mw_get_rhs_assy_vector
external mw_get_rhs_load
external mw_get_rhs_assy_load
!--
! solution vector value, rhs vector value
!--
double precision mw_get_solution_assy_vector_val
double precision mw_get_rhs_assy_vector_val
!--
! solution vector dof, rhs vector dof
!--
integer mw_get_solution_assy_vector_dof
integer mw_get_rhs_assy_vector_dof
!--
! debug : assy_matrix display, rhs_assy_vector display
!--
external mw_dump_assy_matrix
external mw_dump_rhs_assy_vector

!--
! y = A*x : x update, y sumup
!--
external mw_matvec_assy
external mw_matvec

!----
! MG construct (refine)
!----
integer mw_refine
integer mw_mg_construct
!----
! model
!----
! assemble model
integer mw_get_num_of_assemble_model
external mw_select_assemble_model
!
! mesh part
!
integer mw_get_num_of_mesh_part
external mw_select_mesh_part_with_id
external mw_select_mesh_part
!
! mesh_part index && id
!
integer mw_get_mesh_part_id
integer mw_get_mesh_part_index
!
! element
!
external mw_select_element_with_id
external mw_select_element
integer mw_get_element_type
integer mw_get_num_of_element_vert
external mw_get_element_vert_node_id
external mw_get_element_vert_node_index
integer mw_get_num_of_element_edge
external mw_get_element_edge_node_id

integer mw_get_element_face_element_id
integer mw_get_num_Of_element_edge_element
integer mw_get_element_edge_element_id
external mw_setup_neighbors

!
! node
!
external mw_get_node_coord
integer mw_get_dof
integer mw_get_dof_scalar
integer mw_get_dof_vector
!
! node size, element size
!
integer mw_get_num_of_node
integer mw_get_num_of_node_with_mesh
integer mw_get_num_of_element
integer mw_get_num_of_element_with_mesh
!
! id && index
!
integer mw_get_element_id
integer mw_get_node_id
integer mw_get_node_index
integer mw_get_element_index
!
! parent node
!
integer mw_get_num_of_parent_node
integer mw_get_parent_node_id
!----
! node connectivity construct
!----
external mw_construct_node_connect_fem
external mw_get_node_connect_fem_size
external mw_get_node_connect_fem_item
!--
! node around element_id
!--
integer mw_get_num_of_aggregate_element
integer mw_get_aggregate_element_id
!----
! node type
!----
integer mw_nodetype_s
integer mw_nodetype_v
integer mw_nodetype_sv
!----
! element type
!----
integer mw_elemtype_hexa
integer mw_elemtype_hexa2
integer mw_elemtype_tetra
integer mw_elemtype_tetra2
integer mw_elemtype_prism
integer mw_elemtype_prism2
integer mw_elemtype_quad
integer mw_elemtype_quad2
integer mw_elemtype_triangle
integer mw_elemtype_triangle2
integer mw_elemtype_line
integer mw_elemtype_line2
!----
! frontISTR element type
!----
integer mw_fistr_elemtype_hexa
integer mw_fistr_elemtype_hexa2
integer mw_fistr_elemtype_tetra
integer mw_fistr_elemtype_tetra2
integer mw_fistr_elemtype_prism
integer mw_fistr_elemtype_prism2
integer mw_fistr_elemtype_quad
integer mw_fistr_elemtype_quad2
integer mw_fistr_elemtype_triangle
integer mw_fistr_elemtype_triangle2
integer mw_fistr_elemtype_line
integer mw_fistr_elemtype_line2
!----
! FrontISTR element type => MW3 element type
!----
integer mw_fistr_elemtype_to_mw3_elemtype
!----
! MW3 element type　=> FrontISTR element type
!----
integer mw_mw3_elemtype_to_fistr_elemtype
!----
! shape function
!----
integer mw_get_num_of_integereg_pointeger
external mw_shape_function_on_pt
external mw_shape_function_hexa81
external mw_shape_function_hexa82
external mw_shape_function_hexa201
external mw_shape_function_hexa202
external mw_shape_function_hexa203
external mw_shape_function_tetra41
external mw_shape_function_tetra101
external mw_shape_function_tetra104
external mw_shape_function_tetra1015
external mw_shape_function_prism62
external mw_shape_function_prism156
external mw_shape_function_prism159
external mw_shape_function_prism1518
external mw_shape_function_quad41
external mw_shape_function_quad84
external mw_shape_function_quad89
external mw_shape_function_tri31
external mw_shape_function_tri63
external mw_shape_function_line21
external mw_shape_function_line32
!----
! shape function deriv (rst coord)
!----
external mw_dndr
external mw_dndr_hexa81
external mw_dndr_hexa82
external mw_dndr_hexa201
external mw_dndr_hexa202
external mw_dndr_hexa203
external mw_dndr_tetra41
external mw_dndr_tetra101
external mw_dndr_tetra104
external mw_dndr_tetra1015
external mw_dndr_prism62
external mw_dndr_prism156
external mw_dndr_prism159
external mw_dndr_prism1518
external mw_dndr_quad41
external mw_dndr_quad84
external mw_dndr_quad89
external mw_dndr_tri31
external mw_dndr_tri63
external mw_dndr_line21
external mw_dndr_line32
!----
! shape function deriv (xyz coord)
!----
external mw_dndx
external mw_det_jacobian
external mw_weight
!----
! shape function type
!----
integer mw_shapetype_hexa81
integer mw_shapetype_hexa82
integer mw_shapetype_hexa201
integer mw_shapetype_hexa202
integer mw_shapetype_hexa203
integer mw_shapetype_tetra41
integer mw_shapetype_tetra101
integer mw_shapetype_tetra104
integer mw_shapetype_tetra1015
integer mw_shapetype_prism62
integer mw_shapetype_prism156
integer mw_shapetype_prism159
integer mw_shapetype_prism1518
integer mw_shapetype_quad41
integer mw_shapetype_quad84
integer mw_shapetype_quad89
integer mw_shapetype_tri31
integer mw_shapetype_tri63
integer mw_shapetype_line21
integer mw_shapetype_line32
!--
! boundary mesh
!--
! number of boudary_mesh
integer mw_get_num_of_boundary_bnode_mesh
integer mw_get_num_of_boundary_bface_mesh
integer mw_get_num_of_boundary_bedge_mesh
integer mw_get_num_of_boundary_bvolume_mesh
! BND type for each boundary_mesh { Neumann || Dirichlet }
integer mw_get_bnd_type_bnode_mesh
integer mw_get_bnd_type_bface_mesh
integer mw_get_bnd_type_bedge_mesh
integer mw_get_bnd_type_bvolume_mesh
! BND type number
integer mw_get_neumann_type
integer mw_get_dirichlet_type
! number of bnode for each boundary_mesh
integer mw_get_num_of_bnode_in_bnode_mesh
integer mw_get_num_of_bnode_in_bface_mesh
integer mw_get_num_of_bnode_in_bedge_mesh
integer mw_get_num_of_bnode_in_bvolume_mesh
! number of DOF for each boundary_mesh
integer mw_get_num_of_dof_in_bnode_mesh
integer mw_get_num_of_dof_in_bface_mesh
integer mw_get_num_of_dof_in_bedge_mesh
integer mw_get_num_of_dof_in_bvolume_mesh
! DOF number for each boundary_mesh ( DOF_index => DOF Number )
integer mw_get_dof_bnode_mesh
integer mw_get_dof_bface_mesh
integer mw_get_dof_bedge_mesh
integer mw_get_dof_bvolume_mesh
!--
! value of boundary node
!--
double precision mw_get_bnode_value_in_bnode_mesh
double precision mw_get_bnode_value_in_bface_mesh
double precision mw_get_bnode_value_in_bedge_mesh
double precision mw_get_bnode_value_in_bvolume_mesh
integer mw_get_node_id_in_bnode_mesh
integer mw_get_node_id_in_bface_mesh
integer mw_get_node_id_in_bedge_mesh
integer mw_get_node_id_in_bvolume_mesh
!--
! value of boundary face, edge, volume
!--
integer mw_get_num_of_bface
double precision mw_get_bface_value
integer mw_get_num_of_bedge
double precision mw_get_bedge_value
integer mw_get_num_of_bvolume
double precision mw_get_bvolume_value
!--
! node_id, face, edge, volume
!--
integer mw_get_num_of_node_bface
integer mw_get_node_id_bface
integer mw_get_num_of_node_bedge
integer mw_get_node_id_bedge
integer mw_get_num_of_node_bvolume
integer mw_get_node_id_bvolume
!--
! boundary_mesh name
!--
integer mw_get_bnode_mesh_namelength
external mw_get_bnode_mesh_name
integer mw_get_bface_mesh_namelength
external mw_get_bface_mesh_name
integer mw_get_bvolume_mesh_namelength
external mw_get_bvolume_mesh_name
integer mw_get_bedge_mesh_namelength
external mw_get_bedge_mesh_name
!--
! entity_id of boundary_mesh (for FrontISTR)
!--
integer mw_get_edge_id_bedge
integer mw_get_elem_id_bedge
integer mw_get_face_id_bface
integer mw_get_elem_id_bface
integer mw_get_elem_id_bvolume
!--
! mpi
!--
integer  mw_mpi_int
integer  mw_mpi_iint
integer  mw_mpi_uiint
integer  mw_mpi_double
integer  mw_mpi_comm
integer  mw_mpi_sum
integer  mw_mpi_max
integer  mw_mpi_min
integer mw_get_rank
integer mw_get_num_of_process
external mw_allreduce_r
integer mw_barrier
integer mw_abort
integer mw_allgather_r
integer mw_allgather_i
integer mw_gather_r
integer mw_gather_i
integer mw_scatter_r
integer mw_scatter_i
integer mw_bcast_i
integer mw_bcast_r
integer mw_bcast_s
!--
integer mw_get_num_of_neibpe
integer mw_get_transrank
!--
external mw_send_recv_r
external mw_send_recv_i

!--
! contact mesh
!--
integer mw_get_num_of_contact
integer mw_get_contact_id

!--
! Element_Group { after select AssyModel, select Mesh }
!--
integer mw_get_num_of_elementgroup
integer mw_get_num_of_element_id
integer mw_get_elementgroup_name_length
integer mw_get_element_id_with_elementgroup
!----
! logger
!----
external mw_logger_set_mode
external mw_logger_set_device
external mw_logger_info_mssg
external mw_logger_info
!----
! logger parameter
!----
integer mw_get_error_mode
integer mw_get_warn_mode
integer mw_get_info_mode
integer mw_get_debug_mode
integer mw_get_disk_device
integer mw_get_display_device
