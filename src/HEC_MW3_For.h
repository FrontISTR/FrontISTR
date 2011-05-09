!     
! File:   HEC_MW3_For.f
! Author: ktakeda
!
! Created on 2010/10/04, 15:25
!

!----
! HEC_MW3 construct & destruct
!----
integer mw_initialize_1
integer mw_initialize_2
integer mw_finalize

!----
! file i/o API
!----
integer mw_file_read
integer mw_file_write

!----
! linear solver API
!----
integer mw_initialize_matrix
integer mw_initialize_vector

integer mw_matrix_add_elem

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

integer mw_matrix_set_bc
integer mw_rhs_set_bc

integer mw_solve

external mw_store_matrix
external mw_load_matrix


!----
! MG construct (refine)
!----
integereger mw_refine
integer mw_mg_construct


!----
! model
!----
! assemble model
integer mw_get_num_of_assemble_model
external mw_select_assemble_model
! mesh part
integer mw_get_num_of_mesh_part
external mw_select_mesh_part_with_id
external mw_select_mesh_part
! element
external mw_select_element_with_id
external mw_select_element
integer mw_get_element_type
integer mw_get_num_of_element_vert

external mw_get_element_vert_node_id
integer mw_get_num_of_element_edge
external mw_get_element_edge_node_id

! node
external mw_get_node_coord
integer mw_get_dof
integer mw_get_dof_scalar
integer mw_get_dof_vector

external mw_set_node_value
external mw_set_node_value_with_dof
external mw_get_node_value
external mw_get_node_value_with_dof

external mw_set_sv_node_value
external mw_set_sv_node_value_with_dof
external mw_get_sv_node_value
external mw_get_sv_node_value_with_dof

! node size, element size
integer mw_get_num_of_node
integer mw_get_num_of_node_with_mesh
integer mw_get_num_of_element
integer mw_get_num_of_element_with_mesh

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
integer mw_get_num_of_boundary_bnode_mesh
integer mw_get_num_of_boundary_bface_mesh
integer mw_get_num_of_boundary_bedge_mesh
integer mw_get_num_of_boundary_bvolume_mesh
integer mw_get_num_of_bnode_in_bnode_mesh
integer mw_get_num_of_bnode_in_bface_mesh
integer mw_get_num_of_bnode_in_bedge_mesh
integer mw_get_num_of_bnode_in_bvolume_mesh
integer mw_get_num_of_dof_in_bnode_mesh
integer mw_get_num_of_dof_in_bface_mesh
integer mw_get_num_of_dof_in_bedge_mesh
integer mw_get_num_of_dof_in_bvolume_mesh
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
! mpi
!--
integer  mw_mpi_sum
integer  mw_mpi_max
integer  mw_mpi_min
external mw_allreduce_r
external mw_send_recv_r2
external mw_send_recv_r

!----
! logger
!----
external mw_logger_set_mode
external mw_logger_set_device
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



