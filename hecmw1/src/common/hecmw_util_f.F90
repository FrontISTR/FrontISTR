!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief I/O and Utility

module hecmw_util
#ifndef HECMW_SERIAL
  use mpi
#endif
  implicit none
  private :: hecmw_PETOT,hecmw_rank,hecmw_comm,hecmw_group
  public

  integer(kind=4),parameter:: kint  = 4
  integer(kind=4),parameter:: kreal = 8

  integer(kind=kint),parameter :: HECMW_NAME_LEN     =   63
  integer(kind=kint),parameter :: HECMW_HEADER_LEN   =  127
  integer(kind=kint),parameter :: HECMW_MSG_LEN      =  255
  integer(kind=kint),parameter :: HECMW_FILENAME_LEN = 1023

  integer(kind=kint),parameter :: hecmw_sum              = 46801
  integer(kind=kint),parameter :: hecmw_prod             = 46802
  integer(kind=kint),parameter :: hecmw_max              = 46803
  integer(kind=kint),parameter :: hecmw_min              = 46804
  integer(kind=kint),parameter :: hecmw_integer          = 53951
  integer(kind=kint),parameter :: hecmw_single_precision = 53952
  integer(kind=kint),parameter :: hecmw_double_precision = 53953
  integer(kind=kint),parameter :: hecmw_character        = 53954

#ifndef HECMW_SERIAL
  integer(kind=kint),parameter :: hecmw_status_size = MPI_STATUS_SIZE
#else
  integer(kind=kint),parameter :: hecmw_status_size = 1
#endif

  integer(kind=kint) :: hecmw_PETOT,hecmw_rank,hecmw_comm,hecmw_group
  !C
  !C +---------------+
  !C | SECTION info. |
  !C +---------------+
  !C===
  type hecmwST_section
    integer(kind=kint)         :: n_sect
    integer(kind=kint),pointer :: sect_type(:)
    integer(kind=kint),pointer :: sect_opt(:)
    integer(kind=kint),pointer :: sect_mat_ID_index(:)
    integer(kind=kint),pointer :: sect_mat_ID_item(:)
    integer(kind=kint),pointer :: sect_I_index(:)
    integer(kind=kint),pointer :: sect_I_item(:)
    integer(kind=kint),pointer :: sect_R_index(:)
    real(kind=kreal),pointer   :: sect_R_item(:)
    integer(kind=kint),pointer :: sect_orien_ID(:) => null()
  end type hecmwST_section

  !C      for hecmwST_section%sect_type
  integer(kind=kint),parameter :: HECMW_SECT_TYPE_SOLID     = 1
  integer(kind=kint),parameter :: HECMW_SECT_TYPE_SHELL     = 2
  integer(kind=kint),parameter :: HECMW_SECT_TYPE_BEAM      = 3
  integer(kind=kint),parameter :: HECMW_SECT_TYPE_INTERFACE = 4
  !C      for hecmwST_section%sect_opt
  integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRESS      =  0
  integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRAIN      =  1
  integer(kind=kint),parameter :: HECMW_SECT_OPT_ASYMMETRY    =  2
  integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRESS_RI   = 10
  integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRAIN_RI   = 11
  integer(kind=kint),parameter :: HECMW_SECT_OPT_ASYMMETRY_RI = 12
  !C===

  !C
  !C +----------------+
  !C | MATERIAL info. |
  !C +----------------+
  !C===
  type hecmwST_material
    integer(kind=kint)         :: n_mat
    integer(kind=kint)         :: n_mat_item
    integer(kind=kint)         :: n_mat_subitem
    integer(kind=kint)         :: n_mat_table
    character(HECMW_NAME_LEN),pointer :: mat_name(:)
    integer(kind=kint),pointer :: mat_item_index(:)
    integer(kind=kint),pointer :: mat_subitem_index(:)
    integer(kind=kint),pointer :: mat_table_index(:)
    real(kind=kreal),pointer   :: mat_val(:)
    real(kind=kreal),pointer   :: mat_temp(:)
  end type hecmwST_material
  !C===

  !C
  !C +-----------+
  !C | MPC info. |
  !C +-----------+
  !C===
  type hecmwST_mpc
    integer(kind=kint)         :: n_mpc
    integer(kind=kint),pointer :: mpc_index(:)
    integer(kind=kint),pointer :: mpc_item(:)
    integer(kind=kint),pointer :: mpc_dof(:)
    real(kind=kreal),pointer   :: mpc_val(:)
    real(kind=kreal),pointer   :: mpc_const(:)
  end type hecmwST_mpc
  !C===

  !C
  !C +-----------+
  !C | AMPLITUDE |
  !C +-----------+
  !C===
  type hecmwST_amplitude
    integer(kind=kint)         :: n_amp
    character(len=HECMW_NAME_LEN),pointer :: amp_name(:)
    integer(kind=kint),pointer :: amp_type_definition(:)
    integer(kind=kint),pointer :: amp_type_time(:)
    integer(kind=kint),pointer :: amp_type_value(:)
    integer(kind=kint),pointer :: amp_index(:)
    real(kind=kreal),pointer   :: amp_val(:)
    real(kind=kreal),pointer   :: amp_table(:)
  end type hecmwST_amplitude

  !C      for hecmwST_amplitude%amp_type_definition
  integer(kind=kint),parameter :: HECMW_AMP_TYPEDEF_TABULAR  = 1
  !C      for hecmwST_amplitude%amp_type_time
  integer(kind=kint),parameter :: HECMW_AMP_TYPETIME_STEP    = 1
  !C      for hecmwST_amplitude%amp_type_value
  integer(kind=kint),parameter :: HECMW_AMP_TYPEVAL_RELATIVE = 1
  integer(kind=kint),parameter :: HECMW_AMP_TYPEVAL_ABSOLUTE = 2
  !C===

  !C
  !C +-----------+
  !C | NODE grp. |
  !C +-----------+
  !C===
  type hecmwST_node_grp
    integer(kind=kint)         :: n_grp
    integer(kind=kint)         :: n_bc
    character(HECMW_NAME_LEN),pointer :: grp_name(:)
    integer(kind=kint),pointer :: grp_index(:)
    integer(kind=kint),pointer :: grp_item(:)
    integer(kind=kint),pointer :: bc_grp_ID(:)
    integer(kind=kint),pointer :: bc_grp_type(:)
    integer(kind=kint),pointer :: bc_grp_index(:)
    integer(kind=kint),pointer :: bc_grp_dof(:)
    real(kind=kreal),pointer   :: bc_grp_val(:)
  end type hecmwST_node_grp

  !C      for hecmwST_node_grp%bc_grp_type
  integer(kind=kint),parameter :: HECMW_BCGRPTYPE_DESPLACEMENT = 1
  integer(kind=kint),parameter :: HECMW_BCGRPTYPE_FLUX         = 2
  !C===

  !C
  !C +-----------+
  !C | ELEM grp. |
  !C +-----------+
  !C===
  type hecmwST_elem_grp
    integer(kind=kint)         :: n_grp
    integer(kind=kint)         :: n_bc
    character(HECMW_NAME_LEN),pointer :: grp_name(:)
    integer(kind=kint),pointer :: grp_index(:)
    integer(kind=kint),pointer :: grp_item(:)
    integer(kind=kint),pointer :: bc_grp_ID(:)
    integer(kind=kint),pointer :: bc_grp_type(:)
    integer(kind=kint),pointer :: bc_grp_index(:)
    real(kind=kreal),pointer   :: bc_grp_val(:)
  end type hecmwST_elem_grp

  !C      for hecmwST_elem_grp%bc_grp_type
  !C        integer(kind=kint),parameter :: HECMW_BCGRPTYPE_DESPLACEMENT = 1
  !C        integer(kind=kint),parameter :: HECMW_BCGRPTYPE_FLUX         = 2
  !C===

  !C
  !C +-----------+
  !C | SURF grp. |
  !C +-----------+
  !C===
  type hecmwST_surf_grp
    integer(kind=kint)         :: n_grp
    integer(kind=kint)         :: n_bc
    character(HECMW_NAME_LEN),pointer:: grp_name(:)
    integer(kind=kint),pointer :: grp_index(:)
    integer(kind=kint),pointer :: grp_item(:)
    integer(kind=kint),pointer :: bc_grp_ID(:)
    integer(kind=kint),pointer :: bc_grp_type(:)
    integer(kind=kint),pointer :: bc_grp_index(:)
    real(kind=kreal),pointer   :: bc_grp_val(:)
  end type hecmwST_surf_grp

  !C      for hecmwST_surf_grp%bc_grp_type
  !C        integer(kind=kint),parameter :: HECMW_BCGRPTYPE_DESPLACEMENT = 1
  !C        integer(kind=kint),parameter :: HECMW_BCGRPTYPE_FLUX         = 2
  !C

  !C
  !C +---------+
  !C | CONTACT |
  !C +---------+
  !C===
  type hecmwST_contact_pair
    integer(kind=kint)         :: n_pair
    character(HECMW_NAME_LEN),pointer:: name(:)
    integer(kind=kint),pointer :: type(:)
    integer(kind=kint),pointer :: slave_grp_id(:)
    integer(kind=kint),pointer :: master_grp_id(:)
    integer(kind=kint),pointer :: slave_orisgrp_id(:)
  end type hecmwST_contact_pair

  !C      for hecmwST_contact_pair%type
  integer(kind=kint),parameter :: HECMW_CONTACT_TYPE_NODE_SURF = 1
  integer(kind=kint),parameter :: HECMW_CONTACT_TYPE_SURF_SURF = 2
  integer(kind=kint),parameter :: HECMW_CONTACT_TYPE_NODE_ELEM = 3
  !C===

  !C
  !C +----------------+
  !C | REFINE Origin. |
  !C +----------------+
  !C===
  type hecmwST_refine_origin
    integer(kind=kint),pointer :: index(:)
    integer(kind=kint),pointer :: item_index(:)
    integer(kind=kint),pointer :: item_item(:)
  end type hecmwST_refine_origin
  !C===

  !C
  !C +------------------+
  !C | LOCAL MESH info. |
  !C +------------------+
  !C===
  type hecmwST_local_mesh

    !C
    !C-- FILES, GENERAL
    !C
    character(HECMW_FILENAME_LEN)         :: gridfile
    character(HECMW_FILENAME_LEN),pointer :: files(:)
    character(HECMW_HEADER_LEN)           :: header
    integer(kind=kint) :: hecmw_flag_adapt
    integer(kind=kint) :: hecmw_flag_initcon
    integer(kind=kint) :: hecmw_n_file
    integer(kind=kint) :: hecmw_flag_parttype
    integer(kind=kint) :: hecmw_flag_partdepth
    integer(kind=kint) :: hecmw_flag_version
    integer(kind=kint) :: hecmw_flag_partcontact
    real(kind=kreal)   :: zero_temp

    !C
    !C-- NODE
    integer(kind=kint)         :: n_node
    integer(kind=kint)         :: n_node_gross
    !         For parallel contact with hanging slave nodes
    integer(kind=kint)         :: nn_middle
    integer(kind=kint)         :: nn_internal
    integer(kind=kint)         :: n_dof
    integer(kind=kint)         :: n_dof_grp
    integer(kind=kint)         :: n_dof_tot
    real(kind=kreal),pointer   :: node(:)
    integer(kind=kint),pointer :: node_ID(:)
    integer(kind=kint),pointer :: global_node_ID(:)
    integer(kind=kint),pointer :: node_val_index(:)
    real(kind=kreal),pointer   :: node_val_item(:)
    integer(kind=kint),pointer :: node_dof_index(:)
    integer(kind=kint),pointer :: node_dof_item(:)
    integer(kind=kint),pointer :: node_init_val_index(:)
    real(kind=kreal),pointer   :: node_init_val_item(:)
    integer(kind=kint),pointer :: node_internal_list(:)
    !C
    !C-- ELEMENT
    !C
    integer(kind=kint)         :: n_elem
    integer(kind=kint)         :: n_elem_gross
    integer(kind=kint)         :: ne_internal
    integer(kind=kint)         :: n_elem_type
    integer(kind=kint)         :: n_elem_mat_ID
    integer(kind=kint),pointer :: elem_type_index(:)
    integer(kind=kint),pointer :: elem_type_item(:)
    integer(kind=kint),pointer :: elem_type(:)
    integer(kind=kint),pointer :: section_ID(:)
    integer(kind=kint),pointer :: elem_mat_ID_index(:)
    integer(kind=kint),pointer :: elem_mat_ID_item(:)
    integer(kind=kint),pointer :: elem_node_index(:)
    integer(kind=kint),pointer :: elem_node_item(:)
    integer(kind=kint),pointer :: elem_ID(:)
    integer(kind=kint),pointer :: global_elem_ID(:)
    integer(kind=kint),pointer :: elem_internal_list(:)
    integer(kind=kint),pointer :: elem_mat_int_index(:)
    real(kind=kreal),pointer   :: elem_mat_int_val(:)
    integer(kind=kint),pointer :: elem_val_index(:)
    real(kind=kreal),pointer   :: elem_val_item(:)
    !integer(kind=kint)         :: is_33shell
    !integer(kind=kint)         :: is_33beam
    !integer(kind=kint)         :: is_heat
    !C
    !C-- COMMUNICATION
    !C
    integer(kind=kint)         :: zero
    integer(kind=kint)         :: MPI_COMM
    integer(kind=kint)         :: PETOT
    integer(kind=kint)         :: PEsmpTOT
    integer(kind=kint)         :: my_rank
    integer(kind=kint)         :: errnof
    integer(kind=kint)         :: n_subdomain
    integer(kind=kint)         :: n_neighbor_pe
    integer(kind=kint),pointer :: neighbor_pe(:)
    integer(kind=kint),pointer :: import_index(:)
    integer(kind=kint),pointer :: import_item(:)
    integer(kind=kint),pointer :: export_index(:)
    integer(kind=kint),pointer :: export_item(:)
    integer(kind=kint),pointer :: shared_index(:)
    integer(kind=kint),pointer :: shared_item(:)

    !C
    !C-- ADAPTATION
    !C
    integer(kind=kint)         :: coarse_grid_level
    integer(kind=kint)         :: n_adapt
    integer(kind=kint),pointer :: when_i_was_refined_node(:)
    integer(kind=kint),pointer :: when_i_was_refined_elem(:)
    integer(kind=kint),pointer :: adapt_parent_type(:)
    integer(kind=kint),pointer :: adapt_type (:)
    integer(kind=kint),pointer :: adapt_level(:)
    integer(kind=kint),pointer :: adapt_parent(:)
    integer(kind=kint),pointer :: adapt_children_index(:)
    integer(kind=kint),pointer :: adapt_children_item(:)

    integer(kind=kint) :: nn_array, ne_array, nx_array
    integer(kind=kint) :: n_adapt_edge, n_adapt_edge_global
    integer(kind=kint) :: n_adapt_act_node, n_adapt_act_edge
    integer(kind=kint) :: n_adapt_act_elem, n_adapt_act_elem_cur
    integer(kind=kint) :: n_adapt_elem_341, n_adapt_elem_351
    integer(kind=kint) :: n_adapt_elem_341_cur, n_adapt_elem_351_cur
    integer(kind=kint) :: n_adapt_act_elem_341, n_adapt_act_elem_351

    integer(kind=kint) :: n_adapt_node_cur, nn_adapt_internal_cur
    integer(kind=kint) :: n_adapt_node_old, nn_adapt_internal_old
    integer(kind=kint) :: n_adapt_elem_cur, n_adapt_elem_old

    integer(kind=kint), pointer :: adapt_edge_node(:), adapt_mid_edge (:)
    integer(kind=kint), pointer :: adapt_iemb     (:), adapt_edge_home(:)
    integer(kind=kint), pointer :: adapt_act_edge (:)

    integer(kind=kint), pointer ::                                           &
      &           adapt_import_edge_index(:), adapt_import_edge_item (:),&
      &           adapt_export_edge_index(:), adapt_export_edge_item (:),&
      &           adapt_import_elem_index(:), adapt_import_elem_item (:),&
      &           adapt_export_elem_index(:), adapt_export_elem_item (:),&
      &           adapt_import_new_index (:), adapt_import_new_item  (:),&
      &           adapt_export_new_index (:), adapt_export_new_item  (:)
    integer(kind=kint), pointer :: rev_neighbor_pe(:)
    integer(kind=kint), pointer :: adapt_act_elem_341(:)
    integer(kind=kint), pointer :: adapt_act_elem_351(:)
    integer(kind=kint), pointer :: adapt_OLDtoNEW_node(:), adapt_NEWtoOLD_node(:)
    integer(kind=kint), pointer :: adapt_OLDtoNEW_elem(:), adapt_NEWtoOLD_elem(:)
    integer(kind=kint), pointer :: adapt_IWK(:), adapt_children_local(:)

    !C
    !C-- REFINEMENT
    !C
    integer(kind=kint)         :: n_refine
    integer(kind=kint),pointer :: node_old2new(:)
    integer(kind=kint),pointer :: node_new2old(:)
    integer(kind=kint),pointer :: elem_old2new(:)
    integer(kind=kint),pointer :: elem_new2old(:)
    integer(kind=kint),pointer :: n_node_refine_hist(:)

    !C
    !C-- ETC.
    !C
    type (hecmwST_section)   :: section
    type (hecmwST_material)  :: material
    type (hecmwST_mpc)       :: mpc
    type (hecmwST_amplitude) :: amp
    type (hecmwST_node_grp)  :: node_group
    type (hecmwST_elem_grp)  :: elem_group
    type (hecmwST_surf_grp)  :: surf_group
    type (hecmwST_contact_pair):: contact_pair
    type (hecmwST_refine_origin):: refine_origin

  end type hecmwST_local_mesh

  !C      for hecmwST_local_mesh%hecmw_flag_parttype
  integer(kind=kint),parameter :: HECMW_FLAG_PARTTYPE_UNKNOWN   = 0
  integer(kind=kint),parameter :: HECMW_FLAG_PARTTYPE_NODEBASED = 1
  integer(kind=kint),parameter :: HECMW_FLAG_PARTTYPE_ELEMBASED = 2

  !C      for hecmwST_local_mesh%hecmw_flag_partcontact
  integer(kind=kint),parameter :: HECMW_FLAG_PARTCONTACT_UNKNOWN    = 0
  integer(kind=kint),parameter :: HECMW_FLAG_PARTCONTACT_AGGREGATE  = 1
  integer(kind=kint),parameter :: HECMW_FLAG_PARTCONTACT_DISTRIBUTE = 2
  integer(kind=kint),parameter :: HECMW_FLAG_PARTCONTACT_SIMPLE     = 3

  !C
  !C +--------+
  !C | MATRIX |
  !C +--------+
  !C===
  type hecmwST_matrix_comm
    integer(kind=kint)                        :: zero
    integer(kind=kint)                        :: HECMW_COMM
    integer(kind=kint)                        :: PETOT
    integer(kind=kint)                        :: PEsmpTOT
    integer(kind=kint)                        :: my_rank
    integer(kind=kint)                        :: errnof
    integer(kind=kint)                        :: n_subdomain
    integer(kind=kint)                        :: n_neighbor_pe
    integer(kind=kint), dimension(:), pointer :: neighbor_pe
    integer(kind=kint), dimension(:), pointer :: import_index
    integer(kind=kint), dimension(:), pointer :: import_item
    integer(kind=kint), dimension(:), pointer :: export_index
    integer(kind=kint), dimension(:), pointer :: export_item
    integer(kind=kint), dimension(:), pointer :: shared_index
    integer(kind=kint), dimension(:), pointer :: shared_item
  end type hecmwST_matrix_comm

  type hecmwST_index_value_pair
    integer(kind=kint) :: i
    integer(kind=kint) :: j
    real(kind=kreal), dimension(3,3) :: val
  end type hecmwST_index_value_pair

  !> Structure for Lagrange multiplier-related part of stiffness matrix
  !> (Lagrange multiplier-related matrix)
  type hecmwST_matrix_lagrange
    integer(kind=kint) :: num_lagrange = 0 !< total number of Lagrange multipliers
    integer(kind=kint) :: numL_lagrange = 0 !< node-based number of non-zero items in lower triangular half of matrix
    integer(kind=kint) :: numU_lagrange = 0 !< node-based number of non-zero items in upper triangular half of matrix

    integer(kind=kint), pointer  :: indexL_lagrange(:) => null() !< node-based index of first non-zero item of each row in lower part
    integer(kind=kint), pointer  :: indexU_lagrange(:) => null() !< node-based index of first non-zero item of each row in upper part

    integer(kind=kint), pointer  :: itemL_lagrange(:) => null() !< node-based column number of non-zero items in lower part
    integer(kind=kint), pointer  :: itemU_lagrange(:) => null() !< node-based column number of non-zero items in upper part

    real(kind=kreal),    pointer  :: AL_lagrange(:) => null() !< values of non-zero items in lower part
    real(kind=kreal),    pointer  :: AU_lagrange(:) => null() !< values of non-zero items in upper part

    real(kind=kreal),    pointer  :: Lagrange(:) => null() !< values of Lagrange multipliers
  end type hecmwST_matrix_lagrange

  type hecmwST_matrix
    integer(kind=kint) ::  N, NP, NPL, NPU, NDOF
    real(kind=kreal), pointer :: D(:), B(:), X(:), ALU(:)
    real(kind=kreal), pointer :: AL(:), AU(:)
    integer(kind=kint), pointer :: indexL(:), indexU(:)
    integer(kind=kint), pointer ::  itemL(:),  itemU(:)
    integer(kind=kint ), dimension(100) :: Iarray
    real   (kind=kreal), dimension(100) :: Rarray
    logical :: symmetric = .true.
  end type hecmwST_matrix
contains

  !C
  !C***
  !C*** HECMW_INIT
  !C***
  !C
  !C    INIT. HECMW-FEM process's
  !C
  subroutine hecmw_init
    character(len=HECMW_FILENAME_LEN):: ctrlfile = "hecmw_ctrl.dat"
    call hecmw_init_ex(ctrlfile)
  end subroutine hecmw_init


  !C
  !C***
  !C*** HECMW_INIT_EX
  !C***
  !C
  !C    INIT. HECMW-FEM process's
  !C
  subroutine hecmw_init_ex(ctrlfile)
    character(len=HECMW_FILENAME_LEN):: ctrlfile
    integer(kind=kint) :: ierr

#ifndef HECMW_SERIAL
    !call MPI_INIT (ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, hecmw_PETOT, ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, hecmw_rank,  ierr)
    call MPI_COMM_DUP  (MPI_COMM_WORLD, hecmw_comm,  ierr)
    call MPI_COMM_GROUP(MPI_COMM_WORLD, hecmw_group, ierr)
#else
    hecmw_PETOT=1
    hecmw_rank=0
    hecmw_comm=0
    hecmw_group=0
    ierr=0
#endif

    call hecmw_comm_init_if(hecmw_comm, hecmw_PETOT, hecmw_rank, hecmw_group)

    call hecmw_ctrl_init_ex_if(ctrlfile, ierr)
    if(ierr /= 0) then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif
    !      call hecmw_couple_comm_init_if(ierr)
    !      if(ierr /= 0) then
    !          call  hecmw_abort( hecmw_comm_get_comm( ) )
    !      endif

  end subroutine hecmw_init_ex


  !C
  !C***
  !C*** HECMW_FINALIZE
  !C***
  !C
  !C    FINALIZE. HECMW-FEM process's
  !C
  subroutine hecmw_finalize
    integer(kind=kint) :: ierr

    call hecmw_ctrl_finalize_if()

#ifndef HECMW_SERIAL
    call MPI_FINALIZE(ierr)
#endif

  end subroutine hecmw_finalize


  !C******************** MPI WRAPPER SUBROUTINES ************************
  !C
  !C***
  !C*** HECMW_ABORT
  !C***
  !C
  subroutine hecmw_abort(comm)
    integer(kind=kint) :: comm, errorcode, ierror

#ifndef HECMW_SERIAL
    call MPI_ABORT(comm, errorcode, ierror)
#else
    stop
#endif
  end subroutine hecmw_abort

  !C
  !C***
  !C*** HECMW_WTIME
  !C***
  !C
  function hecmw_wtime()
    real(kind=kreal) hecmw_wtime
    external  hecmw_Wtime_fi
    real(kind=kreal) hecmw_Wtime_fi
    hecmw_Wtime = hecmw_Wtime_fi()
  end function hecmw_Wtime
  !C
  !C***
  !C*** HECMW_WTICK
  !C***
  !C
  function hecmw_wtick()
    real(kind=kreal) hecmw_wtick
    external hecmw_Wtick_fi
    real(kind=kreal) hecmw_Wtick_fi
    hecmw_Wtick = hecmw_Wtick_fi()
  end function hecmw_Wtick
  !C
  !C***
  !C*** HECMW_COMM_GET_COMM
  !C***
  !C
  function hecmw_comm_get_comm() result(comm)
    integer(kind=kint) :: comm

    comm = hecmw_comm
  end function hecmw_comm_get_comm

  !C
  !C***
  !C*** HECMW_COMM_GET_RANK
  !C***
  !C
  function hecmw_comm_get_rank() result(rank)
    integer(kind=kint) :: rank

    rank = hecmw_rank
  end function hecmw_comm_get_rank

  !C
  !C***
  !C*** HECMW_COMM_GET_SIZE
  !C***
  !C
  function hecmw_comm_get_size() result(comm_size)
    integer(kind=kint) :: comm_size

    comm_size = hecmw_PETOT
  end function hecmw_comm_get_size


  !C***************  NULL POINTER SETTING UTILITY ****************

  subroutine hecmw_nullify_section( P )
    type( hecmwST_section ) :: P
    nullify( P%sect_type )
    nullify( P%sect_opt )
    nullify( P%sect_mat_ID_index )
    nullify( P%sect_mat_ID_item )
    nullify( P%sect_I_index )
    nullify( P%sect_I_item )
    nullify( P%sect_R_index )
    nullify( P%sect_R_item )
  end subroutine hecmw_nullify_section

  subroutine hecmw_nullify_material( P )
    type( hecmwST_material ) :: P
    nullify( P%mat_name )
    nullify( P%mat_item_index )
    nullify( P%mat_subitem_index )
    nullify( P%mat_table_index )
    nullify( P%mat_val )
    nullify( P%mat_temp )
  end subroutine hecmw_nullify_material

  subroutine hecmw_nullify_mpc( P )
    type( hecmwST_mpc ) :: P
    nullify( P%mpc_index )
    nullify( P%mpc_item )
    nullify( P%mpc_dof )
    nullify( P%mpc_val )
    nullify( P%mpc_const )
  end subroutine hecmw_nullify_mpc

  subroutine hecmw_initialize_mpc( mpc, n_mpc, n_item )
    type( hecmwST_mpc ), intent(inout) :: mpc
    integer(kind=kint), intent(in)    :: n_mpc
    integer(kind=kint), intent(in)    :: n_item

    mpc%n_mpc = n_mpc
    allocate( mpc%mpc_index(0:n_mpc) )
    allocate( mpc%mpc_item(n_item) )
    allocate( mpc%mpc_dof(n_item) )
    allocate( mpc%mpc_val(n_item) )
  end subroutine

  subroutine hecmw_finalize_mpc( P )
    type( hecmwST_mpc ) :: P
    if( associated(P%mpc_index) ) deallocate( P%mpc_index )
    if( associated(P%mpc_item) )  deallocate( P%mpc_item )
    if( associated(P%mpc_dof) )   deallocate( P%mpc_dof )
    if( associated(P%mpc_val) )   deallocate( P%mpc_val )
  end subroutine hecmw_finalize_mpc

  subroutine hecmw_nullify_amplitude( P )
    type( hecmwST_amplitude ) :: P
    nullify( P%amp_name )
    nullify( P%amp_type_definition )
    nullify( P%amp_type_time )
    nullify( P%amp_type_value )
    nullify( P%amp_index )
    nullify( P%amp_val )
    nullify( P%amp_table )
  end subroutine hecmw_nullify_amplitude

  subroutine hecmw_nullify_node_grp( P )
    type( hecmwST_node_grp ) :: P
    nullify( P%grp_name )
    nullify( P%grp_index )
    nullify( P%grp_item )
    nullify( P%bc_grp_ID )
    nullify( P%bc_grp_type )
    nullify( P%bc_grp_index )
    nullify( P%bc_grp_dof )
    nullify( P%bc_grp_val )
  end subroutine hecmw_nullify_node_grp

  subroutine hecmw_nullify_elem_grp( P )
    type( hecmwST_elem_grp ) :: P
    nullify( P%grp_name )
    nullify( P%grp_index )
    nullify( P%grp_item )
    nullify( P%bc_grp_ID )
    nullify( P%bc_grp_type )
    nullify( P%bc_grp_index )
    nullify( P%bc_grp_val )
  end subroutine hecmw_nullify_elem_grp

  subroutine hecmw_nullify_surf_grp( P )
    type( hecmwST_surf_grp ) :: P
    nullify( P%grp_name )
    nullify( P%grp_index )
    nullify( P%grp_item )
    nullify( P%bc_grp_ID )
    nullify( P%bc_grp_type )
    nullify( P%bc_grp_index )
    nullify( P%bc_grp_val )
  end subroutine hecmw_nullify_surf_grp

  subroutine hecmw_nullify_contact_pair( P )
    type( hecmwST_contact_pair ) :: P
    nullify( P%name )
    nullify( P%type )
    nullify( P%slave_grp_id )
    nullify( P%slave_orisgrp_id )
    nullify( P%master_grp_id )
  end subroutine hecmw_nullify_contact_pair

  subroutine hecmw_nullify_refine_origin( P )
    type( hecmwST_refine_origin ) :: P
    nullify( P%index )
    nullify( P%item_index )
    nullify( P%item_item )
  end subroutine hecmw_nullify_refine_origin

  subroutine hecmw_nullify_mesh( P )
    type( hecmwST_local_mesh ) :: P
    nullify( P%files )
    nullify( P%node )
    nullify( P%node_ID )
    nullify( P%global_node_ID )
    nullify( P%node_val_index )
    nullify( P%node_val_item )
    nullify( P%node_dof_index )
    nullify( P%node_dof_item )
    nullify( P%node_init_val_index )
    nullify( P%node_init_val_item )
    nullify( P%node_internal_list )
    nullify( P%elem_type_index )
    nullify( P%elem_type_item )
    nullify( P%elem_type )
    nullify( P%section_ID )
    nullify( P%elem_mat_ID_index )
    nullify( P%elem_mat_ID_item )
    nullify( P%elem_node_index )
    nullify( P%elem_node_item )
    nullify( P%elem_ID )
    nullify( P%global_elem_ID )
    nullify( P%elem_internal_list )
    nullify( P%elem_mat_int_index )
    nullify( P%elem_mat_int_val )
    nullify( P%elem_val_index )
    nullify( P%elem_val_item )
    nullify( P%neighbor_pe )
    nullify( P%import_index )
    nullify( P%import_item )
    nullify( P%export_index )
    nullify( P%export_item )
    nullify( P%shared_index )
    nullify( P%shared_item )
    nullify( P%when_i_was_refined_node )
    nullify( P%when_i_was_refined_elem )
    nullify( P%adapt_parent_type )
    nullify( P%adapt_type  )
    nullify( P%adapt_level )
    nullify( P%adapt_parent )
    nullify( P%adapt_children_index )
    nullify( P%adapt_children_item )
    nullify( P%adapt_edge_node )
    nullify( P%adapt_mid_edge  )
    nullify( P%adapt_iemb      )
    nullify( P%adapt_edge_home )
    nullify( P%adapt_act_edge  )
    nullify( P%adapt_import_edge_index )
    nullify( P%adapt_import_edge_item  )
    nullify( P%adapt_export_edge_index )
    nullify( P%adapt_export_edge_item  )
    nullify( P%adapt_import_elem_index )
    nullify( P%adapt_import_elem_item  )
    nullify( P%adapt_export_elem_index )
    nullify( P%adapt_export_elem_item  )
    nullify( P%adapt_import_new_index  )
    nullify( P%adapt_import_new_item   )
    nullify( P%adapt_export_new_index  )
    nullify( P%adapt_export_new_item   )
    nullify( P%rev_neighbor_pe )
    nullify( P%adapt_act_elem_341 )
    nullify( P%adapt_act_elem_351 )
    nullify( P%adapt_OLDtoNEW_node )
    nullify( P%adapt_NEWtoOLD_node )
    nullify( P%adapt_OLDtoNEW_elem )
    nullify( P%adapt_NEWtoOLD_elem )
    nullify( P%adapt_IWK )
    nullify( P%adapt_children_local )
    nullify( P%node_old2new )
    nullify( P%node_new2old )
    nullify( P%elem_old2new )
    nullify( P%elem_new2old )
    nullify( P%n_node_refine_hist )

    call hecmw_nullify_section( P%section )
    call hecmw_nullify_material( P%material )
    call hecmw_nullify_mpc( P%mpc )
    call hecmw_nullify_amplitude( P%amp )
    call hecmw_nullify_node_grp( P%node_group )
    call hecmw_nullify_elem_grp( P%elem_group )
    call hecmw_nullify_surf_grp( P%surf_group )
    call hecmw_nullify_contact_pair( P%contact_pair )
    call hecmw_nullify_refine_origin( P%refine_origin )

  end subroutine hecmw_nullify_mesh


  subroutine hecmw_nullify_matrix_comm( P )
    type( hecmwST_matrix_comm ) :: P
    nullify( P%neighbor_pe )
    nullify( P%import_index )
    nullify( P%import_item )
    nullify( P%export_index )
    nullify( P%export_item )
    nullify( P%shared_index )
    nullify( P%shared_item )
  end subroutine hecmw_nullify_matrix_comm

  subroutine hecmw_nullify_matrix( P )
    type( hecmwST_matrix ) :: P
    nullify( P%D )
    nullify( P%B )
    nullify( P%X )
    nullify( P%ALU )
    nullify( P%AL )
    nullify( P%AU )
    nullify( P%indexL )
    nullify( P%indexU )
    nullify( P%itemL )
    nullify( P%itemU )
  end subroutine hecmw_nullify_matrix

  subroutine hecmw_print_matrix( fname, P )
    character(len=*), intent(in)       :: fname
    type( hecmwST_matrix ), intent(in) :: P

    integer :: i, nf, nBlock
    nf = 777
    nBlock = P%NDOF * P%NDOF
    open( unit=nf, file=fname)
    write( nf, * ) P%N,P%NP,P%NPL,P%NPU,P%NDOF
    !---- index
    do i=0, P%NP
      write( nf,* )  P%indexL(i), P%indexU(i)
    enddo
    !---- itemL, AL
    do i=1, P%NPL
      write( nf,* ) P%itemL(i)
    enddo
    do i=1,nBlock * P%NPL
      write( nf,* ) P%AL(i)
    enddo
    !---- itemU, AU
    do i=1,P%NPU
      write( nf,* ) P%itemU(i)
    enddo
    do i=1,nBlock * P%NPU
      write( nf,* ) P%AU(i)
    enddo
    !---- D
    do i=1,nBlock * P%NP
      write( nf, * ) P%D(i)
    enddo
    !---- B
    do i=1,P%NDOF * P%NP
      write( nf, * ) P%B(i)
    enddo

    close( nf )
  end subroutine

  subroutine hecmw_read_matrix( fname, P )
    character(len=*), intent(in)        :: fname
    type( hecmwST_matrix ), intent(out) :: P

    integer :: i, nf, nBlock, istat
    nf = 777
    open( unit=nf, file=fname, status='old', iostat= istat)
    if(istat /= 0) then
      print *, "cannot open file ",fname
      stop
    endif
    read( nf, * ) P%N,P%NP,P%NPL,P%NPU,P%NDOF
    nBlock = P%NDOF * P%NDOF

    !----It is supposing of array are not allocated yet
    allocate( P%indexL(0:P%NP), P%indexU(0:P%NP) )
    allocate( P%itemL(P%NPL), P%itemU(P%NPU) )
    allocate( P%AL(P%NPL*nBlock), P%AU(P%NPU*nBlock) )
    allocate( P%D(P%NP*nBlock) )
    allocate( P%B(P%NDOF*P%NP) )
    !---- index
    do i=0, P%NP
      read( nf,* )  P%indexL(i), P%indexU(i)
    enddo
    !---- itemL, AL
    do i=1, P%NPL
      read( nf,* ) P%itemL(i)
    enddo
    do i=1,nBlock * P%NPL
      read( nf,* ) P%AL(i)
    enddo
    !---- itemU, AU
    do i=1,P%NPU
      read( nf,* ) P%itemU(i)
    enddo
    do i=1,nBlock * P%NPU
      read( nf,* ) P%AU(i)
    enddo
    !---- D
    do i=1,nBlock * P%NP
      read( nf, * ) P%D(i)
    enddo
    !---- B
    do i=1,P%NDOF * P%NP
      read( nf, * ) P%B(i)
    enddo

    close( nf )
  end subroutine

  subroutine hecmw_clone_matrix(hecMATorig,hecMAT)
    type (hecmwST_matrix    ) :: hecMATorig
    type (hecmwST_matrix    ),pointer :: hecMAT
    allocate(hecMAT)
    call hecmw_nullify_matrix( hecMAT )

    hecMAT%B => hecMATorig%B
    hecMAT%X => hecMATorig%X
    hecMAT%D => hecMATorig%D
    hecMAT%AL => hecMATorig%AL
    hecMAT%AU => hecMATorig%AU
    hecMAT%indexL => hecMATorig%indexL
    hecMAT%indexU => hecMATorig%indexU
    hecMAT%itemL  => hecMATorig%itemL
    hecMAT%itemU  => hecMATorig%itemU
    hecMAT%N    = hecMATorig%N
    hecMAT%NP   = hecMATorig%NP
    hecMAT%NPL  = hecMATorig%NPL
    hecMAT%NPU  = hecMATorig%NPU
    hecMAT%NDOF = hecMATorig%NDOF
    hecMAT%Iarray = hecMATorig%Iarray
    hecMAT%Rarray = hecMATorig%Rarray
  end subroutine hecmw_clone_matrix

  subroutine hecmw_copy_matrix(hecMATorig,hecMAT)
    type (hecmwST_matrix    ) :: hecMATorig
    type (hecmwST_matrix    ),pointer :: hecMAT
    integer(kind=kint) NDOF,NDOF2,N,NP,NPL,NPU

    allocate(hecMAT)
    call hecmw_nullify_matrix( hecMAT )
    N = hecMATorig%N
    NDOF = hecMATorig%NDOF
    NDOF2 = NDOF*NDOF
    NP = hecMATorig%NP
    NPL = hecMATorig%NPL
    NPU = hecMATorig%NPU
    hecMAT%N    = N
    hecMAT%NP   = NP
    hecMAT%NPL  = NPL
    hecMAT%NPU  = NPU
    hecMAT%NDOF = NDOF
    allocate(hecMAT%B(NDOF*NP))
    allocate(hecMAT%X(NDOF*NP))
    allocate(hecMAT%D(NDOF2*NP))
    allocate(hecMAT%AU(NDOF2*NPU))
    allocate(hecMAT%AL(NDOF2*NPL))
    allocate(hecMAT%indexL(0:NP), hecMAT%indexU(0:NP), hecMAT%itemL(NPL), hecMAT%itemU(NPU))
    hecMAT%B       = hecMATorig%B
    hecMAT%X       = hecMATorig%X
    hecMAT%D       = hecMATorig%D
    hecMAT%AU      = hecMATorig%AU
    hecMAT%AL      = hecMATorig%AL
    hecMAT%indexL  = hecMATorig%indexL
    hecMAT%indexU  = hecMATorig%indexU
    hecMAT%itemL   = hecMATorig%itemL
    hecMAT%itemU   = hecMATorig%itemU
    hecMAT%Iarray  = hecMATorig%Iarray
    hecMAT%Rarray  = hecMATorig%Rarray
  end subroutine hecmw_copy_matrix

  subroutine hecmw_blockmatrix_expand(hecMATorig,hecMAT,NDOF)
    type (hecmwST_matrix    ) :: hecMATorig
    type (hecmwST_matrix    ),pointer :: hecMAT
    integer(kind=kint) NDOF,NDOF2,oNDOF,oNDOF2,i,j,k

    NDOF2 = NDOF*NDOF
    oNDOF = hecMATorig%NDOF
    oNDOF2 = oNDOF*oNDOF
    allocate(hecMAT)
    call hecmw_nullify_matrix( hecMAT )

    allocate(hecMAT%B(NDOF*hecMATorig%NP))
    allocate(hecMAT%X(NDOF*hecMATorig%NP))
    allocate(hecMAT%D(NDOF2*hecMATorig%NP))
    allocate(hecMAT%AL(NDOF2*hecMATorig%NPL))
    allocate(hecMAT%AU(NDOF2*hecMATorig%NPU))
    hecMAT%indexL => hecMATorig%indexL
    hecMAT%indexU => hecMATorig%indexU
    hecMAT%itemL  => hecMATorig%itemL
    hecMAT%itemU  => hecMATorig%itemU
    hecMAT%N    = hecMATorig%N
    hecMAT%NP   = hecMATorig%NP
    hecMAT%NPL  = hecMATorig%NPL
    hecMAT%NPU  = hecMATorig%NPU
    hecMAT%NDOF = NDOF
    hecMAT%Iarray = hecMATorig%Iarray
    hecMAT%Rarray = hecMATorig%Rarray
    hecMAT%X = 0.0d0
    do i = 1, hecMATorig%NP
      do j = 1, NDOF
        do k = 1, NDOF
          if (j<=oNDOF .and. k<=oNDOF) then
            hecMAT%D(NDOF2*(i-1) + (j-1)*NDOF + k) = hecMATorig%D(oNDOF2*(i-1) + (j-1)*oNDOF + k)
          else
            if(j==k) then
              hecMAT%D(NDOF2*(i-1) + (j-1)*NDOF + k)=1
            else
              hecMAT%D(NDOF2*(i-1) + (j-1)*NDOF + k)=0
            end if
          end if
        end do
        if (j<=oNDOF) then
          hecMAT%B(NDOF*(i-1) + j) = hecMATorig%B(oNDOF*(i-1) + j)
        else
          hecMAT%B(NDOF*(i-1) + j)=0
        end if
      end do
    end do
    do i = 1, hecMATorig%NPL
      do j = 1, NDOF
        do k = 1, NDOF
          if (j<=oNDOF .and. k<=oNDOF) then
            hecMAT%AL(NDOF2*(i-1) + (j-1)*NDOF + k) = hecMATorig%AL(oNDOF2*(i-1) + (j-1)*oNDOF + k)
          else
            hecMAT%AL(NDOF2*(i-1) + (j-1)*NDOF + k) = 0
          end if
        end do
      end do
    end do
    do i = 1, hecMATorig%NPU
      do j = 1, NDOF
        do k = 1, NDOF
          if (j<=oNDOF .and. k<=oNDOF) then
            hecMAT%AU(NDOF2*(i-1) + (j-1)*NDOF + k) = hecMATorig%AU(oNDOF2*(i-1) + (j-1)*oNDOF + k)
          else
            hecMAT%AU(NDOF2*(i-1) + (j-1)*NDOF + k) = 0
          end if
        end do
      end do
    end do
  end subroutine hecmw_blockmatrix_expand
  subroutine hecmw_vector_contract(hecMATorig,hecMAT,NDOF)
    type (hecmwST_matrix    ) :: hecMATorig
    type (hecmwST_matrix    ),pointer :: hecMAT
    integer(kind=kint) NDOF,NDOF2,oNDOF,i,j
    NDOF2 = NDOF*NDOF
    oNDOF = hecMATorig%NDOF
    do i = 1, hecMATorig%NP
      do j = 1, oNDOF
        hecMATorig%X(oNDOF*(i-1) + j) = hecMAT%X(NDOF*(i-1) + j)
      end do
    end do
    hecMATorig%Iarray = hecMAT%Iarray
    hecMATorig%Rarray = hecMAT%Rarray
    deallocate(hecMAT%B)
    deallocate(hecMAT%D)
    deallocate(hecMAT%X)
    deallocate(hecMAT%AL)
    deallocate(hecMAT%AU)
    deallocate(hecMAT)
  end subroutine hecmw_vector_contract
end module hecmw_util
