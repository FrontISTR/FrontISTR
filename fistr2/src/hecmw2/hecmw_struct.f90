!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.3                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Noboru Imai (Univ. of Tokyo)                   !
!                       Kazuaki Sakane (RIST)                          !
!                       Shin'ichi Ezure (RIST)                         !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!> \brief  This module provides basic data structure definition 
Module hecmw

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

        integer(kind=kint) :: hecmw_PETOT,hecmw_rank,hecmw_comm,hecmw_group

!C
!C +---------------+
!C | SECTION info. |
!C +---------------+
!C===
        type tSection
          character(HECMW_NAME_LEN)  :: mat_name
          character(HECMW_NAME_LEN)  :: egroup_name
          integer(kind=kint)         :: sect_type      ! 0: solid  1: shell
          integer                    :: sect_opt
          integer(kind=kint)         :: sect_mat_ID
      !    integer(kind=kint),pointer :: sect_I_item(:) => null()
          real(kind=kreal)           :: sect_R_item(4)
        end type tSection
		
        type( tSection ), pointer :: MWSections(:) => null()
		
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
        end type hecmwST_section

!C      for hecmwST_section%sect_type
        integer(kind=kint),parameter :: HECMW_SECT_TYPE_SOLID     = 1
        integer(kind=kint),parameter :: HECMW_SECT_TYPE_SHELL     = 2
        integer(kind=kint),parameter :: HECMW_SECT_TYPE_BEAM      = 3
        integer(kind=kint),parameter :: HECMW_SECT_TYPE_INTERFACE = 4
!C      for hecmwST_section%sect_opt( 2-dimensional problems)
        integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRESS      =  0
        integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRAIN      =  1
        integer(kind=kint),parameter :: HECMW_SECT_OPT_ASYMMETRY    =  2
        integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRESS_RI   = 10
        integer(kind=kint),parameter :: HECMW_SECT_OPT_PSTRAIN_RI   = 11
        integer(kind=kint),parameter :: HECMW_SECT_OPT_ASYMMETRY_RI = 12
!C      for hecmwST_section%sect_opt( solid element )
        integer(kind=kint),parameter :: HECMW_SECT_OPT_FULL         = 0
        integer(kind=kint),parameter :: HECMW_SECT_OPT_SRI          = 21
        integer(kind=kint),parameter :: HECMW_SECT_OPT_URI          = 22
        integer(kind=kint),parameter :: HECMW_SECT_OPT_IC           = 23
!C      for hecmwST_section%sect_opt( shell element )
        integer(kind=kint),parameter :: HECMW_SECT_OPT_MITC         = 0

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
        type tAmplitude
          character(len=HECMW_NAME_LEN) :: apname
          integer(kind=kint)            :: itype
          integer(kind=kint)            :: type_time
          integer(kind=kint)            :: type_value
          real(kind=kreal),pointer      :: table(:,:)
        end type tAmplitude
		
        type( tAmplitude ), pointer :: MWAmplitudes(:) => null()
		
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

        type fstr_boundary_grp
          integer                    :: gid
          integer                    :: part_id
          character(HECMW_NAME_LEN)  :: grp_name
          character(HECMW_NAME_LEN)  :: amp
          integer                    :: dof(6)
          real(kind=kreal)           :: fval
        end type
		
        type fstr_dload_grp
          integer                    :: gid
          integer                    :: part_id
          character(HECMW_NAME_LEN)  :: grp_name
          character(HECMW_NAME_LEN)  :: amp
          integer                    :: itype
          real(kind=kreal)           :: fval(7)
        end type
		
        type fstr_ndscalar_grp
          integer                    :: gid
          integer                    :: part_id
          character(HECMW_NAME_LEN)  :: grp_name
          character(HECMW_NAME_LEN)  :: amp
          character(HECMW_NAME_LEN)  :: fname
          real(kind=kreal)           :: fval
        end type
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
        end type hecmwST_contact_pair

!C      for hecmwST_contact_pair%type
        integer(kind=kint),parameter :: HECMW_CONTACT_TYPE_NODE_SURF = 1
        integer(kind=kint),parameter :: HECMW_CONTACT_TYPE_SURF_SURF = 2
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
          real(kind=kreal)   :: zero_temp

!C
!C-- NODE
          integer(kind=kint)         :: n_node
          integer(kind=kint)         :: n_node_gross
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

        end type hecmwST_local_mesh


!C
!C +--------+
!C | MATRIX |
!C +--------+
!C===

        type hecmwST_index_value_pair
          integer(kind=kint) :: i
          integer(kind=kint) :: j
          real(kind=kreal), dimension(3,3) :: val
        end type hecmwST_index_value_pair

        type hecmwST_matrix_contact
          integer(kind=kint)                      :: n_val
          integer(kind=kint)                      :: max_val
          type(hecmwST_index_value_pair), pointer :: pair(:)
          logical :: checked
          logical :: sorted
          integer(kind=kint) :: max_row
          integer(kind=kint) :: max_col
        end type hecmwST_matrix_contact

        type hecmwST_matrix
          integer(kind=kint)  N, NP, NPL, NPU, NDOF

          real(kind=kreal), pointer :: D(:), B(:), X(:), ALU(:)
          real(kind=kreal), pointer :: AL(:), AU(:)
          integer(kind=kint), pointer :: indexL(:), indexU(:)
          integer(kind=kint), pointer ::  itemL(:),  itemU(:)
          integer(kind=kint ), dimension(100) :: Iarray
          real   (kind=kreal), dimension(100) :: Rarray

          type(hecmwST_matrix_contact) :: cmat
        end type hecmwST_matrix

      contains

!C
!C***
!C*** HECMW_WTIME
!C***
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
        nullify( P%master_grp_id )
        end subroutine hecmw_nullify_contact_pair

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


        call hecmw_nullify_section( P%section )
        call hecmw_nullify_material( P%material )
        call hecmw_nullify_mpc( P%mpc )
        call hecmw_nullify_amplitude( P%amp )
        call hecmw_nullify_node_grp( P%node_group )
        call hecmw_nullify_elem_grp( P%elem_group )
        call hecmw_nullify_surf_grp( P%surf_group )
        call hecmw_nullify_contact_pair( P%contact_pair )

        end subroutine hecmw_nullify_mesh


        subroutine hecmw_nullify_matrix_contact( P )
        type( hecmwST_matrix_contact ) :: P
        nullify( P%pair )
        end subroutine hecmw_nullify_matrix_contact

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

		!--- cmat
        write( nf, * ) P%cmat%n_val, P%cmat%max_val, P%cmat%max_row,  &
                       P%cmat%max_col, P%cmat%checked, P%cmat%sorted
        do i=1,P%cmat%n_val
          write( nf, * ) P%cmat%pair(i)%i,P%cmat%pair(i)%j
          write( nf, * ) P%cmat%pair(i)%val(1,1), P%cmat%pair(i)%val(1,2), P%cmat%pair(i)%val(1,3)
          write( nf, * ) P%cmat%pair(i)%val(2,1), P%cmat%pair(i)%val(2,2), P%cmat%pair(i)%val(2,3)
          write( nf, * ) P%cmat%pair(i)%val(3,1), P%cmat%pair(i)%val(3,2), P%cmat%pair(i)%val(3,3)
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
        allocate( P%X(P%NDOF*P%NP) )
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
		!---- X
		P%X = 0.d0

		!--- cmat
        read( nf, * ) P%cmat%n_val, P%cmat%max_val, P%cmat%max_row,  &
                       P%cmat%max_col, P%cmat%checked, P%cmat%sorted
        allocate( P%cmat%pair( p%cmat%max_val ) )
        do i=1,P%cmat%n_val
          read( nf, * ) P%cmat%pair(i)%i,P%cmat%pair(i)%j
          read( nf, * ) P%cmat%pair(i)%val(1,1), P%cmat%pair(i)%val(1,2), P%cmat%pair(i)%val(1,3)
          read( nf, * ) P%cmat%pair(i)%val(2,1), P%cmat%pair(i)%val(2,2), P%cmat%pair(i)%val(2,3)
          read( nf, * ) P%cmat%pair(i)%val(3,1), P%cmat%pair(i)%val(3,2), P%cmat%pair(i)%val(3,3)
        enddo

        close( nf )
        end subroutine
		
		
! ---- Following management of hecmwST_matrix_contact --	
  subroutine hecmw_cmat_init( cmat )
    type(hecmwST_matrix_contact) :: cmat

    nullify( cmat%pair )
    call hecmw_cmat_clear( cmat)
  end subroutine hecmw_cmat_init

  subroutine hecmw_cmat_finalize( cmat )
    type(hecmwST_matrix_contact) :: cmat

    if( cmat%max_val > 0 ) deallocate( cmat%pair )
    call hecmw_cmat_init( cmat )
  end subroutine hecmw_cmat_finalize

  subroutine hecmw_cmat_clear( cmat )
    type(hecmwST_matrix_contact) :: cmat

    cmat%n_val = 0
    cmat%checked = .true.
    cmat%sorted = .true.
    cmat%max_row = 0
    cmat%max_col = 0
  end subroutine hecmw_cmat_clear
  
  subroutine hecmw_cmat_copy( lhs, rhs )
    type(hecmwST_matrix_contact) :: lhs, rhs
	integer :: i

    lhs%n_val = rhs%n_val
    lhs%checked = rhs%checked
    lhs%sorted = rhs%sorted
    lhs%max_row = rhs%max_row
    lhs%max_col = rhs%max_col
    if( associated( lhs%pair )) deallocate(lhs%pair)
	if( lhs%max_val > 0 ) then
      allocate( lhs%pair(rhs%n_val) )
      do i = 1, lhs%n_val
        lhs%pair(i) = rhs%pair(i)
      enddo
    endif
  end subroutine hecmw_cmat_copy
  
! ---- End of management of hecmwST_matrix_contact --
		

! ---- Following management of hecmwST_matrix --	
!C ----------------------------------------------------------------------------
        subroutine hecMAT_init( idbg, hecMAT )
! ----------------------------------------------------------------------------
!  Purpose: Memeory allocation 
!           July.6, 2009  YUAN Xi
!  Notes: 1. Not available for mixed-dof problems
! ----------------------------------------------------------------------------
        integer :: idbg
        type( hecmwST_matrix ) :: hecMAT
        integer ::  ndof, nn, ierror
        ndof = hecMAT%NDOF
        nn = ndof*ndof
        allocate (hecMAT%AL(nn*hecMAT%NPL)        ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%AU(nn*hecMAT%NPU)        ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%B(ndof*hecMAT%NP)          ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        hecMAT%B(:)=0.d0
        allocate (hecMAT%D(nn*hecMAT%NP)          ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%X(ndof*hecMAT%NP)          ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%ALU(nn*hecMAT%N)         ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        endif
        call hecmw_cmat_init( hecMAT%cmat )
        end subroutine hecMAT_init
! ----------------------------------------------------------------------------
        subroutine hecMAT_finalize( idbg, hecMAT )
! ----------------------------------------------------------------------------
!  Purpose: Memeory allocation 
!           July.6, 2009  YUAN Xi
!  Notes: 1. Not available for mixed-dof problems
! ----------------------------------------------------------------------------
        integer :: idbg
        type( hecmwST_matrix ) :: hecMAT
        integer ::  ierror
        if( associated(hecMAT%AL) ) then
            deallocate(hecMAT%AL                  ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%AU) ) then
            deallocate(hecMAT%AU                  ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%B) ) then
            deallocate(hecMAT%B                   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
             call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%D) ) then
            deallocate(hecMAT%D                   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(HECMAT%X) ) then
            deallocate(hecMAT%X                   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%ALU) ) then 
            deallocate(hecMAT%ALU                 ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        call hecmw_cmat_finalize( hecMAT%cmat )
        end subroutine hecMAT_finalize

        !> copy hecMAT to another one		
        subroutine hecMAT_copy( idbg, lhs, rhs )
        integer :: idbg
        type( hecmwST_matrix ) :: lhs, rhs
        integer ::  ndof, nn, ierror
        ndof = rhs%NDOF
        nn = ndof*ndof
        if( .not. associated( lhs%AL ) ) then
          allocate (lhs%AL(nn*rhs%NPL)        ,STAT=ierror )
          if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
        endif
        lhs%AL = rhs%AL
        if( .not. associated( lhs%AU ) ) then
          allocate (lhs%AU(nn*rhs%NPU)        ,STAT=ierror )
          if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
        endif
        lhs%AU = rhs%AU
        if( .not. associated( lhs%B ) ) then
          allocate (lhs%B(ndof*rhs%NP)          ,STAT=ierror )
          if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
        endif
        lhs%B = rhs%B
        if( .not. associated( lhs%D ) ) then
          allocate (lhs%D(nn*rhs%NP)          ,STAT=ierror )
          if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
        endif
        lhs%D = rhs%D
        if( .not. associated( lhs%X ) ) then
          allocate (lhs%X(ndof*rhs%NP)          ,STAT=ierror )
          if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
        endif
        lhs%X = rhs%X
        if( .not. associated( lhs%ALU ) ) then
          allocate (lhs%ALU(nn*rhs%N)         ,STAT=ierror )
          if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
          endif
        endif
        lhs%ALU = rhs%ALU
        call hecmw_cmat_copy( lhs%cmat, rhs%cmat )
        end subroutine hecMAT_copy
		
        !<  clear hecMAT matrix 
        subroutine hecMAT_clear( hecMAT )
        type( hecmwST_matrix ) :: hecMAT
        hecMAT%AL= 0.d0
        hecMAT%AU= 0.d0
        hecMAT%X = 0.d0
        hecMAT%D = 0.d0
        call hecmw_cmat_clear( hecMAT%cmat )
        end subroutine hecMAT_clear

        subroutine hecmw_abort(comm)
          integer(kind=kint) :: comm

          !call MPI_ABORT(comm)
          stop
        end subroutine hecmw_abort
		
        subroutine init_boundary_grp( bgrp )
          type(fstr_boundary_grp), intent(inout) :: bgrp
          bgrp%gid = 0
          bgrp%part_id = 0
          bgrp%grp_name = ""
          bgrp%amp = ""
          bgrp%dof = 0
          bgrp%fval = 0.d0
        end subroutine
		
        subroutine copy_boundary_grp( lhs, rhs )
          type(fstr_boundary_grp), intent(in)  :: lhs
          type(fstr_boundary_grp), intent(out) :: rhs
          rhs%gid = lhs%gid
          rhs%part_id = lhs%part_id
          rhs%grp_name = lhs%grp_name
          rhs%amp = lhs%amp
          rhs%dof = lhs%dof
          rhs%fval = lhs%fval
        end subroutine
		
        subroutine init_dload_grp( bgrp )
          type(fstr_dload_grp), intent(inout) :: bgrp
          bgrp%gid = 0
          bgrp%part_id = 0
          bgrp%grp_name = ""
          bgrp%amp = ""
          bgrp%itype = -1
          bgrp%fval = 0.d0
        end subroutine
		
        subroutine copy_dload_grp( lhs, rhs )
          type(fstr_dload_grp), intent(in)  :: lhs
          type(fstr_dload_grp), intent(out) :: rhs
          rhs%gid = lhs%gid
          rhs%part_id = lhs%part_id
          rhs%grp_name = lhs%grp_name
          rhs%amp = lhs%amp
          rhs%itype = lhs%itype
          rhs%fval = lhs%fval
        end subroutine
		
        subroutine init_ndscalar_grp( bgrp )
          type(fstr_ndscalar_grp), intent(inout) :: bgrp
          bgrp%gid = 0
          bgrp%part_id = 0
          bgrp%grp_name = ""
          bgrp%amp = ""
          bgrp%fname = ""
          bgrp%fval = 0.d0
        end subroutine
		
        subroutine append_boundary_grp( oldgrp, newgrp )
          type(fstr_boundary_grp), pointer    :: oldgrp(:)
          type(fstr_boundary_grp), intent(in) :: newgrp(:)
		  
          type(fstr_boundary_grp), pointer    :: dummygrp(:)
		  
          integer :: i,j, nso, nsn
          nsn = size( newgrp )
          nso = 0
          if( associated(oldgrp) ) then
            nso=size(oldgrp)
            allocate( dummygrp(nso) )
            do i=1,nso
              call copy_boundary_grp( oldgrp(i), dummygrp(i) )
            enddo
            deallocate( oldgrp )
          endif
          allocate( oldgrp(nso+nsn) )
          do i=1,nso
              call copy_boundary_grp( dummygrp(i), oldgrp(i) )
          enddo
          do i=1,nsn
              call copy_boundary_grp( newgrp(i), oldgrp(nso+i) )
          enddo
          if( associated(dummygrp) ) deallocate(dummygrp)
        end subroutine
		
        subroutine append_dload_grp( oldgrp, newgrp )
          type(fstr_dload_grp), pointer    :: oldgrp(:)
          type(fstr_dload_grp), intent(in) :: newgrp(:)
		  
          type(fstr_dload_grp), pointer    :: dummygrp(:)
		  
          integer :: i,j, nso, nsn
          nsn = size( newgrp )
          nso = 0
          if( associated(oldgrp) ) then
            nso=size(oldgrp)
            allocate( dummygrp(nso) )
            do i=1,nso
              call copy_dload_grp( oldgrp(i), dummygrp(i) )
            enddo
            deallocate( oldgrp )
          endif
          allocate( oldgrp(nso+nsn) )
          do i=1,nso
              call copy_dload_grp( dummygrp(i), oldgrp(i) )
          enddo
          do i=1,nsn
              call copy_dload_grp( newgrp(i), oldgrp(nso+i) )
          enddo
          if( associated(dummygrp) ) deallocate(dummygrp)
        end subroutine
		
end module


