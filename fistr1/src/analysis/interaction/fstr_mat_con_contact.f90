!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions of reconstructing
!         stiffness matrix structure for the contact analysis
!         employing standard Lagrange multiplier algorithm

module fstr_matrix_con_contact

  use m_fstr
  use elementInfo
  use m_fstr_contact_damping, only: is_damping_enabled
  use m_fstr_contact_elem_alag, only: get_unique_map
  use hecmw_matrix_ass, only: hecmw_mat_profile_has_node

  implicit none
  private
  public :: fstr_get_num_lagrange_pernode
  public :: hecmwST_matrix_lagrange
  public :: fstr_save_originalMatrixStructure
  public :: fstr_mat_con_contact
  public :: fstr_s2s_profile_needs_refresh
  public :: fstr_is_matrixStruct_symmetric
  public :: fstr_is_contactALag_symmetric
  public :: fstr_is_material_symmetric
  public :: fstr_set_lagrange_diagonal
  public :: fstr_get_lagrange_diagonal

  integer(kind=kint), save         :: NPL_org, NPU_org !< original number of non-zero items
  type(nodeRelated), pointer, save :: list_nodeRelated_org(:) => null() !< original structure of matrix

  type(nodeRelated), pointer       :: list_nodeRelated(:) => null() !< current structure of matrix

  logical                          :: permission = .false.

contains

  integer(kind=kint) function fstr_get_num_lagrange_pernode(algtype)
    integer(kind=kint) :: algtype !< current loading step
    if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
      fstr_get_num_lagrange_pernode = 1
    else if( algtype == CONTACTTIED ) then
      fstr_get_num_lagrange_pernode = 3
    endif
  end function

  !> \brief This subroutine saves original matrix structure constructed originally by hecMW_matrix
  subroutine fstr_save_originalMatrixStructure(hecMAT)

    type(hecmwST_matrix) :: hecMAT !< type hecmwST_matrix

    if( associated(list_nodeRelated_org) ) return
    call hecmw_construct_nodeRelated_from_hecMAT(hecMAT, NPL_org, NPU_org, list_nodeRelated_org)

  end subroutine fstr_save_originalMatrixStructure

  !> \brief this subroutine reconstructs node-based (stiffness) matrix structure
  !> \corresponding to contact state
  subroutine fstr_mat_con_contact(cstep,contact_algo,hecMAT,fstrSOLID,hecLagMAT,infoCTChange,conMAT,is_contact_active_flag)

    integer(kind=kint)                   :: cstep !< current loading step
    integer(kind=kint)                   :: contact_algo !< current loading step
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange) :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(fstr_info_contactChange)        :: infoCTChange !< type fstr_contactChange

    integer(kind=kint)                   :: num_lagrange !< number of Lagrange multipliers
    integer(kind=kint)                   :: countNon0LU_node, countNon0LU_lagrange !< counter of node-based number of non-zero items
    integer(kind=kint)                   :: numNon0_node, numNon0_lagrange !< node-based number of displacement-related non-zero items in half of the matrix
    !< node-based number of Lagrange multiplier-related non-zero items in half of the matrix
    type (hecmwST_matrix)                :: conMAT
    logical, intent(in)                  :: is_contact_active_flag

    integer(kind=kint)                   :: i, j, grpid
    integer(kind=kint)                   :: count_n2s, count_s2s
    integer(kind=kint)                   :: nlag !< number of Lagrange multipliers per node

    count_n2s = 0
    count_s2s = 0
    do i = 1, fstrSOLID%n_contacts
      if( fstrSOLID%contacts(i)%method == CONTACTN2S ) count_n2s = count_n2s + 1
      if( fstrSOLID%contacts(i)%method == CONTACTS2S ) count_s2s = count_s2s + 1
    enddo
    count_n2s = count_n2s + fstrSOLID%n_embeds

    ! Lagrange rows are reserved for the standard-Lagrange algorithm and, additionally, for
    ! SURF-SURF mortar pairs: their sparsity reservation is keyed on one Lagrange row per
    ! active slave node (get_lag_node_list / register_pair_to_sparsity). A pure NODE-SURF
    ! augmented-Lagrange analysis keeps num_lagrange = 0, exactly as before.
    num_lagrange = 0
    if( contact_algo == kcaSLagrange .or. count_s2s > 0 ) then
      do i = 1, fstrSOLID%n_contacts
        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
        nlag = fstr_get_num_lagrange_pernode(fstrSOLID%contacts(i)%algtype)
        do j = 1, size(fstrSOLID%contacts(i)%slave)
          if( .not. is_contact_active(fstrSOLID%contacts(i)%states(j)%state) ) cycle
          num_lagrange = num_lagrange + nlag
        enddo
      enddo

      do i = 1, fstrSOLID%n_embeds
        grpid = fstrSOLID%embeds(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
        nlag = 3
        do j = 1, size(fstrSOLID%embeds(i)%slave)
          if( .not. is_contact_active(fstrSOLID%embeds(i)%states(j)%state) ) cycle
          num_lagrange = num_lagrange + nlag
        enddo
      enddo
    endif

    ! Get original list of related nodes
    call hecmw_init_nodeRelated_from_org(hecMAT%NP,num_lagrange,is_contact_active_flag,list_nodeRelated_org,list_nodeRelated)
    if( count_s2s > 0 ) call get_lag_node_list(fstrSOLID, hecLagMAT, hecMAT%NP)

    ! Construct new list of related nodes and Lagrange multipliers
    countNon0LU_node = NPL_org + NPU_org
    countNon0LU_lagrange = 0
    if( is_contact_active_flag )then
      if( count_n2s > 0 ) &
        call getNewListOFrelatednodesANDLagrangeMultipliers(cstep,contact_algo, &
        &  hecMAT%NP,fstrSOLID,countNon0LU_node,countNon0LU_lagrange,list_nodeRelated)
      if( count_s2s > 0 ) &
        call getNewListOFrelatednodesANDLagrangeMultipliers_ss(cstep,contact_algo, &
        &  hecMAT%NP,fstrSOLID,countNon0LU_node,countNon0LU_lagrange,list_nodeRelated,hecLagMAT%lag_node_table)
    endif

    ! Construct new matrix structure(hecMAT&hecLagMAT)
    numNon0_node = countNon0LU_node/2
    numNon0_lagrange = countNon0LU_lagrange/2
    call hecmw_construct_hecMAT_from_nodeRelated(hecMAT%N, hecMAT%NP, hecMAT%NDOF, &
    & numNon0_node, num_lagrange, list_nodeRelated, hecMAT)
    call hecmw_construct_hecMAT_from_nodeRelated(hecMAT%N, hecMAT%NP, hecMAT%NDOF, &
    & numNon0_node, num_lagrange, list_nodeRelated, conMAT)
    if( contact_algo == kcaSLagrange ) call hecmw_construct_hecLagMAT_from_nodeRelated(hecMAT%NP, &
    & hecMAT%NDOF, num_lagrange, numNon0_lagrange, is_contact_active_flag, list_nodeRelated, hecLagMAT)
    call hecmw_finalize_nodeRelated(list_nodeRelated)

    ! Copy Lagrange multipliers
    if( is_contact_active_flag .and. contact_algo == kcaSLagrange ) &
      call fstr_copy_lagrange_contact(fstrSOLID,hecLagMAT)

  end subroutine fstr_mat_con_contact

  subroutine get_lag_node_list(fstrSOLID, hecLagMAT, np)
    type(fstr_solid)                        :: fstrSOLID                !< type fstr_solid
    type(hecmwST_matrix_lagrange)          :: hecLagMAT            !< hecmwST_matrix_lagrange
    integer (kind=kint)                    :: np                      !< total number of nodes
    integer (kind=kint)                    :: id_lagrange, algtype, i, j, nlag, slave_node, ierr

    ! init lag_node_table
    if( associated(hecLagMAT%lag_node_table) ) deallocate(hecLagMAT%lag_node_table)
    allocate(hecLagMAT%lag_node_table(np), stat=ierr)
    if ( ierr /= 0) stop " Allocation error, hecLagMAT%lag_node_table "
    hecLagMAT%lag_node_table = 0
    id_lagrange = 0

    do i = 1, fstrSOLID%n_contacts

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      do j = 1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
        slave_node = fstrSOLID%contacts(i)%slave(j)
        hecLagMAT%lag_node_table(slave_node) = id_lagrange + 1
        id_lagrange = id_lagrange + nlag
      enddo
    enddo

    do i = 1, fstrSOLID%n_embeds
      nlag = 3
      do j = 1, size(fstrSOLID%embeds(i)%slave)
        if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle
        slave_node = fstrSOLID%embeds(i)%slave(j)
        hecLagMAT%lag_node_table(slave_node) = id_lagrange + 1
        id_lagrange = id_lagrange + nlag
      enddo
    enddo

  end subroutine get_lag_node_list

  !> Construct new list of related nodes and Lagrange multipliers. Here, a procedure similar to HEC_MW is used.
  subroutine getNewListOFrelatednodesANDLagrangeMultipliers( &
  & cstep, contact_algo, np, fstrSOLID, countNon0LU_node, countNon0LU_lagrange, list_nodeRelated )
    integer(kind=kint),intent(in)             :: cstep !< current loading step
    integer(kind=kint),intent(in)             :: contact_algo !< contact algo
    integer(kind=kint),intent(in)             :: np !< total number of nodes
    type(fstr_solid),intent(in)               :: fstrSOLID !< type fstr_solid
    integer(kind=kint), intent(inout)         :: countNon0LU_node, countNon0LU_lagrange !< counters of node-based number of non-zero items
    type(nodeRelated), pointer, intent(inout) :: list_nodeRelated(:) !< nodeRelated structure of matrix

    integer(kind=kint)            :: grpid !< contact pairs group ID
    integer(kind=kint)            :: count_lagrange !< counter of Lagrange multiplier
    integer(kind=kint)            :: ctsurf, etype, nnode, ndLocal(l_max_surface_node + 1) !< contents of type tContact
    integer(kind=kint)            :: i, j, k, nlag, algtype
    real(kind=kreal)              :: fcoeff !< friction coefficient
    logical                       :: necessary_to_insert_node, necessary_to_insert_node_pair
    logical                       :: is_contact_active_flag, is_damping_active_flag

    count_lagrange = 0
    do i = 1, fstrSOLID%n_contacts
      ! Process only NODE-SURF pairs here (mortar pairs are handled by getNewListOF..._ss).
      ! For a mortar pair states(j)%surface is -1, so master(ctsurf)%etype below would be an
      ! out-of-bounds read on a mixed deck.
      if( fstrSOLID%contacts(i)%method /= CONTACTN2S ) cycle

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      fcoeff = fstrSOLID%contacts(i)%fcoeff
      necessary_to_insert_node = ( fcoeff /= 0.0d0 .or. contact_algo == kcaALagrange )

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)
      if( contact_algo == kcaALagrange ) nlag = 1
      if( algtype == CONTACTTIED ) permission = .true.

      do j = 1, size(fstrSOLID%contacts(i)%slave)
        ! stick or sliding contact is active
        is_contact_active_flag = is_contact_active(fstrSOLID%contacts(i)%states(j)%state)
        ! damping is active
        is_damping_active_flag = fstrSOLID%contacts(i)%states(j)%state == CONTACTNEAR .and. &
          &  is_damping_enabled(fstrSOLID%contacts(i))
        
        if( is_contact_active_flag .or. is_damping_active_flag ) then

          ctsurf = fstrSOLID%contacts(i)%states(j)%surface
          etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
          if( etype/=fe_tri3n .and. etype/=fe_quad4n ) &
            stop " ##Error: This element type is not supported in contact analysis !!! "
          nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
          ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
          ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)

          ! For CONTACTNEAR damping (especially S-Lagrange + frictionless),
          ! we still need slave-master connectivity to assemble damping terms.
          necessary_to_insert_node_pair = necessary_to_insert_node .or. is_damping_active_flag

          if( is_contact_active_flag ) then
            do k=1,nlag
              if( contact_algo == kcaSLagrange ) count_lagrange = count_lagrange + 1
              call hecmw_ass_nodeRelated_from_contact_pair(np, nnode, ndLocal, count_lagrange, permission, &
              & necessary_to_insert_node_pair, list_nodeRelated_org, list_nodeRelated, countNon0LU_node, countNon0LU_lagrange )
            enddo
          else
            ! NEAR damping only: no Lagrange multiplier, insert connectivity once
            call hecmw_ass_nodeRelated_from_contact_pair(np, nnode, ndLocal, 0, permission, &
            & necessary_to_insert_node_pair, list_nodeRelated_org, list_nodeRelated, countNon0LU_node, countNon0LU_lagrange )
          endif
              
        end if

      enddo

    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      necessary_to_insert_node = ( contact_algo == kcaALagrange )

      nlag = 3
      if( contact_algo == kcaALagrange ) nlag = 1
      permission = .true.

      do j = 1, size(fstrSOLID%embeds(i)%slave)

        if( .not. is_contact_active(fstrSOLID%embeds(i)%states(j)%state) ) cycle
        ctsurf = fstrSOLID%embeds(i)%states(j)%surface
        etype = fstrSOLID%embeds(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%embeds(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%embeds(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%embeds(i)%master(ctsurf)%nodes(1:nnode)

        do k=1,nlag
          if( contact_algo == kcaSLagrange ) count_lagrange = count_lagrange + 1
          call hecmw_ass_nodeRelated_from_contact_pair(np, nnode, ndLocal, count_lagrange, permission, &
          & necessary_to_insert_node, list_nodeRelated_org, list_nodeRelated, countNon0LU_node, countNon0LU_lagrange )
        enddo
      enddo

    enddo

  end subroutine getNewListOFrelatednodesANDLagrangeMultipliers

  !> Copy Lagrange multipliers
  subroutine fstr_copy_lagrange_contact(fstrSOLID,hecLagMAT)

    type(fstr_solid)                        :: fstrSOLID                !< type fstr_solid
    type(hecmwST_matrix_lagrange)          :: hecLagMAT            !< hecmwST_matrix_lagrange
    integer (kind=kint)                    :: id_lagrange, algtype, i, j, k, nlag, slave_node

    id_lagrange = 0

    do i = 1, fstrSOLID%n_contacts

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      do j = 1, size(fstrSOLID%contacts(i)%slave)
        if( .not. is_contact_active(fstrSOLID%contacts(i)%states(j)%state) ) cycle
        slave_node = fstrSOLID%contacts(i)%slave(j)
        hecLagMAT%lag_node_table(slave_node) = id_lagrange + 1
        do k=1,nlag
          id_lagrange = id_lagrange + 1
          hecLagMAT%Lagrange(id_lagrange)=fstrSOLID%contacts(i)%states(j)%multiplier(k)
        enddo
      enddo
    enddo

    do i = 1, fstrSOLID%n_embeds
      nlag = 3
      do j = 1, size(fstrSOLID%embeds(i)%slave)
        if( .not. is_contact_active(fstrSOLID%embeds(i)%states(j)%state) ) cycle
        slave_node = fstrSOLID%embeds(i)%slave(j)
        hecLagMAT%lag_node_table(slave_node) = id_lagrange + 1
        do k=1,nlag
          id_lagrange = id_lagrange + 1
          hecLagMAT%Lagrange(id_lagrange)=fstrSOLID%embeds(i)%states(j)%multiplier(k)
        enddo
      enddo
    enddo

  end subroutine fstr_copy_lagrange_contact

  !> \brief this function judges whether sitiffness matrix is symmetric or not
  logical function fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)

    type(fstr_solid )        :: fstrSOLID
    type(hecmwST_local_mesh) :: hecMESH
    integer (kind=kint)      :: is_in_contact

    is_in_contact = 0
    if( fstrSOLID%n_contacts>0 ) then
      if( any(fstrSOLID%contacts(:)%fcoeff /= 0.0d0) )  is_in_contact = 1
    endif
    call hecmw_allreduce_I1(hecMESH, is_in_contact, HECMW_MAX)
    if( is_in_contact == 0 .and. hecMESH%n_dof /= 4 .and. fstr_is_material_symmetric(fstrSOLID,hecMESH) ) then
      fstr_is_matrixStruct_symmetric = .true.
    else
      fstr_is_matrixStruct_symmetric = .false.
    endif

  end function fstr_is_matrixStruct_symmetric

  !> \brief this function judges whether the ALagrange contact tangent is symmetric or not
  logical function fstr_is_contactALag_symmetric(fstrSOLID,hecMESH)

    type(fstr_solid )        :: fstrSOLID
    type(hecmwST_local_mesh) :: hecMESH
    integer (kind=kint)      :: is_unsymmetric

    ! the ALagrange contact terms are symmetric while the friction cone radius stays frozen
    ! at the multiplier; !CONTACT_ALGO, FRICTION_CONE=FOLLOW clears contact%symmetric and
    ! the tangent then carries the coupling block of getContactStiffness_Alag
    is_unsymmetric = 0
    if( fstrSOLID%n_contacts>0 ) then
      if( any( fstrSOLID%contacts(:)%fcoeff /= 0.0d0 .and. .not.fstrSOLID%contacts(:)%symmetric ) ) is_unsymmetric = 1
    endif
    call hecmw_allreduce_I1(hecMESH, is_unsymmetric, HECMW_MAX)
    fstr_is_contactALag_symmetric = ( is_unsymmetric == 0 ) .and. fstr_is_material_symmetric(fstrSOLID,hecMESH)

  end function fstr_is_contactALag_symmetric

  !> \brief this function judges whether all materials yield a symmetric tangent stiffness
  logical function fstr_is_material_symmetric(fstrSOLID,hecMESH)

    type(fstr_solid )        :: fstrSOLID
    type(hecmwST_local_mesh) :: hecMESH
    integer (kind=kint)      :: is_unsymmetric, i, ytype

    ! non-associated flow (dilatancy angle psi /= friction angle phi) makes the consistent
    ! tangent unsymmetric; when psi is omitted the parser stores psi=phi, so the exact
    ! comparison below classifies that case as associated (symmetric)
    is_unsymmetric = 0
    if( associated(fstrSOLID%materials) ) then
      do i = 1, size(fstrSOLID%materials)
        ytype = getYieldFunction( fstrSOLID%materials(i)%mtype )
        if( ytype == 1 ) then      ! Mohr-Coulomb: PLCONST3=phi, PLCONST4=psi
          if( fstrSOLID%materials(i)%variables(M_PLCONST3) /= fstrSOLID%materials(i)%variables(M_PLCONST4) ) &
            is_unsymmetric = 1
        elseif( ytype == 2 ) then  ! Drucker-Prager: PLCONST3=eta(phi), PLCONST5=etabar(psi)
          if( fstrSOLID%materials(i)%variables(M_PLCONST3) /= fstrSOLID%materials(i)%variables(M_PLCONST5) ) &
            is_unsymmetric = 1
        endif
      enddo
    endif
    call hecmw_allreduce_I1(hecMESH, is_unsymmetric, HECMW_MAX)
    fstr_is_material_symmetric = (is_unsymmetric == 0)

  end function fstr_is_material_symmetric

  !> \brief Set diagonal component value for specified Lagrange multiplier
  subroutine fstr_set_lagrange_diagonal(hecLagMAT, ilag, value)
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< hecmwST_matrix_lagrange
    integer(kind=kint), intent(in) :: ilag !< Lagrange multiplier index (1-based)
    real(kind=kreal), intent(in) :: value !< diagonal component value to set

    if (ilag < 1 .or. ilag > hecLagMAT%num_lagrange) then
      write(*,*) 'Error in fstr_set_lagrange_diagonal: invalid Lagrange multiplier index', ilag
      stop
    endif

    if (.not. associated(hecLagMAT%D_lagrange)) then
      write(*,*) 'Error in fstr_set_lagrange_diagonal: D_lagrange not allocated'
      stop
    endif

    hecLagMAT%D_lagrange(ilag) = value

  end subroutine fstr_set_lagrange_diagonal

  !> \brief Get diagonal component value for specified Lagrange multiplier
  real(kind=kreal) function fstr_get_lagrange_diagonal(hecLagMAT, ilag)
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< hecmwST_matrix_lagrange
    integer(kind=kint), intent(in) :: ilag !< Lagrange multiplier index (1-based)

    if (ilag < 1 .or. ilag > hecLagMAT%num_lagrange) then
      write(*,*) 'Error in fstr_get_lagrange_diagonal: invalid Lagrange multiplier index', ilag
      fstr_get_lagrange_diagonal = 0.0d0
      return
    endif

    if (.not. associated(hecLagMAT%D_lagrange)) then
      write(*,*) 'Error in fstr_get_lagrange_diagonal: D_lagrange not allocated'
      fstr_get_lagrange_diagonal = 0.0d0
      return
    endif

    fstr_get_lagrange_diagonal = hecLagMAT%D_lagrange(ilag)

  end function fstr_get_lagrange_diagonal


  subroutine getNewListOFrelatednodesANDLagrangeMultipliers_ss( &
      & cstep, contact_algo, np, fstrSOLID, countNon0LU_node, countNon0LU_lagrange, list_nodeRelated, lag_node_table )
    integer(kind=kint),intent(in)             :: cstep !< current loading step
    integer(kind=kint),intent(in)             :: contact_algo !< contact algo
    integer(kind=kint),intent(in)             :: np !< total number of nodes
    type(fstr_solid),intent(in)               :: fstrSOLID !< type fstr_solid
    integer(kind=kint), intent(inout)         :: countNon0LU_node, countNon0LU_lagrange !< counters of node-based non-zero items
    type(nodeRelated), pointer, intent(inout) :: list_nodeRelated(:) !< nodeRelated structure of matrix
    integer(kind=kint), intent(in)         :: lag_node_table(:) !< table of Lagrange multipliers
    integer(kind=kint)            :: grpid !< contact pairs group ID
    integer(kind=kint)            :: ctsurf, nsurf !< contents of type tContact
    integer(kind=kint)            :: i, j, m
    integer(kind=kint)            :: g, unique_count !< unique (slave_surf, master) pair iteration
    integer(kind=kint), allocatable :: maplist(:), master_idxs(:) !< unique master mapping from get_unique_map
    real(kind=kreal)              :: fcoeff !< friction coefficient
    logical                       :: necessary_to_insert_node
    permission = .true.

    do i = 1, fstrSOLID%n_contacts
      if( fstrSOLID%contacts(i)%method /= CONTACTS2S ) cycle
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      fcoeff = fstrSOLID%contacts(i)%fcoeff
      necessary_to_insert_node = ( fcoeff /= 0.0d0 .or. contact_algo == kcaALagrange )

      do j = 1, size(fstrSOLID%contacts(i)%slave_surf)
        if( fstrSOLID%contacts(i)%slave_surf(j)%state == CONTACTFREE ) cycle

        if( contact_algo == kcaALagrange ) then
          ! ALagrange SS: register sparsity once per unique (slave_surf, master) pair.
          ! sparsity_expansion == SPARSITY_NONE:     current master element only (matrix rebuilt on contact2neighbor)
          ! sparsity_expansion == SPARSITY_NEIGHBOR: current + neighbor master elements (matrix rebuilt only on contact2beyond)
          call get_unique_map( fstrSOLID%contacts(i)%slave_surf(j), maplist, master_idxs, unique_count )
          do g = 1, unique_count
            ctsurf = master_idxs(g)
            ! Register current master element
            call register_pair_to_sparsity( np, fstrSOLID%contacts(i)%slave_surf(j)%nodes, &
              fstrSOLID%contacts(i)%master(ctsurf), lag_node_table, necessary_to_insert_node, &
              countNon0LU_node, countNon0LU_lagrange, list_nodeRelated )

            ! SPARSITY_NEIGHBOR: also register neighbor master elements
            if( fstrSOLID%contacts(i)%sparsity_expansion == SPARSITY_NEIGHBOR ) then
              do m = 1, fstrSOLID%contacts(i)%master(ctsurf)%n_neighbor
                nsurf = fstrSOLID%contacts(i)%master(ctsurf)%neighbor(m)
                call register_pair_to_sparsity( np, fstrSOLID%contacts(i)%slave_surf(j)%nodes, &
                  fstrSOLID%contacts(i)%master(nsurf), lag_node_table, necessary_to_insert_node, &
                  countNon0LU_node, countNon0LU_lagrange, list_nodeRelated )
              enddo
            endif
          enddo
          deallocate( maplist, master_idxs )
        endif
      enddo
    enddo

  end subroutine getNewListOFrelatednodesANDLagrangeMultipliers_ss

  !> \brief Register sparsity coupling of one (slave_surf, master) pair.
  !> For each slave node, register coupling with the other slave nodes AND the master
  !> nodes. Slave-slave cross terms are required for consistent mortar tangent stiffness.
  subroutine register_pair_to_sparsity( np, slave_nodes, master_surf, lag_node_table, &
      & necessary_to_insert_node, countNon0LU_node, countNon0LU_lagrange, list_nodeRelated )
    integer(kind=kint), intent(in)            :: np !< total number of nodes
    integer(kind=kint), intent(in)            :: slave_nodes(:)  !< slave surface node ids
    type(tSurfElement), intent(in)            :: master_surf !< master surface id
    integer(kind=kint), intent(in)            :: lag_node_table(:) !< table of Lagrange multipliers
    logical, intent(in)                       :: necessary_to_insert_node
    integer(kind=kint), intent(inout)         :: countNon0LU_node, countNon0LU_lagrange
    type(nodeRelated), pointer, intent(inout) :: list_nodeRelated(:)

    integer(kind=kint) :: nnode_s, nnode_m, l, m, idx, id_lag, nnode_pair, etype
    integer(kind=kint) :: ndLocal(2*l_max_surface_node + 1)

    etype = master_surf%etype
    if( etype/=fe_tri3n .and. etype/=fe_quad4n ) stop " ##Error: This element type is not supported in contact analysis !!! "
    nnode_s = size(slave_nodes)
    nnode_m = size(master_surf%nodes)
    do l = 1, nnode_s
      ndLocal(1) = slave_nodes(l)
      id_lag = lag_node_table( ndLocal(1) ) - 1
      ! ndLocal: slave_l + other_slave_nodes + master_nodes
      idx = 1
      do m = 1, nnode_s
        if( m == l ) cycle
        idx = idx + 1
        ndLocal(idx) = slave_nodes(m)
      enddo
      ndLocal(idx+1:idx+nnode_m) = master_surf%nodes(:)
      nnode_pair = nnode_s - 1 + nnode_m
      call hecmw_ass_nodeRelated_from_contact_pair( np, nnode_pair, ndLocal, id_lag, permission, &
        & necessary_to_insert_node, list_nodeRelated_org, list_nodeRelated, countNon0LU_node, countNon0LU_lagrange )
    enddo
  end subroutine register_pair_to_sparsity

  !> \brief S2S mortar profile-invariant check (read-only companion to the reservation in
  !> getNewListOFrelatednodesANDLagrangeMultipliers_ss / register_pair_to_sparsity).
  !>
  !> Returns .true. if any currently-active (slave_surf, master) mortar pair has a
  !> slave-node x master-node coupling that is NOT present in the matrix profile of conMAT.
  !> Such a missing coupling means the active master has drifted beyond the reserved 1-ring
  !> via a silent facet-hop / re-association (one that does not increment the structural
  !> change counters, so fstr_is_matrixStructure_changed stays false), and the next contact
  !> stiffness assembly (calcu_contact_stiffness_SurfSurf -> hecmw_mat_ass_elem) would hit
  !> an out-of-profile connectivity and abort in hecmw_mat_add_node. Reporting it here lets
  !> the caller force a matrix rebuild first, restoring the invariant.
  !>
  !> Only slave x master CROSS couplings are checked: slave-slave and master-master
  !> couplings are base FE element couplings (a contact surface element is a face of one
  !> solid element, so its nodes are mutually coupled in list_nodeRelated_org and can never
  !> drift), whereas the cross terms are exactly what register_pair_to_sparsity adds and the
  !> only ones that can be missing. This is necessary and sufficient and avoids any spurious
  !> every-iteration refresh. The active master set is obtained from get_unique_map, the same
  !> enumeration getIntGap uses inside the assembly, so the checked footprint matches exactly.
  !>
  !> method-gated to CONTACTS2S + ALag (N2S and SLagrange paths untouched). Per-rank (local)
  !> decision, mirroring the existing per-rank fstr_is_matrixStructure_changed gate; the
  !> collective solver re-init rides the existing contact_changed_global allreduce.
  logical function fstr_s2s_profile_needs_refresh( cstep, contact_algo, fstrSOLID, conMAT )
    integer(kind=kint), intent(in)   :: cstep         !< current loading step
    integer(kind=kint), intent(in)   :: contact_algo  !< contact algorithm (kcaALagrange/kcaSLagrange)
    type(fstr_solid), intent(in)     :: fstrSOLID     !< type fstr_solid
    type(hecmwST_matrix), intent(in) :: conMAT        !< contact matrix (S2S stiffness assembly target)

    integer(kind=kint) :: i, j, g, ctsurf, grpid, unique_count
    integer(kind=kint), allocatable :: maplist(:), master_idxs(:)

    fstr_s2s_profile_needs_refresh = .false.
    if( contact_algo /= kcaALagrange ) return

    do i = 1, fstrSOLID%n_contacts
      if( fstrSOLID%contacts(i)%method /= CONTACTS2S ) cycle
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      do j = 1, size(fstrSOLID%contacts(i)%slave_surf)
        if( fstrSOLID%contacts(i)%slave_surf(j)%state == CONTACTFREE ) cycle

        call get_unique_map( fstrSOLID%contacts(i)%slave_surf(j), maplist, master_idxs, unique_count )
        do g = 1, unique_count
          ctsurf = master_idxs(g)
          if( .not. pair_cross_profile_complete( fstrSOLID%contacts(i)%slave_surf(j)%nodes, &
              &  fstrSOLID%contacts(i)%master(ctsurf), conMAT ) ) then
            fstr_s2s_profile_needs_refresh = .true.
            deallocate( maplist, master_idxs )
            return
          endif
        enddo
        deallocate( maplist, master_idxs )
      enddo
    enddo
  end function fstr_s2s_profile_needs_refresh

  !> \brief True if every slave-node x master-node coupling of one (slave_surf, master)
  !> mortar pair is present in the conMAT profile. See fstr_s2s_profile_needs_refresh for
  !> why only cross couplings are checked (slave-slave / master-master are base FE terms).
  logical function pair_cross_profile_complete( slave_nodes, master_surf, conMAT )
    integer(kind=kint), intent(in)   :: slave_nodes(:)  !< slave surface node ids
    type(tSurfElement), intent(in)   :: master_surf     !< master surface element
    type(hecmwST_matrix), intent(in) :: conMAT          !< contact matrix

    integer(kind=kint) :: s, m

    pair_cross_profile_complete = .true.
    do s = 1, size(slave_nodes)
      do m = 1, size(master_surf%nodes)
        if( .not. hecmw_mat_profile_has_node( conMAT, slave_nodes(s), master_surf%nodes(m) ) ) then
          pair_cross_profile_complete = .false.
          return
        endif
      enddo
    enddo
  end function pair_cross_profile_complete

end module fstr_matrix_con_contact
