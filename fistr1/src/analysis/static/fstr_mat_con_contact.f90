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

  implicit none
  private
  public :: fstr_get_num_lagrange_pernode
  public :: hecmwST_matrix_lagrange
  public :: fstr_save_originalMatrixStructure
  public :: fstr_mat_con_contact
  public :: fstr_is_matrixStruct_symmetric

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
  subroutine fstr_mat_con_contact(cstep,contact_algo,hecMAT,fstrSOLID,hecLagMAT,infoCTChange,conMAT,is_contact_active)

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
    logical, intent(in)                  :: is_contact_active

    integer(kind=kint)                   :: i, j, grpid
    integer(kind=kint)                   :: nlag !< number of Lagrange multipliers per node

    num_lagrange = 0
    if( contact_algo == kcaSLagrange ) then
      do i = 1, fstrSOLID%n_contacts
        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
        nlag = fstr_get_num_lagrange_pernode(fstrSOLID%contacts(i)%algtype)
        do j = 1, size(fstrSOLID%contacts(i)%slave)
          if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
          num_lagrange = num_lagrange + nlag
        enddo
      enddo

      do i = 1, fstrSOLID%n_embeds
        grpid = fstrSOLID%embeds(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
        nlag = 3
        do j = 1, size(fstrSOLID%embeds(i)%slave)
          if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle
          num_lagrange = num_lagrange + nlag
        enddo
      enddo
  endif

    ! Get original list of related nodes
    call hecmw_init_nodeRelated_from_org(hecMAT%NP,num_lagrange,is_contact_active,list_nodeRelated_org,list_nodeRelated)

    ! Construct new list of related nodes and Lagrange multipliers
    countNon0LU_node = NPL_org + NPU_org
    countNon0LU_lagrange = 0
    if( is_contact_active ) call getNewListOFrelatednodesANDLagrangeMultipliers(cstep,contact_algo, &
       &  hecMAT%NP,fstrSOLID,countNon0LU_node,countNon0LU_lagrange,list_nodeRelated)

    ! Construct new matrix structure(hecMAT&hecLagMAT)
    numNon0_node = countNon0LU_node/2
    numNon0_lagrange = countNon0LU_lagrange/2
    call hecmw_construct_hecMAT_from_nodeRelated(hecMAT%N, hecMAT%NP, hecMAT%NDOF, &
      & numNon0_node, num_lagrange, list_nodeRelated, hecMAT)
    call hecmw_construct_hecMAT_from_nodeRelated(hecMAT%N, hecMAT%NP, hecMAT%NDOF, &
      & numNon0_node, num_lagrange, list_nodeRelated, conMAT)
    if( contact_algo == kcaSLagrange ) call hecmw_construct_hecLagMAT_from_nodeRelated(hecMAT%NP, &
      & hecMAT%NDOF, num_lagrange, numNon0_lagrange, is_contact_active, list_nodeRelated, hecLagMAT)
    call hecmw_finalize_nodeRelated(list_nodeRelated)

    ! Copy Lagrange multipliers
    if( is_contact_active .and. contact_algo == kcaSLagrange ) &
      call fstr_copy_lagrange_contact(fstrSOLID,hecLagMAT)

  end subroutine fstr_mat_con_contact

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
    logical                       :: necessary_to_insert_node

    count_lagrange = 0
    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      fcoeff = fstrSOLID%contacts(i)%fcoeff
      necessary_to_insert_node = ( fcoeff /= 0.0d0 .or. contact_algo == kcaALagrange )

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)
      if( contact_algo == kcaALagrange ) nlag = 1
      if( algtype == CONTACTTIED ) permission = .true.

      do j = 1, size(fstrSOLID%contacts(i)%slave)

        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
        ctsurf = fstrSOLID%contacts(i)%states(j)%surface
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
        if( etype/=fe_tri3n .and. etype/=fe_quad4n ) &
          stop " ##Error: This element type is not supported in contact analysis !!! "
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)

        do k=1,nlag
          if( contact_algo == kcaSLagrange ) count_lagrange = count_lagrange + 1
          call hecmw_ass_nodeRelated_from_contact_pair(np, nnode, ndLocal, count_lagrange, permission, &
            & necessary_to_insert_node, list_nodeRelated_org, list_nodeRelated, countNon0LU_node, countNon0LU_lagrange )
        enddo
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

        if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle
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
    integer (kind=kint)                    :: id_lagrange, algtype, i, j, k, nlag

    id_lagrange = 0

    do i = 1, fstrSOLID%n_contacts

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      do j = 1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
        do k=1,nlag
          id_lagrange = id_lagrange + 1
          hecLagMAT%Lagrange(id_lagrange)=fstrSOLID%contacts(i)%states(j)%multiplier(k)
        enddo
      enddo
    enddo

    do i = 1, fstrSOLID%n_embeds
      nlag = 3
      do j = 1, size(fstrSOLID%embeds(i)%slave)
        if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle
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
    fstr_is_matrixStruct_symmetric = (is_in_contact == 0)

  end function fstr_is_matrixStruct_symmetric


end module fstr_matrix_con_contact
