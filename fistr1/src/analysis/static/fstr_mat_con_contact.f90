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
  public :: fstrST_matrix_contact_lagrange
  public :: fstr_save_originalMatrixStructure
  public :: fstr_mat_con_contact
  public :: fstr_is_matrixStruct_symmetric

  !> Structure for Lagrange multiplier-related part of stiffness matrix
  !> (Lagrange multiplier-related matrix)
  type fstrST_matrix_contact_lagrange
    integer(kind=kint) :: num_lagrange = 0 !< total number of Lagrange multipliers
    integer(kind=kint) :: numL_lagrange = 0, numU_lagrange = 0 !< node-based number of non-zero items in lower triangular half of matrix
    !< node-based number of non-zero items in upper triangular half of matrix
    integer(kind=kint), pointer  :: indexL_lagrange(:) => null(), &
      indexU_lagrange(:) => null() !< node-based index of first non-zero item of each row
    integer(kind=kint), pointer  :: itemL_lagrange(:) => null(), &
      itemU_lagrange(:) => null() !< node-based column number of non-zero items
    real(kind=kreal),    pointer  :: AL_lagrange(:) => null(), &
      AU_lagrange(:) => null() !< values of non-zero items
    real(kind=kreal),    pointer  :: Lagrange(:) => null() !< values of Lagrange multipliers
  end type fstrST_matrix_contact_lagrange

  integer(kind=kint), save         :: NPL_org, NPU_org !< original number of non-zero items
  type(nodeRelated), pointer, save :: list_nodeRelated_org(:) => null() !< original structure of matrix

  type(nodeRelated), pointer       :: list_nodeRelated(:) => null() !< current structure of matrix

  logical                          :: permission = .false.

contains


  !> \brief This subroutine saves original matrix structure constructed originally by hecMW_matrix
  subroutine fstr_save_originalMatrixStructure(hecMAT)

    type(hecmwST_matrix) :: hecMAT !< type hecmwST_matrix

    if( associated(list_nodeRelated_org) ) return
    call hecmw_construct_nodeRelated_from_hecMAT(hecMAT, NPL_org, NPU_org, list_nodeRelated_org)

  end subroutine fstr_save_originalMatrixStructure

  !> \brief this subroutine reconstructs node-based (stiffness) matrix structure
  !> \corresponding to contact state
  subroutine fstr_mat_con_contact(cstep,hecMAT,fstrSOLID,fstrMAT,infoCTChange,conMAT,is_contact_active)

    integer(kind=kint)                   :: cstep !< current loading step
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(fstrST_matrix_contact_lagrange) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    type(fstr_info_contactChange)        :: infoCTChange !< type fstr_contactChange

    integer(kind=kint)                   :: num_lagrange !< number of Lagrange multipliers
    integer(kind=kint)                   :: countNon0LU_node, countNon0LU_lagrange !< counter of node-based number of non-zero items
    integer(kind=kint)                   :: numNon0_node, numNon0_lagrange !< node-based number of displacement-related non-zero items in half of the matrix
    !< node-based number of Lagrange multiplier-related non-zero items in half of the matrix
    type (hecmwST_matrix)                :: conMAT
    logical, intent(in)                  :: is_contact_active

    num_lagrange = infoCTChange%contactNode_current

    ! Get original list of related nodes
    call hecmw_init_nodeRelated_from_org(hecMAT%NP,num_lagrange,is_contact_active,list_nodeRelated_org,list_nodeRelated)

    ! Construct new list of related nodes and Lagrange multipliers
    countNon0LU_node = NPL_org + NPU_org
    countNon0LU_lagrange = 0
    if( is_contact_active ) call getNewListOFrelatednodesANDLagrangeMultipliers(cstep,hecMAT%NP,fstrSOLID, &
       &  countNon0LU_node,countNon0LU_lagrange,list_nodeRelated)

    ! Construct new matrix structure(hecMAT&fstrMAT)
    numNon0_node = countNon0LU_node/2
    numNon0_lagrange = countNon0LU_lagrange/2
    call hecmw_construct_hecMAT_from_nodeRelated(hecMAT%N, hecMAT%NP, hecMAT%NDOF, &
      & numNon0_node, num_lagrange, list_nodeRelated, hecMAT)
    call hecmw_construct_hecMAT_from_nodeRelated(hecMAT%N, hecMAT%NP, hecMAT%NDOF, &
      & numNon0_node, num_lagrange, list_nodeRelated, conMAT)
    call construct_fstrMAT_from_nodeRelated(hecMAT%NP, hecMAT%NDOF, num_lagrange, numNon0_lagrange, &
      & is_contact_active, list_nodeRelated, fstrMAT)
    call hecmw_finalize_nodeRelated(list_nodeRelated)

    ! Copy Lagrange multipliers
    if( is_contact_active ) &
      call fstr_copy_lagrange_contact(fstrSOLID,fstrMAT)

  end subroutine fstr_mat_con_contact

  !> Construct new list of related nodes and Lagrange multipliers. Here, a procedure similar to HEC_MW is used.
  subroutine getNewListOFrelatednodesANDLagrangeMultipliers( &
      & cstep, np, fstrSOLID, countNon0LU_node, countNon0LU_lagrange, list_nodeRelated )
    integer(kind=kint),intent(in)             :: cstep !< current loading step
    integer(kind=kint),intent(in)             :: np !< total number of nodes
    type(fstr_solid),intent(in)               :: fstrSOLID !< type fstr_solid
    integer(kind=kint), intent(inout)         :: countNon0LU_node, countNon0LU_lagrange !< counters of node-based number of non-zero items
    type(nodeRelated), pointer, intent(inout) :: list_nodeRelated(:) !< nodeRelated structure of matrix

    integer(kind=kint)            :: grpid !< contact pairs group ID
    integer(kind=kint)            :: count_lagrange !< counter of Lagrange multiplier
    integer(kind=kint)            :: ctsurf, etype, nnode, ndLocal(l_max_surface_node + 1) !< contants of type tContact
    integer(kind=kint)            :: i, j
    real(kind=kreal)              :: fcoeff !< friction coefficient
    logical                       :: necessary_to_insert_node

    count_lagrange = 0
    do i = 1, size(fstrSOLID%contacts)

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      fcoeff = fstrSOLID%contacts(i)%fcoeff
      necessary_to_insert_node = ( fcoeff /= 0.0d0 )

      do j = 1, size(fstrSOLID%contacts(i)%slave)

        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
        ctsurf = fstrSOLID%contacts(i)%states(j)%surface
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
        if( etype/=fe_tri3n .and. etype/=fe_quad4n ) &
          stop " ##Error: This element type is not supported in contact analysis !!! "
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)

        count_lagrange = count_lagrange + 1

        call hecmw_ass_nodeRelated_from_contact_pair(np, nnode, ndLocal, count_lagrange, permission, &
          & necessary_to_insert_node, list_nodeRelated_org, list_nodeRelated, countNon0LU_node, countNon0LU_lagrange )
      enddo

    enddo

  end subroutine getNewListOFrelatednodesANDLagrangeMultipliers

  !> Construct fstrMAT structure
  subroutine construct_fstrMAT_from_nodeRelated(np,ndof,num_lagrange,numNon0_lagrange,is_contact_active,list_nodeRelated,fstrMAT)
    integer(kind=kint), intent(in)                      :: np
    integer(kind=kint), intent(in)                      :: ndof
    integer(kind=kint), intent(in)                      :: num_lagrange
    integer(kind=kint), intent(in)                      :: numNon0_lagrange
    logical, intent(in)                                 :: is_contact_active
    type(nodeRelated), pointer, intent(in)              :: list_nodeRelated(:) !< nodeRelated structure of matrix
    type(fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT             !< type fstrST_matrix_contact_lagrange

    integer(kind=kint)  :: countNon0U_lagrange, countNon0L_lagrange !< counters of node-based number of non-zero items
    integer(kind=kint)  :: numI_lagrange
    integer(kind=kint)  :: i, j, ierr

    fstrMAT%num_lagrange = num_lagrange
    fstrMAT%numL_lagrange = numNon0_lagrange
    fstrMAT%numU_lagrange = numNon0_lagrange

    if(associated(fstrMAT%indexL_lagrange)) deallocate(fstrMAT%indexL_lagrange)
    if(associated(fstrMAT%indexU_lagrange)) deallocate(fstrMAT%indexU_lagrange)
    if(associated(fstrMAT%itemL_lagrange)) deallocate(fstrMAT%itemL_lagrange)
    if(associated(fstrMAT%itemU_lagrange)) deallocate(fstrMAT%itemU_lagrange)
    if(associated(fstrMAT%AL_lagrange)) deallocate(fstrMAT%AL_lagrange)
    if(associated(fstrMAT%AU_lagrange)) deallocate(fstrMAT%AU_lagrange)
    if(associated(fstrMAT%Lagrange)) deallocate(fstrMAT%Lagrange)

    if( is_contact_active ) then
      ! init indexU_lagrange
      allocate(fstrMAT%indexU_lagrange(0:np), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, fstrMAT%indexU_lagrange "
      fstrMAT%indexU_lagrange = 0

      ! init itemU_lagrange
      allocate(fstrMAT%itemU_lagrange(numNon0_lagrange), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, fstrMAT%itemU_lagrange "
      fstrMAT%itemU_lagrange = 0

      ! setup upper lagrange CRS matrix
      countNon0U_lagrange = 0
      do i = 1, np
        numI_lagrange = list_nodeRelated(i)%num_lagrange
        do j = 1, numI_lagrange
          countNon0U_lagrange = countNon0U_lagrange + 1
          fstrMAT%itemU_lagrange(countNon0U_lagrange) = list_nodeRelated(i)%id_lagrange(j)
        enddo
        fstrMAT%indexU_lagrange(i) = countNon0U_lagrange
      end do

      ! init indexL_lagrange
      allocate(fstrMAT%indexL_lagrange(0:num_lagrange), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, fstrMAT%indexL_lagrange "
      fstrMAT%indexL_lagrange = 0

      ! init itemL_lagrange
      allocate(fstrMAT%itemL_lagrange(numNon0_lagrange), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, fstrMAT%itemL_lagrange "
      fstrMAT%itemL_lagrange = 0

      ! setup lower lagrange CRS matrix
      countNon0L_lagrange = 0
      do i = 1, num_lagrange
        numI_lagrange = list_nodeRelated(np+i)%num_lagrange
        do j = 1, numI_lagrange
          countNon0L_lagrange = countNon0L_lagrange + 1
          fstrMAT%itemL_lagrange(countNon0L_lagrange) = list_nodeRelated(np+i)%id_lagrange(j)
        enddo
        fstrMAT%indexL_lagrange(i) = countNon0L_lagrange
      enddo

      ! init AU_lagrange
      allocate(fstrMAT%AU_lagrange(ndof*numNon0_lagrange), stat=ierr)
      if ( ierr /= 0 ) stop " Allocation error, fstrMAT%AU_lagrange "
      fstrMAT%AU_lagrange = 0.0D0

      ! init AL_lagrange
      allocate(fstrMAT%AL_lagrange(ndof*numNon0_lagrange), stat=ierr)
      if ( ierr /= 0 ) stop " Allocation error, fstrMAT%AL_lagrange "
      fstrMAT%AL_lagrange = 0.0D0

      ! init Lagrange
      allocate(fstrMAT%Lagrange(num_lagrange))
      fstrMAT%Lagrange = 0.0D0
    endif

  end subroutine

  !> Copy Lagrange multipliers
  subroutine fstr_copy_lagrange_contact(fstrSOLID,fstrMAT)

    type(fstr_solid)                        :: fstrSOLID                !< type fstr_solid
    type(fstrST_matrix_contact_lagrange)    :: fstrMAT                  !< fstrST_matrix_contact_lagrange
    integer (kind=kint)                    :: id_lagrange, i, j

    id_lagrange = 0

    do i = 1, size(fstrSOLID%contacts)
      do j = 1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
        id_lagrange = id_lagrange + 1
        fstrMAT%Lagrange(id_lagrange)=fstrSOLID%contacts(i)%states(j)%multiplier(1)
      enddo
    enddo

  end subroutine fstr_copy_lagrange_contact

  !> \brief this function judges whether sitiffness matrix is symmetric or not
  logical function fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)

    type(fstr_solid )        :: fstrSOLID
    type(hecmwST_local_mesh) :: hecMESH
    integer (kind=kint)      :: is_in_contact

    is_in_contact = 0
    if( any(fstrSOLID%contacts(:)%fcoeff /= 0.0d0) )  &
      is_in_contact = 1
    call hecmw_allreduce_I1(hecMESH, is_in_contact, HECMW_MAX)
    fstr_is_matrixStruct_symmetric = (is_in_contact == 0)

  end function fstr_is_matrixStruct_symmetric


end module fstr_matrix_con_contact
