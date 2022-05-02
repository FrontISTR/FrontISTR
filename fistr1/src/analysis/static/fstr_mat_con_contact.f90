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
    if( is_contact_active ) &
      call getNewListOFrelatednodesANDLagrangeMultipliers(cstep,hecMAT%NP,fstrSOLID,countNon0LU_node,countNon0LU_lagrange)

    ! Construct new matrix structure(hecMAT&fstrMAT)
    numNon0_node = countNon0LU_node/2
    numNon0_lagrange = countNon0LU_lagrange/2
    call constructNewMatrixStructure(hecMAT,fstrMAT,num_lagrange,numNon0_node,numNon0_lagrange,conMAT,is_contact_active)

    ! Copy Lagrange multipliers
    if( is_contact_active ) &
      call fstr_copy_lagrange_contact(fstrSOLID,fstrMAT)

  end subroutine fstr_mat_con_contact

  !> Construct new list of related nodes and Lagrange multipliers. Here, a procedure similar to HEC_MW is used.
  subroutine getNewListOFrelatednodesANDLagrangeMultipliers(cstep,np,fstrSOLID,countNon0LU_node,countNon0LU_lagrange)
    integer(kind=kint),intent(in)     :: cstep !< current loading step
    integer(kind=kint),intent(in)     :: np !< total number of nodes
    type(fstr_solid),intent(in)       :: fstrSOLID !< type fstr_solid
    integer(kind=kint), intent(inout) :: countNon0LU_node, countNon0LU_lagrange !< counters of node-based number of non-zero items

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
  subroutine construct_fstrMAT(np,ndof,num_lagrange,numNon0_lagrange,is_contact_active,list_nodeRelated,fstrMAT)
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

    if(associated(fstrMAT%indexL_lagrange).and.associated(fstrMAT%indexU_lagrange)) &
      deallocate(fstrMAT%indexL_lagrange,fstrMAT%indexU_lagrange)
    if(associated(fstrMAT%itemL_lagrange).and.associated(fstrMAT%itemU_lagrange)) &
      deallocate(fstrMAT%itemL_lagrange,fstrMAT%itemU_lagrange)
    if( is_contact_active ) then
      allocate(fstrMAT%indexL_lagrange(0:num_lagrange), fstrMAT%indexU_lagrange(0:np), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, fstrMAT%indexL_lagrange-fstrMAT%indexU_lagrange "
      fstrMAT%indexL_lagrange = 0 ; fstrMAT%indexU_lagrange = 0
      allocate(fstrMAT%itemL_lagrange(numNon0_lagrange), fstrMAT%itemU_lagrange(numNon0_lagrange), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, fstrMAT%itemL_lagrange-fstrMAT%itemU_lagrange "
      fstrMAT%itemL_lagrange = 0 ; fstrMAT%itemU_lagrange = 0
    endif

    if( is_contact_active ) then
      ! upper
      countNon0U_lagrange = 0
      do i = 1, np
        numI_lagrange = list_nodeRelated(i)%num_lagrange
        do j = 1, numI_lagrange
          countNon0U_lagrange = countNon0U_lagrange + 1
          fstrMAT%itemU_lagrange(countNon0U_lagrange) = list_nodeRelated(i)%id_lagrange(j)
        enddo
        fstrMAT%indexU_lagrange(i) = countNon0U_lagrange
      end do

      ! lower
      countNon0L_lagrange = 0
      do i = 1, num_lagrange
        numI_lagrange = list_nodeRelated(np+i)%num_lagrange
        do j = 1, numI_lagrange
          countNon0L_lagrange = countNon0L_lagrange + 1
          fstrMAT%itemL_lagrange(countNon0L_lagrange) = list_nodeRelated(np+i)%id_lagrange(j)
        enddo
        fstrMAT%indexL_lagrange(i) = countNon0L_lagrange
      enddo
    endif

    if(associated(fstrMAT%AL_lagrange)) deallocate(fstrMAT%AL_lagrange)
    if(associated(fstrMAT%AU_lagrange)) deallocate(fstrMAT%AU_lagrange)
    if(associated(fstrMAT%Lagrange)) deallocate(fstrMAT%Lagrange)

    if( is_contact_active ) then
      allocate(fstrMAT%AL_lagrange(ndof*numNon0_lagrange), stat=ierr)
      if ( ierr /= 0 ) stop " Allocation error, fstrMAT%AL_lagrange "
      fstrMAT%AL_lagrange = 0.0D0
      allocate(fstrMAT%AU_lagrange(ndof*numNon0_lagrange), stat=ierr)
      if ( ierr /= 0 ) stop " Allocation error, fstrMAT%AU_lagrange "
      fstrMAT%AU_lagrange = 0.0D0
      allocate(fstrMAT%Lagrange(num_lagrange))
      fstrMAT%Lagrange = 0.0D0
    endif

  end subroutine

  !> Construct new stiffness matrix structure
  subroutine constructNewMatrixStructure(hecMAT,fstrMAT,num_lagrange,numNon0_node,numNon0_lagrange,conMAT,is_contact_active)

    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstrST_matrix_contact_lagrange) :: fstrMAT !< type fstrST_matrix_contact_lagrange

    integer(kind=kint), intent(in)       :: num_lagrange
    integer(kind=kint)                   :: numNon0_node, numNon0_lagrange !< node-based number of non-zero items in half of the matrix
    integer(kind=kint)                   :: countNon0L_node, countNon0U_node !< counters of node-based number ofnon-zero items
    integer(kind=kint)                   :: i, j, ierr
    integer(kind=kint)                   :: numI_node
    integer(kind=kint)                   :: ndof, nn
    type(hecmwST_matrix)                 :: conMAT
    logical, intent(in)                  :: is_contact_active

    ndof = hecMAT%NDOF
    nn = ndof*ndof

    conMAT%N  = hecMAT%N
    conMAT%NP = hecMAT%NP
    conMAT%ndof = hecMAT%ndof
    if(associated(conMAT%indexL).and.associated(conMAT%indexU))deallocate(conMAT%indexL,conMAT%indexU)
    allocate(conMAT%indexL(0:conMAT%NP), conMAT%indexU(0:conMAT%NP), stat=ierr)
    if ( ierr /= 0) stop " Allocation error, conMAT%indexL-conMAT%indexU "
    conMAT%indexL = 0 ; conMAT%indexU = 0
    if(associated(conMAT%itemL).and.associated(conMAT%itemU))deallocate(conMAT%itemL,conMAT%itemU)
    allocate(conMAT%itemL(numNon0_node), conMAT%itemU(numNon0_node), stat=ierr)
    if ( ierr /= 0) stop " Allocation error, conMAT%itemL-conMAT%itemU "
    conMAT%itemL = 0 ; conMAT%itemU = 0
    !
    conMAT%NPL = numNon0_node
    conMAT%NPU = numNon0_node

    if(associated(hecMAT%indexL).and.associated(hecMAT%indexU))deallocate(hecMAT%indexL,hecMAT%indexU)
    allocate(hecMAT%indexL(0:hecMAT%NP), hecMAT%indexU(0:hecMAT%NP), stat=ierr)
    if ( ierr /= 0) stop " Allocation error, hecMAT%indexL-hecMAT%indexU "
    hecMAT%indexL = 0 ; hecMAT%indexU = 0
    if(associated(hecMAT%itemL).and.associated(hecMAT%itemU))deallocate(hecMAT%itemL,hecMAT%itemU)
    allocate(hecMAT%itemL(numNon0_node), hecMAT%itemU(numNon0_node), stat=ierr)
    if ( ierr /= 0) stop " Allocation error, hecMAT%itemL-hecMAT%itemU "
    hecMAT%itemL = 0 ; hecMAT%itemU = 0

    call construct_fstrMAT(hecMAT%NP,ndof,num_lagrange,numNon0_lagrange,is_contact_active,list_nodeRelated,fstrMAT)

    hecMAT%NPL = numNon0_node
    hecMAT%NPU = numNon0_node


    countNon0L_node = 0
    countNon0U_node = 0
    do i = 1, hecMAT%NP
      list_nodeRelated(i)%num_node = count(list_nodeRelated(i)%id_node /= 0)
      numI_node = list_nodeRelated(i)%num_node

      do j = 1, numI_node
        if( list_nodeRelated(i)%id_node(j) < i ) then
          countNon0L_node = countNon0L_node + 1
          hecMAT%itemL(countNon0L_node) = list_nodeRelated(i)%id_node(j)
        elseif( list_nodeRelated(i)%id_node(j) > i ) then
          countNon0U_node = countNon0U_node + 1
          hecMAT%itemU(countNon0U_node) = list_nodeRelated(i)%id_node(j)
        endif
      enddo
      hecMAT%indexL(i) = countNon0L_node
      hecMAT%indexU(i) = countNon0U_node
    end do

    do i = 1, hecMAT%NP
      deallocate(list_nodeRelated(i)%id_node)
      if(associated(list_nodeRelated(i)%id_lagrange)) deallocate(list_nodeRelated(i)%id_lagrange)
    end do

    conMAT%itemL(:)   = hecMAT%itemL(:)
    conMAT%indexL(:)  = hecMAT%indexL(:)
    conMAT%itemU(:)   = hecMAT%itemU(:)
    conMAT%indexU(:)  = hecMAT%indexU(:)

    if( is_contact_active ) then
      do i = 1, fstrMAT%num_lagrange
        deallocate(list_nodeRelated(hecMAT%NP+i)%id_lagrange)
      enddo
    endif

    deallocate(list_nodeRelated)

    if(associated(hecMAT%AL)) deallocate(hecMAT%AL)
    allocate(hecMAT%AL(nn*hecMAT%NPL), stat=ierr)
    if ( ierr /= 0 ) stop " Allocation error, hecMAT%AL "
    hecMAT%AL = 0.0D0

    if(associated(hecMAT%AU)) deallocate(hecMAT%AU)
    allocate(hecMAT%AU(nn*hecMAT%NPU), stat=ierr)
    if ( ierr /= 0 ) stop " Allocation error, hecMAT%AU "
    hecMAT%AU = 0.0D0

    if(associated(hecMAT%B)) deallocate(hecMAT%B)
    allocate(hecMAT%B(hecMAT%NP*ndof+fstrMAT%num_lagrange))
    hecMAT%B = 0.0D0

    if(associated(hecMAT%X)) deallocate(hecMAT%X)
    allocate(hecMAT%X(hecMAT%NP*ndof+fstrMAT%num_lagrange))
    hecMAT%X = 0.0D0

    if(associated(hecMAT%D)) deallocate(hecMAT%D)
    allocate(hecMAT%D(hecMAT%NP*ndof**2+fstrMAT%num_lagrange))
    hecMAT%D = 0.0D0


    if(associated(conMAT%AL)) deallocate(conMAT%AL)
    allocate(conMAT%AL(nn*conMAT%NPL), stat=ierr)
    if ( ierr /= 0 ) stop " Allocation error, conMAT%AL "
    conMAT%AL = 0.0D0

    if(associated(conMAT%AU)) deallocate(conMAT%AU)
    allocate(conMAT%AU(nn*conMAT%NPU), stat=ierr)
    if ( ierr /= 0 ) stop " Allocation error, conMAT%AU "
    conMAT%AU = 0.0D0

    if(associated(conMAT%B)) deallocate(conMAT%B)
    allocate(conMAT%B(conMAT%NP*ndof+fstrMAT%num_lagrange))
    conMAT%B = 0.0D0

    if(associated(conMAT%X)) deallocate(conMAT%X)
    allocate(conMAT%X(conMAT%NP*ndof+fstrMAT%num_lagrange))
    conMAT%X = 0.0D0

    if(associated(conMAT%D)) deallocate(conMAT%D)
    allocate(conMAT%D(conMAT%NP*ndof**2+fstrMAT%num_lagrange))
    conMAT%D = 0.0D0

  end subroutine ConstructNewMatrixStructure

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
