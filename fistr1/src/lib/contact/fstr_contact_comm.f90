!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module m_fstr_contact_comm

  use hecmw

  implicit none

  private
  public :: fstrST_contact_comm
  public :: fstr_contact_comm_init
  public :: fstr_contact_comm_finalize
  public :: fstr_contact_comm_reduce_r
  public :: fstr_contact_comm_bcast_r
  public :: fstr_contact_comm_reduce_i
  public :: fstr_contact_comm_bcast_i
  public :: fstr_contact_comm_allreduce_r
  public :: fstr_contact_comm_allreduce_i

  type fstrST_contact_comm
    private
    integer(kind=kint) :: n_neighbor_pe
    integer(kind=kint), pointer :: neighbor_pe(:)
    integer(kind=kint) :: MPI_COMM
    integer(kind=kint), pointer :: ext_index(:)
    integer(kind=kint), pointer :: ext_item(:)
    integer(kind=kint), pointer :: int_index(:)
    integer(kind=kint), pointer :: int_item(:)
  end type fstrST_contact_comm

  integer(kind=kint), parameter :: op_overwrite = 46810

  integer(kind=kint), parameter :: DEBUG = 0

contains

  subroutine fstr_contact_comm_init(conComm, hecMESH, ndof, n_contact_dof, contact_dofs)
    implicit none
    type (fstrST_contact_comm), intent(inout) :: conComm
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof, n_contact_dof
    integer(kind=kint), intent(in) :: contact_dofs(:)
    integer(kind=kint), pointer :: ext_index(:), ext_item(:), int_index(:), int_item(:)
    integer(kind=kint), allocatable :: n_ext_per_dom(:), n_int_per_dom(:), ext_item_remote(:)
    integer(kind=kint), allocatable :: statuses(:,:), requests(:)
    integer(kind=kint) :: nn_int, np, ilag, icontact, irow, irank, idom, tag, idof, idx, irow_remote
    integer(kind=kint) :: n_send, is, ie, len
    if (hecMESH%n_neighbor_pe == 0) then
      conComm%n_neighbor_pe = 0
      return
    endif
    nn_int = hecMESH%nn_internal
    np = hecMESH%n_node
    ! count external contact dofs
    allocate(n_ext_per_dom(hecMESH%n_neighbor_pe))
    n_ext_per_dom(:) = 0
    do ilag = 1, n_contact_dof
      icontact = contact_dofs(ilag)
      if (icontact <= nn_int*ndof) cycle   ! skip internal contact dof
      irow = (icontact+ndof-1) / ndof
      irank = hecMESH%node_ID(2*irow)
      call rank_to_idom(hecMESH, irank, idom)
      n_ext_per_dom(idom) = n_ext_per_dom(idom) + 1
    enddo
    ! send external / recv internal contact dofs
    allocate(statuses(HECMW_STATUS_SIZE, hecMESH%n_neighbor_pe))
    allocate(requests(hecMESH%n_neighbor_pe))
    do idom = 1, hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      tag = 5001
      call HECMW_ISEND_INT(n_ext_per_dom(idom), 1, irank, tag, &
           hecMESH%MPI_COMM, requests(idom))
    enddo
    allocate(n_int_per_dom(hecMESH%n_neighbor_pe))
    do idom = 1, hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      tag = 5001
      call HECMW_RECV_INT(n_int_per_dom(idom), 1, irank, tag, &
           hecMESH%MPI_COMM, statuses(:,1))
    enddo
    call HECMW_Waitall(hecMESH%n_neighbor_pe, requests, statuses)
    ! make index
    allocate(ext_index(0:hecMESH%n_neighbor_pe))
    allocate(int_index(0:hecMESH%n_neighbor_pe))
    ext_index(0) = 0
    int_index(0) = 0
    do idom = 1, hecMESH%n_neighbor_pe
      ext_index(idom) = ext_index(idom-1) + n_ext_per_dom(idom)
      int_index(idom) = int_index(idom-1) + n_int_per_dom(idom)
    enddo
    ! make ext_item
    allocate(ext_item(ext_index(hecMESH%n_neighbor_pe)))
    allocate(ext_item_remote(ext_index(hecMESH%n_neighbor_pe)))
    n_ext_per_dom(:) = 0
    do ilag = 1, n_contact_dof
      icontact = contact_dofs(ilag)
      if (icontact <= nn_int*ndof) cycle   ! skip internal contact dof
      irow = (icontact+ndof-1) / ndof
      idof = icontact - ndof*(irow-1)
      irank = hecMESH%node_ID(2*irow)
      call rank_to_idom(hecMESH, irank, idom)
      n_ext_per_dom(idom) = n_ext_per_dom(idom) + 1
      idx = ext_index(idom-1)+n_ext_per_dom(idom)
      ext_item(idx) = icontact
      irow_remote = hecMESH%node_ID(2*irow-1)
      ext_item_remote(idx) = ndof*(irow_remote-1)+idof
    enddo
    deallocate(n_ext_per_dom)
    deallocate(n_int_per_dom)
    ! send ext_item_remote and recv int_item
    n_send = 0
    do idom = 1, hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      is = ext_index(idom-1)+1
      ie = ext_index(idom)
      len = ie-is+1
      if (len == 0) cycle
      n_send = n_send + 1
      tag = 5002
      call HECMW_ISEND_INT(ext_item_remote(is:ie), len, irank, tag, &
           hecMESH%MPI_COMM, requests(n_send))
    enddo
    allocate(int_item(int_index(hecMESH%n_neighbor_pe)))
    do idom = 1, hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      is = int_index(idom-1)+1
      ie = int_index(idom)
      len = ie-is+1
      if (len == 0) cycle
      tag = 5002
      call HECMW_RECV_INT(int_item(is:ie), len, irank, tag, &
           hecMESH%MPI_COMM, statuses(:,1))
    enddo
    call HECMW_Waitall(n_send, requests, statuses)
    deallocate(statuses, requests)
    if (DEBUG >= 2) then
      write(0,*) '  DEBUG2: ext_index',ext_index(:)
      write(0,*) '  DEBUG2: ext_item',ext_item(:)
      write(0,*) '  DEBUG2: ext_item_remote',ext_item_remote(:)
      write(0,*) '  DEBUG2: int_index',int_index(:)
      write(0,*) '  DEBUG2: int_item',int_item(:)
    endif
    deallocate(ext_item_remote)
    !
    conComm%n_neighbor_pe = hecMESH%n_neighbor_pe
    allocate(conComm%neighbor_pe(conComm%n_neighbor_pe))
    conComm%neighbor_pe(:) = hecMESH%neighbor_pe(:)
    conComm%MPI_COMM = hecMESH%MPI_COMM
    conComm%ext_index => ext_index
    conComm%ext_item => ext_item
    conComm%int_index => int_index
    conComm%int_item => int_item
  end subroutine fstr_contact_comm_init

  subroutine fstr_contact_comm_finalize(conComm)
    implicit none
    type (fstrST_contact_comm), intent(inout) :: conComm
    if (conComm%n_neighbor_pe == 0) return
    if (associated(conComm%neighbor_pe)) deallocate(conComm%neighbor_pe)
    if (associated(conComm%ext_index)) deallocate(conComm%ext_index)
    if (associated(conComm%ext_item)) deallocate(conComm%ext_item)
    if (associated(conComm%int_index)) deallocate(conComm%int_index)
    if (associated(conComm%int_item)) deallocate(conComm%int_item)
    conComm%n_neighbor_pe = 0
    conComm%MPI_COMM = 0
  end subroutine fstr_contact_comm_finalize

  subroutine fstr_contact_comm_reduce_r(conComm, vec, op)
    implicit none
    type (fstrST_contact_comm), intent(in) :: conComm
    real(kind=kreal), intent(inout) :: vec(:)
    integer(kind=kint), intent(in) :: op
    if (conComm%n_neighbor_pe == 0) return
    call send_recv_contact_info_r(conComm%n_neighbor_pe, conComm%neighbor_pe, conComm%MPI_COMM, &
         conComm%ext_index, conComm%ext_item, conComm%int_index, conComm%int_item, vec, op)
  end subroutine fstr_contact_comm_reduce_r

  subroutine fstr_contact_comm_bcast_r(conComm, vec)
    implicit none
    type (fstrST_contact_comm), intent(in) :: conComm
    real(kind=kreal), intent(inout) :: vec(:)
    integer(kind=kint) :: op
    if (conComm%n_neighbor_pe == 0) return
    op = op_overwrite
    call send_recv_contact_info_r(conComm%n_neighbor_pe, conComm%neighbor_pe, conComm%MPI_COMM, &
         conComm%int_index, conComm%int_item, conComm%ext_index, conComm%ext_item, vec, op)
  end subroutine fstr_contact_comm_bcast_r

  subroutine fstr_contact_comm_reduce_i(conComm, vec, op)
    implicit none
    type (fstrST_contact_comm), intent(in) :: conComm
    integer(kind=kint), intent(inout) :: vec(:)
    integer(kind=kint), intent(in) :: op
    if (conComm%n_neighbor_pe == 0) return
    call send_recv_contact_info_i(conComm%n_neighbor_pe, conComm%neighbor_pe, conComm%MPI_COMM, &
         conComm%ext_index, conComm%ext_item, conComm%int_index, conComm%int_item, vec, op)
  end subroutine fstr_contact_comm_reduce_i

  subroutine fstr_contact_comm_bcast_i(conComm, vec)
    implicit none
    type (fstrST_contact_comm), intent(in) :: conComm
    integer(kind=kint), intent(inout) :: vec(:)
    integer(kind=kint) :: op
    if (conComm%n_neighbor_pe == 0) return
    op = op_overwrite
    call send_recv_contact_info_i(conComm%n_neighbor_pe, conComm%neighbor_pe, conComm%MPI_COMM, &
         conComm%int_index, conComm%int_item, conComm%ext_index, conComm%ext_item, vec, op)
  end subroutine fstr_contact_comm_bcast_i

  subroutine fstr_contact_comm_allreduce_r(conComm, vec, op)
    implicit none
    type (fstrST_contact_comm), intent(in) :: conComm
    real(kind=kreal), intent(inout) :: vec(:)
    integer(kind=kint), intent(in) :: op
    call fstr_contact_comm_reduce_r(conComm, vec, op)
    call fstr_contact_comm_bcast_r(conComm, vec)
  end subroutine fstr_contact_comm_allreduce_r

  subroutine fstr_contact_comm_allreduce_i(conComm, vec, op)
    implicit none
    type (fstrST_contact_comm), intent(in) :: conComm
    integer(kind=kint), intent(inout) :: vec(:)
    integer(kind=kint), intent(in) :: op
    call fstr_contact_comm_reduce_i(conComm, vec, op)
    call fstr_contact_comm_bcast_i(conComm, vec)
  end subroutine fstr_contact_comm_allreduce_i

  !
  ! private subroutines
  !

  subroutine rank_to_idom(hecMESH, rank, idom)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: rank
    integer(kind=kint), intent(out) :: idom
    integer(kind=kint) :: i
    do i = 1, hecMESH%n_neighbor_pe
      if (hecMESH%neighbor_pe(i) == rank) then
        idom = i
        return
      endif
    enddo
    stop 'ERROR: exp_rank not found in neighbor_pe'
  end subroutine rank_to_idom

  subroutine send_recv_contact_info_r(n_neighbor_pe, neighbor_pe, MPI_COMM, &
       send_index, send_item, recv_index, recv_item, vec, op)
    implicit none
    integer(kind=kint), intent(in) :: n_neighbor_pe
    integer(kind=kint), intent(in) :: neighbor_pe(:)
    integer(kind=kint), intent(in) :: MPI_COMM
    integer(kind=kint), pointer, intent(in) :: send_index(:), send_item(:), recv_index(:), recv_item(:)
    real(kind=kreal), intent(inout) :: vec(:)
    integer(kind=kint), intent(in) :: op
    real(kind=kreal), allocatable :: send_buf(:), recv_buf(:)
    integer(kind=kint) :: i, n_send, idom, irank, is, ie, len, tag
    integer(kind=kint), allocatable :: requests(:), statuses(:,:)
    if (n_neighbor_pe == 0) return
    allocate(requests(n_neighbor_pe))
    allocate(statuses(HECMW_STATUS_SIZE, n_neighbor_pe))
    allocate(send_buf(send_index(n_neighbor_pe)))
    allocate(recv_buf(recv_index(n_neighbor_pe)))
    do i = 1, send_index(n_neighbor_pe)
      send_buf(i) = vec(send_item(i))
    enddo
    n_send = 0
    do idom = 1, n_neighbor_pe
      irank = neighbor_pe(idom)
      is = send_index(idom-1)+1
      ie = send_index(idom)
      len = ie-is+1
      if (len == 0) cycle
      n_send = n_send + 1
      tag = 5011
      call HECMW_ISEND_R(send_buf(is:ie), len, irank, tag, &
           MPI_COMM, requests(n_send))
    enddo
    do idom = 1, n_neighbor_pe
      irank = neighbor_pe(idom)
      is = recv_index(idom-1)+1
      ie = recv_index(idom)
      len = ie-is+1
      if (len == 0) cycle
      tag = 5011
      call HECMW_RECV_R(recv_buf(is:ie), len, irank, tag, &
           MPI_COMM, statuses(:,1))
    enddo
    call HECMW_Waitall(n_send, requests, statuses)
    if (op == HECMW_SUM) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = vec(recv_item(i)) + recv_buf(i)
      enddo
    elseif (op == HECMW_PROD) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = vec(recv_item(i)) * recv_buf(i)
      enddo
    elseif (op == HECMW_MAX) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = max(vec(recv_item(i)), recv_buf(i))
      enddo
    elseif (op == HECMW_MIN) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = min(vec(recv_item(i)), recv_buf(i))
      enddo
    else  ! overwrite
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = recv_buf(i)
      enddo
    endif
    deallocate(requests)
    deallocate(statuses)
    if (DEBUG >= 2) then
      write(0,*) '  DEBUG2: send_buf',send_buf(:)
      write(0,*) '  DEBUG2: recv_buf',recv_buf(:)
    endif
    deallocate(send_buf)
    deallocate(recv_buf)
  end subroutine send_recv_contact_info_r

  subroutine send_recv_contact_info_i(n_neighbor_pe, neighbor_pe, MPI_COMM, &
       send_index, send_item, recv_index, recv_item, vec, op)
    implicit none
    integer(kind=kint), intent(in) :: n_neighbor_pe
    integer(kind=kint), intent(in) :: neighbor_pe(:)
    integer(kind=kint), intent(in) :: MPI_COMM
    integer(kind=kint), pointer, intent(in) :: send_index(:), send_item(:), recv_index(:), recv_item(:)
    integer(kind=kint), intent(inout) :: vec(:)
    integer(kind=kint), intent(in) :: op
    integer(kind=kint), allocatable :: send_buf(:), recv_buf(:)
    integer(kind=kint) :: i, n_send, idom, irank, is, ie, len, tag
    integer(kind=kint), allocatable :: requests(:), statuses(:,:)
    if (n_neighbor_pe == 0) return
    allocate(requests(n_neighbor_pe))
    allocate(statuses(HECMW_STATUS_SIZE, n_neighbor_pe))
    allocate(send_buf(send_index(n_neighbor_pe)))
    allocate(recv_buf(recv_index(n_neighbor_pe)))
    do i = 1, send_index(n_neighbor_pe)
      send_buf(i) = vec(send_item(i))
    enddo
    n_send = 0
    do idom = 1, n_neighbor_pe
      irank = neighbor_pe(idom)
      is = send_index(idom-1)+1
      ie = send_index(idom)
      len = ie-is+1
      if (len == 0) cycle
      n_send = n_send + 1
      tag = 5011
      call HECMW_ISEND_INT(send_buf(is:ie), len, irank, tag, &
           MPI_COMM, requests(n_send))
    enddo
    do idom = 1, n_neighbor_pe
      irank = neighbor_pe(idom)
      is = recv_index(idom-1)+1
      ie = recv_index(idom)
      len = ie-is+1
      if (len == 0) cycle
      tag = 5011
      call HECMW_RECV_INT(recv_buf(is:ie), len, irank, tag, &
           MPI_COMM, statuses(:,1))
    enddo
    call HECMW_Waitall(n_send, requests, statuses)
    if (op == HECMW_SUM) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = vec(recv_item(i)) + recv_buf(i)
      enddo
    elseif (op == HECMW_PROD) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = vec(recv_item(i)) * recv_buf(i)
      enddo
    elseif (op == HECMW_MAX) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = max(vec(recv_item(i)), recv_buf(i))
      enddo
    elseif (op == HECMW_MIN) then
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = min(vec(recv_item(i)), recv_buf(i))
      enddo
    else  ! overwrite
      do i = 1, recv_index(n_neighbor_pe)
        vec(recv_item(i)) = recv_buf(i)
      enddo
    endif
    deallocate(requests)
    deallocate(statuses)
    if (DEBUG >= 2) then
      write(0,*) '  DEBUG2: send_buf',send_buf(:)
      write(0,*) '  DEBUG2: recv_buf',recv_buf(:)
    endif
    deallocate(send_buf)
    deallocate(recv_buf)
  end subroutine send_recv_contact_info_i

end module m_fstr_contact_comm
